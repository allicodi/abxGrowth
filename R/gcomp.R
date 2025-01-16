#'
#' Function to do g-computation for antibiotics growth analysis
#' 
#' @param data dataframe containing dataset to use for gcomp
#' @param laz_var_name name of growth outcome variable
#' @param abx_var_name name of binary antibiotic variable
#' @param infection_var_name name of binary infection variable
#' @param severity_list character vector containing names of severity-related covariates (post-infection). If NULL, perform traditional gcomputation. Else, perform longitudinal gcomputation
#' @param covariate_list character vector containing names of baseline covariates
#' @param site_var_name name of site covariate in dataset (if applicable, else null)
#' @param site_interaction TRUE or FALSE indicating interaction between site and antibiotics 
#' @param age_var_name name of age covariate in dataset (if applicable, else NULL)
#' @param outcome_type gaussian or binomial for continuous or binomial outcome
#' 
#' @import splines
#' @export
#' 
#' @returns 
#' \describe{
#'  List of containing the following:
#'  \item{\code{effect_inf_no_abx}}{numeric effect of infection on growth in subgroup that did *not* receive antibiotics}
#'  \item{\code{effect_inf_abx}}{numeric effect of infection on growth in subgroup that received antibiotics}
#'  \item{\code{abx_0_inf_1}}{expected growth outcome in infected subgroup who did not receive abx}
#'  \item{\code{abx_0_inf_0}}{expected growth outcome in uninfected subgroup who did not receive abx}
#'  \item{\code{abx_1_inf_1}}{expected growth outcome in infected subgroup who received abx}
#'  \item{\code{abx_1_inf_0}}{expected growth outcome in uninfected subgroup who receieved abx}
#'  }
abx_growth_gcomp <- function(data, 
                             laz_var_name = "mo3_haz",
                             abx_var_name = "who_rec_abx",
                             infection_var_name = "tac_shigella_attributable", # causal exposure
                             severity_list = c(
                               "enroll_diar_blood",
                               "enroll_diar_vom_days",
                               "enroll_diar_fever_days",
                               "enroll_diar_vom_num"
                             ),
                             covariate_list = c( #separate covariate lists for exposed and unexposed 
                               "sex",
                               "enr_age_months",
                               "enr_haz",
                               "final_quintile",
                               "enroll_site"
                             ),
                             site_var_name = "enroll_site",
                             site_interaction = TRUE,
                             age_var_name = "enr_age_months",
                             outcome_type = "gaussian"){
  
  # if age in covariate_list, change to spline with 3 knots to increase flexibility
  if(age_var_name %in% covariate_list){
    covariate_list <- covariate_list[covariate_list != age_var_name]
    covariate_list <- c(covariate_list, paste0("splines::ns(", age_var_name, ", df = 4)"))
  }
  
  if(!is.null(severity_list)){
    
    # -----------------------------------
    # Longitudinal G-Computation 
    # -----------------------------------
    
    # Get model formula 
    if(site_interaction == "TRUE"){
      
      # if site in covariate list, remove for model1 (only want interaction terms)
      if(length(site_var_name) == 1){
        if(site_var_name %in% covariate_list){
          covariate_list <- covariate_list[covariate_list != site_var_name]
        }
      } else if (any(site_var_name %in% covariate_list)){
        covariate_list <- covariate_list[!covariate_list %in% site_var_name]
      }
      
      model_1_formula <- stats::as.formula(paste(laz_var_name, "~", 
                                          abx_var_name, "+", infection_var_name, "+",
                                          paste(covariate_list, collapse = "+"), "+", # union of both covariate lists -- collapse 
                                          paste(severity_list, collapse = "+"), "+",
                                          paste(site_var_name, "*", abx_var_name, collapse = "+"), "+",
                                          paste0(infection_var_name, "*", abx_var_name))) 
      
      # put site back in covariate list
      if(!is.null(site_var_name)){
        covariate_list <- c(covariate_list, site_var_name)
      }
      
    } else {
      
      # if site was listed but is not already in covariate list, add
      if(!is.null(site_var_name) & (!(site_var_name %in% covariate_list))){
        covariate_list <- c(covariate_list, site_var_name)
      }
      
      model_1_formula <- stats::as.formula(paste(laz_var_name, "~", 
                                          abx_var_name, "+", infection_var_name, "+", 
                                          paste(covariate_list, collapse = "+"), "+",
                                          paste(severity_list, collapse = "+"), "+",
                                          paste0(infection_var_name, "*", abx_var_name)))
    }
    
    # Step 1: regress LAZ on abx, infection, all severity / non-mediating variables
    model_1 <- stats::glm(model_1_formula,
                   data = data, 
                   family = outcome_type)
    
    # Step 2: predict from model setting abx = 0, infection = 1
    data_01 <- data
    data_01[[infection_var_name]] <- 1
    data_01[[abx_var_name]] <- 0
    
    # predict from the subset model (that no longer includes abx or infection)
    yhat_01 <- stats::predict(model_1, newdata = data_01)
    
    # Step 3: regress yhat_01 on all other non-mediating variables in subset with infection = 1, call this model2
    data$yhat_01 <- yhat_01
    sub_inf_1 <- data[data[[infection_var_name]] == 1,]
    
    # covariate list for the exposed (may or may not include the pathogens)
    model2 <- stats::glm(stats::as.formula(paste("yhat_01", "~", paste(covariate_list, collapse = "+"))),
                  data = sub_inf_1,
                  family = outcome_type)
    
    ## NEW: GLMM style avg over child first
    
    # Step 4a: predict from model2 on everyone
    data$ybar_01_preds <- stats::predict(model2, newdata = data)
    
    # Step 4b: average over child (if child_id in data) then avg over all children
    if("child_id" %in% colnames(data)){
      # multiple enrollments avg together
      ybar_01_child <- aggregate(ybar_01_preds ~ child_id, data = data, FUN = mean, na.rm=TRUE)$ybar_01_preds
      ybar_01 <- mean(ybar_01_child, na.rm = TRUE)
    } else {
      # assume all children unique, skip to mean of preds directly
      ybar_01 <- mean(data$ybar_01_preds, na.rm = TRUE)
    }
    
    # Old:
    # ybar_01 <- mean(stats::predict(model2, newdata = data), na.rm = TRUE)
    
    # Step 5: predict from model1 setting infection = 0 and abx = 0, call that yhat_00
    data_00 <- data
    data_00[[infection_var_name]] <- 0
    data_00[[abx_var_name]] <- 0
    
    yhat_00 <- stats::predict(model_1, newdata = data_00)
    
    # Step 6: regress yhat_00 on all other non-mediating variables in subset with infection = 0, call this model3
    data$yhat_00 <- yhat_00
    sub_inf_0 <- data[data[[infection_var_name]] == 0,]
    
    # include quant shigella in model for uninfected
    model3 <- stats::glm(stats::as.formula(paste("yhat_00", "~", paste(covariate_list, collapse = "+"))),
                  data = sub_inf_0,
                  family = outcome_type)
    
    ## NEW: GLMM style avg over child first
    
    # Step 7a: stats::predict from model3 on everyone
    data$ybar_00_preds <- stats::predict(model3, newdata = data)
    
    # Step 7b: average over child (if child_id in data) then avg over all children
    if("child_id" %in% colnames(data)){
      # multiple enrollments avg together
      ybar_00_child <- aggregate(ybar_00_preds ~ child_id, data = data, FUN = mean, na.rm = TRUE)$ybar_00_preds
      ybar_00 <- mean(ybar_00_child, na.rm = TRUE)
    } else {
      # assume all children unique, skip to mean of preds directly
      ybar_00 <- mean(data$ybar_00_preds, na.rm = TRUE)
    }
    
    # Old:
    #ybar_00 <- mean(stats::predict(model3, newdata = data), na.rm = TRUE)
    
    # Step 8: effect of infection without abx = ybar_01 - ybar_00
    inf_no_abx <- ybar_01 - ybar_00
    
    # Step 9: stats::predict from model1 setting infection = 1 and abx = 1, call that yhat_11
    data_11 <- data
    data_11[[infection_var_name]] <- 1
    data_11[[abx_var_name]] <- 1
    
    yhat_11 <- stats::predict(model_1, newdata = data_11)
    
    # Step 10: regress yhat_11 on all other non-mediating variables in subset with infection = 1, call this model4
    data$yhat_11 <- yhat_11
    sub_inf_1 <- data[data[[infection_var_name]] == 1,]
    
    # same as above -- don't include quant shig here
    model4 <- stats::glm(stats::as.formula(paste("yhat_11", "~", paste(covariate_list, collapse = "+"))),
                  data = sub_inf_1,
                  family = outcome_type)
    
    
    
    ## NEW: GLMM style avg over child first
    
    # Step 11a: predict from model4 on everyone
    data$ybar_11_preds <- stats::predict(model4, newdata = data)
    
    # Step 11b: average over child (if child_id in data) then avg over all children
    if("child_id" %in% colnames(data)){
      # multiple enrollments avg together
      ybar_11_child <- aggregate(ybar_11_preds ~ child_id, data = data, FUN = mean, na.rm=TRUE)$ybar_11_preds
      ybar_11 <- mean(ybar_11_child, na.rm = TRUE)
    } else {
      # assume all children unique, skip to mean of preds directly
      ybar_11 <- mean(data$ybar_11_preds, na.rm = TRUE)
    }
    
    # Old:
    # ybar_11 <- mean(stats::predict(model4, newdata = data), na.rm = TRUE)
    
    # Step 12: stats::predict from model1 setting infection = 0 and abx = 1, call that yhat_10
    data_10 <- data
    
    data_10[[infection_var_name]] <- 0
    data_10[[abx_var_name]] <- 1 
    
    yhat_10 <- stats::predict(model_1, newdata = data_10)
    
    # Step 13: regress yhat_10 on all other blah blah in subset with infection = 0, call this model5
    data$yhat_10 <- yhat_10
    sub_inf_0 <- data[data[[infection_var_name]] == 0,]
    
    # same as above -- include quant shig here
    model5 <- stats::glm(stats::as.formula(paste("yhat_10", "~", paste(covariate_list, collapse = "+"))),
                  data = sub_inf_0,
                  family = outcome_type)
    
    # Step 14a: predict from model4 on everyone
    data$ybar_10_preds <- stats::predict(model5, newdata = data)
    
    # Step 14b: average over child (if child_id in data) then avg over all children
    if("child_id" %in% colnames(data)){
      # multiple enrollments avg together
      ybar_10_child <- aggregate(ybar_10_preds ~ child_id, data = data, FUN = mean, na.rm=TRUE)$ybar_10_preds
      ybar_10 <- mean(ybar_10_child, na.rm = TRUE)
    } else {
      # assume all children unique, skip to mean of preds directly
      ybar_10 <- mean(data$ybar_10_preds, na.rm = TRUE)
    }
    
    # Old:
    # ybar_10 <- mean(stats::predict(model5, newdata = data), na.rm = TRUE)
    
    # Step 15: effect of infection with abx = ybar_11 - ybar_10
    inf_abx <- ybar_11 - ybar_10
    
    # RETURN RESULTS
    return(list(effect_inf_no_abx = inf_no_abx,
                effect_inf_abx = inf_abx,
                abx_0_inf_1 = ybar_01,
                abx_0_inf_0 = ybar_00,
                abx_1_inf_1 = ybar_11,
                abx_1_inf_0 = ybar_10))
  } else {
    
    # ----------------------------------------
    # Traditional G-Computation 
    # ----------------------------------------
    if(site_interaction == "TRUE"){
      model_1_formula <- stats::as.formula(paste(laz_var_name, "~",
                                          abx_var_name, "+", infection_var_name, "+",
                                          paste(covariate_list, collapse = "+"), "+",
                                          paste0(site_var_name, "*", abx_var_name), "+",
                                          paste0(infection_var_name, "*", abx_var_name)))
    } else {
      model_1_formula <- stats::as.formula(paste(laz_var_name, "~",
                                          abx_var_name, "+", infection_var_name, "+",
                                          paste(covariate_list, collapse = "+"), "+",
                                          paste0(infection_var_name, "*", abx_var_name)))
    }
    
    # Fit model with LAZ ~ infection + abx + covariates + infection*abx
    model_1 <- stats::glm(model_1_formula,
                   data = data,
                   family = outcome_type)
    
    # Predict Abx 0, infection 1
    data_01 <- data
    data_01[[infection_var_name]] <- 1
    data_01[[abx_var_name]] <- 0
    
    ## NEW: GLMM style avg over child first
    
    data$yhat_01 <- stats::predict(model_1, newdata = data_01)
    if("child_id" %in% colnames(data)){
      yhat_01_child <- aggregate(yhat_01 ~ child_id, data = data, FUN = mean, na.rm=TRUE)$yhat_01
      ybar_01 <- mean(yhat_01_child, na.rm = TRUE)
    } else {
      ybar_01 <- mean(data$yhat_01, na.rm = TRUE)
    }
    
    #ybar_01 <- mean(stats::predict(model_1, newdata = data_01), na.rm = TRUE)
    
    # Predict Abx 0, infection 0
    data_00 <- data
    data_00[[infection_var_name]] <- 0
    data_00[[abx_var_name]] <- 0
    
    ## NEW: GLMM style avg over child first
    
    data$yhat_00 <- stats::predict(model_1, newdata = data_00)
    if("child_id" %in% colnames(data)){
      yhat_00_child <- aggregate(yhat_00 ~ child_id, data = data, FUN = mean, na.rm=TRUE)$yhat_00
      ybar_00 <- mean(yhat_00_child, na.rm = TRUE)
    } else {
      ybar_01 <- mean(data$yhat_01, na.rm = TRUE)
    }
    
    #ybar_00 <- mean(stats::predict(model_1, newdata = data_00), na.rm = TRUE)
    
    # Difference
    inf_no_abx <- ybar_01 - ybar_00
    
    # Predict Abx 1, infection 1
    data_11 <- data
    data_11[[infection_var_name]] <- 1
    data_11[[abx_var_name]] <- 1
    
    
    ## NEW: GLMM style avg over child first
    
    data$yhat_11 <- stats::predict(model_1, newdata = data_11)
    if("child_id" %in% colnames(data)){
      yhat_11_child <- aggregate(yhat_11 ~ child_id, data = data, FUN = mean, na.rm=TRUE)$yhat_11
      ybar_11 <- mean(yhat_11_child, na.rm = TRUE)
    } else {
      ybar_11 <- mean(data$yhat_11, na.rm = TRUE)
    }
    
    #ybar_11 <- mean(stats::predict(model_1, newdata = data_11), na.rm = TRUE)
    
    # Predict Abx 1, infection 0
    data_10 <- data
    data_10[[infection_var_name]] <- 0
    data_10[[abx_var_name]] <- 1
    
    ## NEW: GLMM style avg over child first
    
    data$yhat_10 <- stats::predict(model_1, newdata = data_10)
    if("child_id" %in% colnames(data)){
      yhat_10_child <- aggregate(yhat_10 ~ child_id, data = data, FUN = mean, na.rm=TRUE)$yhat_10
      ybar_10 <- mean(yhat_10_child, na.rm = TRUE)
    } else {
      ybar_10 <- mean(data$yhat_10, na.rm = TRUE)
    }
    
    #ybar_10 <- mean(stats::predict(model_1, newdata = data_10), na.rm = TRUE)
    
    # Difference
    inf_abx <- ybar_11 - ybar_10
    
    # return list with results
    return(list(effect_inf_no_abx = inf_no_abx,
                effect_inf_abx = inf_abx, 
                abx_0_inf_1 = ybar_01,
                abx_0_inf_0 = ybar_00,
                abx_1_inf_1 = ybar_11,
                abx_1_inf_0 = ybar_10))
  }
  
}


#'
#' Function to do g-computation for antibiotics growth analysis for case-control data
#' 
#' @param data dataframe containing dataset to use for gcomp
#' @param laz_var_name name of growth outcome variable
#' @param abx_var_name name of binary antibiotic variable
#' @param case_var_name name of binary exposure/case variable 
#' @param severity_list character vector containing names of severity-related covariates (post-infection). 
#' @param covariate_list character vector containing names of baseline covariates
#' @param covariate_list_control character vector containing names of baseline covariates for controls. If NULL, same covariate_list as cases
#' @param site_var_name name of site covariate in dataset (if applicable, else null)
#' @param site_interaction TRUE or FALSE indicating interaction between site and antibiotics 
#' @param age_var_name name of age covariate in dataset (if applicable, else NULL)
#' @param outcome_type gaussian or binomial for continuous or binomial outcome
#' 
#' @import splines
#' @export
#' 
#' @returns 
#' \describe{
#'  List of containing the following:
#'  \item{\code{effect_inf_no_abx}}{effect of infection on growth in subgroup that did *not* receive antibiotics vs heathy controls}
#'  \item{\code{effect_inf_abx}}{effect of infection on growth in subgroup that received antibiotics vs healthy controls}
#'  \item{\code{abx_0_inf_1}}{expected growth outcome in infected subgroup who did not receive abx}
#'  \item{\code{abx_0_inf_0}}{expected growth outcome in healthy controls}
#'  \item{\code{abx_1_inf_1}}{expected growth outcome in infected subgroup who received abx}
#'  }
abx_growth_gcomp_case_control <- function(data, 
                             laz_var_name = "mo3_haz",
                             abx_var_name = "who_rec_abx",
                             case_var_name = "case", # binary variable = 1 if shigella attributable diarrhea / 0 if enrolled as control
                             severity_list = c(
                               "enroll_diar_blood",
                               "enroll_diar_vom_days",
                               "enroll_diar_fever_days",
                               "enroll_diar_vom_num"
                             ),
                             covariate_list = c(
                               "sex",
                               "enr_age_months",
                               "enr_haz",
                               "final_quintile",
                               "enroll_site"
                             ),
                             covariate_list_control = NULL,
                             site_var_name = "enroll_site",
                             site_interaction = TRUE,
                             age_var_name = "enr_age_months",
                             outcome_type = "gaussian"){
  
  # Estimands of interest:
  # E[Growth(Shigella diar = 1, Antibiotics = a) - Growth(Shigella diar = 0) | control ] = 
  #    E[ E [ Growth | Shigella diar = 1, Antibiotics = a, severity, baseline] | Shigella diar = 1, baseline ] | control ] -
  #    E[ Growth | control ]
  
  # ---------------------------------------------------------------------------------------------------------------------------
  # Part 1: E[ E [ Growth | Shigella diar = 1, Antibiotics = a, severity, baseline] | Shigella diar = 1, baseline ] | control ]
  # ---------------------------------------------------------------------------------------------------------------------------
  
  # get case vs control data
  case_data <- data[data[[case_var_name]] == 1,]
  control_data <- data[data[[case_var_name]] == 0,]
  
  # healthy controls only (no abx)
  control_data <- control_data[control_data[[abx_var_name]] == 0,]
  
  # if null, same list as cases
  if(is.null(covariate_list_control)){
    covariate_list_control <- covariate_list
  } 
  
  # if age in covariate_list, change to spline with 3 knots to increase flexibility
  if(age_var_name %in% covariate_list){
    covariate_list <- covariate_list[covariate_list != age_var_name]
    covariate_list <- c(covariate_list, paste0("splines::ns(", age_var_name, ", df = 4)"))
  }
  
  # Get model formula
  if(site_interaction == "TRUE"){

    # if site in covariate list, remove for model1 (only want interaction terms)
    if(length(site_var_name) == 1){
      if(site_var_name %in% covariate_list){
        covariate_list <- covariate_list[covariate_list != site_var_name]
      }
    } else if (any(site_var_name %in% covariate_list)){
      covariate_list <- covariate_list[!covariate_list %in% site_var_name]
    }
    
    model_1_formula_case <- stats::as.formula(
      paste(laz_var_name, "~",
       abx_var_name, "+",
       paste(covariate_list, collapse = "+"), "+",
       paste(severity_list, collapse = "+"), "+",
       paste(site_var_name, "*", abx_var_name, collapse = "+"))
    )
    
    # put site back in covariate list
    if(!is.null(site_var_name)){
      covariate_list <- c(covariate_list, site_var_name)
    }

  } else {

    # if site was listed but is not already in covariate list, add
    if(!is.null(site_var_name) & (!(site_var_name %in% covariate_list))){
      covariate_list <- c(covariate_list, site_var_name)
    }

    model_1_formula_case <- stats::as.formula(paste(laz_var_name, "~",
                                               abx_var_name, "+", case_var_name, "+",
                                               paste(covariate_list, collapse = "+"), "+",
                                               paste(severity_list, collapse = "+")))
  }

  # Regress LAZ on abx, infection, all severity / non-mediating variables
  model_1_case <- stats::glm(model_1_formula_case,
                        data = case_data,
                        family = outcome_type)

  # Predict from model setting abx = 0
  data_01 <- case_data
  data_01[[abx_var_name]] <- 0
  yhat_01 <- stats::predict(model_1_case, newdata = data_01, type = "response")

  # Regress yhat_01 on all other non-mediating variables 
  case_data$yhat_01 <- yhat_01

  model2 <- stats::glm(stats::as.formula(paste("yhat_01", "~", paste(covariate_list, collapse = "+"))),
                       data = case_data,
                       family = outcome_type)

  # Predict from model2 on controls
  ybar_01 <- mean(stats::predict(model2, newdata = control_data, type = "response"), na.rm = TRUE)
  
  # stats::predict from model1 setting infection = 1 and abx = 1, call that yhat_11
  data_11 <- case_data
  data_11[[abx_var_name]] <- 1

  yhat_11 <- stats::predict(model_1_case, newdata = data_11, type = "response")

  # Regress yhat_11 on all other non-mediating variables, call this model4
  case_data$yhat_11 <- yhat_11

  model4 <- stats::glm(stats::as.formula(paste("yhat_11", "~", paste(covariate_list, collapse = "+"))),
                       data = case_data,
                       family = outcome_type)

  # stats::predict from model4 on controls, average stats::predictions, call avg ybar_11
  ybar_11 <- mean(stats::predict(model4, newdata = control_data, type = "response"), na.rm = TRUE)

  # ------------------------------------------------------------
  # Part 2: E[Growth | Control]
  # ------------------------------------------------------------
  
  ybar_00 <- mean(control_data[[laz_var_name]], na.rm = TRUE)
  
  # -----------------------------------------------------------
  # Part 3: Take difference to get estimands of interest
  # -----------------------------------------------------------
  
  # Effect of infection with abx in reference to healthy controls
  effect_abx <- ybar_11 - ybar_00
  
  # Effect of infection without abx in reference to healthy controls
  effect_no_abx <- ybar_01 - ybar_00

  return(list(effect_inf_no_abx = effect_no_abx,
              effect_inf_abx = effect_abx,
              abx_0_inf_1 = ybar_01,
              abx_0_inf_0 = ybar_00,
              abx_1_inf_1 = ybar_11))
}

