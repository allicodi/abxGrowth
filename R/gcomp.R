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
#' @param att boolean if effect should be estimated among people who would naturally get infection, default TRUE
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
                             outcome_type = "gaussian",
                             att = TRUE){
  
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
      
      model1_formula <- stats::as.formula(paste(laz_var_name, "~", 
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
      
      model1_formula <- stats::as.formula(paste(laz_var_name, "~", 
                                          abx_var_name, "+", infection_var_name, "+", 
                                          paste(covariate_list, collapse = "+"), "+",
                                          paste(severity_list, collapse = "+"), "+",
                                          paste0(infection_var_name, "*", abx_var_name)))
    }
    
    # Step 1: regress LAZ on abx, infection, all severity / non-mediating variables
    model1 <- stats::glm(model1_formula,
                   data = data, 
                   family = outcome_type)
    
    
    # Step 2: For each level of abx, predict setting abx = x, infection = 1 & abx = x, infection = 0
    abx_levels <- unique(data[[abx_var_name]])[!is.na(unique(data[[abx_var_name]]))]
    results_df <- data.frame(abx_levels = abx_levels,
                             abx_level_inf_1 = vector("numeric", length = length(abx_levels)),
                             abx_level_inf_0 = vector("numeric", length = length(abx_levels)),
                             effect_inf_abx_level = vector("numeric", length = length(abx_levels)))
    
    for(i in 1:length(abx_levels)){
      
      abx_level <- abx_levels[i]
      
      # Estimate outcome for abx = abx_level, infection = 1
      data_level_1 <- data
      data_level_1[[infection_var_name]] <- 1
      data_level_1[[abx_var_name]] <- abx_level
      
      yhat_level_1 <- stats::predict(model1, newdata = data_level_1, type = "response")
      
      if(!att){
        # Regress yhat_level_1 on all other non-mediating variables in subset with infection = 1, call this model2
        data$yhat_level_1 <- yhat_level_1
        sub_inf_1 <- data[data[[infection_var_name]] == 1,]
        
        # covariate list for the exposed (may or may not include the pathogens)
        model2 <- stats::glm(stats::as.formula(paste("yhat_level_1", "~", paste(covariate_list, collapse = "+"))),
                             data = sub_inf_1,
                             family = outcome_type)
        
        # Predict from model2 on everyone and average
        data$ybar_level_1_preds <- stats::predict(model2, newdata = data, type = "response")
        ybar_level_1 <- mean(data$ybar_level_1_preds, na.rm = TRUE)
        
      }else{
        ybar_level_1 <- mean(yhat_level_1[data[[infection_var_name]] == 1], na.rm = TRUE)
      }
      
      # Estimate outcome for abx = abx_level, infection = 0
      data_level_0 <- data
      data_level_0[[infection_var_name]] <- 0
      data_level_1[[abx_var_name]] <- abx_level
      
      yhat_level_0 <- stats::predict(model1, newdata = data_level_0, type = "response")
      
      data$yhat_level_0 <- yhat_level_0
      sub_inf_0 <- data[data[[infection_var_name]] == 0,]
      
      model3 <- stats::glm(stats::as.formula(paste("yhat_level_0", "~", paste(covariate_list, collapse = "+"))),
                           data = sub_inf_0,
                           family = outcome_type)
      
      if(!att){
        data$ybar_level_0_preds <- stats::predict(model3, newdata = data, type = "response")
        ybar_level_0 <- mean(data$ybar_level_0_preds, na.rm = TRUE)
      }else{
        ybar_level_0 <- mean(stats::predict(model3, newdata = data[data[[infection_var_name]] == 1,], type = "response"), na.rm = TRUE)
      }
      
      inf_abx_level <- ybar_level_1 - ybar_level_0
      
      # Add results to dataframe
      results_df[i, "abx_level_inf_1"] <- ybar_level_1
      results_df[i, "abx_level_inf_0"] <- ybar_level_0
      results_df[i, "effect_inf_abx_level"] <- inf_abx_level
      
    }
    
    # Return results dataframe
    return(results_df)
    
  } else {
    
    # ----------------------------------------
    # Traditional G-Computation 
    # ----------------------------------------
    if(site_interaction == "TRUE"){
      model1_formula <- stats::as.formula(paste(laz_var_name, "~",
                                          abx_var_name, "+", infection_var_name, "+",
                                          paste(covariate_list, collapse = "+"), "+",
                                          paste0(site_var_name, "*", abx_var_name), "+",
                                          paste0(infection_var_name, "*", abx_var_name)))
    } else {
      model1_formula <- stats::as.formula(paste(laz_var_name, "~",
                                          abx_var_name, "+", infection_var_name, "+",
                                          paste(covariate_list, collapse = "+"), "+",
                                          paste0(infection_var_name, "*", abx_var_name)))
    }
    
    # Fit model with LAZ ~ infection + abx + covariates + infection*abx
    model1 <- stats::glm(model1_formula,
                   data = data,
                   family = outcome_type)
    
    
    # For each level of abx, predict setting abx = x, infection = 1 & abx = x, infection = 0
    abx_levels <- unique(data[[abx_var_name]])
    results_df <- data.frame(abx_levels = abx_levels,
                             abx_level_inf_1 = vector("numeric", length = length(abx_levels)),
                             abx_level_inf_0 = vector("numeric", length = length(abx_levels)),
                             effect_inf_abx_level = vector("numeric", length = length(abx_levels)))
    
    
    for(i in 1:length(abx_levels)){
      abx_level <- abx_levels[i]
      
      # Abx = abx_level, Infection = 1
      data_level_1 <- data
      data_level_1[[infection_var_name]] <- 1
      data_level_1[[abx_var_name]] <- abx_level
      
      if(!att){
        ybar_level_1 <- mean(stats::predict(model1, newdata = data_01, type = "response"), na.rm = TRUE)
      } else{
        data$yhat_level_1 <- stats::predict(model1, newdata = data_01, type = "response")
        ybar_level_1 <- mean(data$yhat_level_1[data[[infection_var_name]] == 1], na.rm = TRUE)
      }
      
      # Abx = abx_level, Infection = 0
      data_level_0 <- data
      data_level_0[[infection_var_name]] <- 0
      data_level_0[[abx_var_name]] <- abx_level
      
      if(!att){
        ybar_level_0 <- mean(stats::predict(model_1, newdata = data_level_0), na.rm = TRUE)
      } else {
        data$yhat_level_0 <- stats::predict(model_1, newdata = data_level_0, type = "response")
        ybar_level_0 <- mean(data$yhat_level_0[data[[infection_var_name]] == 1], na.rm = TRUE) 
      }
      
      inf_abx_level <- ybar_level_1 - ybar_level_0
      
      # Add results to dataframe
      results_df[i, "abx_level_inf_1"] <- ybar_level_1
      results_df[i, "abx_level_inf_0"] <- ybar_level_0
      results_df[i, "effect_inf_abx_level"] <- inf_abx_level
      
    }
    
    # Return results dataframe
    return(results_df)
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
#' @param att boolean if effect should be estimated among people who would naturally get infection, default TRUE
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
                             outcome_type = "gaussian",
                             att = TRUE){
  
  # Estimands of interest: (!ATT)
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
  # Now doing in data cleaning step
  # control_data <- control_data[control_data[[abx_var_name]] == 0,]
  
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
    
    model1_formula_case <- stats::as.formula(
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

    model1_formula_case <- stats::as.formula(paste(laz_var_name, "~",
                                               abx_var_name, "+", case_var_name, "+",
                                               paste(covariate_list, collapse = "+"), "+",
                                               paste(severity_list, collapse = "+")))
  }

  # Regress LAZ on abx, infection, all severity / non-mediating variables
  model1_case <- stats::glm(model1_formula_case,
                        data = case_data,
                        family = outcome_type)
  
  abx_levels <- unique(case_data[[abx_var_name]])
  
  results_df <- data.frame(abx_levels = abx_levels,
                           abx_level_case = vector("numeric", length = length(abx_levels)),
                           abx_level_control = vector("numeric", length = length(abx_levels)),
                           effect_inf_abx_level = vector("numeric", length = length(abx_levels)))
  
  for(i in 1:length(abx_levels)){
    
    abx_level <- abx_levels[i]
    
    data_abx_level <- case_data
    data_abx_level[[abx_var_name]] <- abx_level
    yhat_abx_level <- stats::predict(model1_case, newdata = data_abx_level, type = "response")
    
    if(!att){
      # Regress yhat_01 on all other non-mediating variables 
      case_data$yhat_abx_level <- yhat_abx_level
      
      model2 <- stats::glm(stats::as.formula(paste("yhat_abx_level", "~", paste(covariate_list, collapse = "+"))),
                           data = case_data,
                           family = outcome_type)
      
      # Predict from model2 on controls
      ybar_abx_level <- mean(stats::predict(model2, newdata = control_data, type = "response"), na.rm = TRUE)
    } else {
      ybar_abx_level <- mean(yhat_abx_level, na.rm = TRUE)
    }
    
    results_df[i,"abx_level_case"] <- ybar_abx_level
    
  }
  
  # ------------------------------------------------------------
  # Part 2: E[Growth | Control] 
  # ------------------------------------------------------------
  if(!att){
    ybar_control <- mean(control_data[[laz_var_name]], na.rm = TRUE)
  }else{
    # fit model of laz_var_name ~ covariate_list in the controls
    # predict from model in cases, avg predictions
    
    # NOTE in this case covariate lists should be the same but leave as option
    model_1_formula_control <- stats::as.formula(
      paste(laz_var_name, "~",
            paste(covariate_list_control, collapse = "+"))
    )
    
    # Regress LAZ on covariates
    model_1_control <- stats::glm(model_1_formula_control,
                                  data = control_data,
                                  family = outcome_type)
    
    yhat_control <- stats::predict(model_1_control,
                              newdata = case_data,
                              type = "response")
    
    ybar_control <- mean(yhat_control, na.rm = TRUE)
  }
  
  results_df$abx_level_control <- rep(ybar_control, nrow(results_df))
  
  # Get effect estimates
  results_df$effect_inf_abx_level <- results_df$abx_level_case - results_df$abx_level_control
  
  return(results_df)
  
}

