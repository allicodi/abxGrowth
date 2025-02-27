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
                             no_etiology_var_name = NA,
                             pathogen_quantity_list = c(
                               "rotavirus_attributable"
                             ),
                             pathogen_attributable_list = NULL,
                             outcome_type = "gaussian",
                             sl.library.outcome = NULL, 
                             sl.library.treatment = NULL,
                             sl.library.missingness = NULL,
                             v_folds = 5,
                             att = TRUE,
                             return_models = FALSE){
  
  if(!is.null(severity_list)){
    
    # -----------------------------------
    # Longitudinal G-Computation 
    # -----------------------------------
    
    ### STEP 0: Create subsets of data for model fitting
    
    # will use to merge later
    data$idx <- 1:nrow(data)
    
    ## Subset 0a: subset shigella (or other infection variable) attr cases
    sub_inf_attr <- data[which(data[[infection_var_name]] == 1),]
    
    # Full data (including missing outcome)
    I_Y_inf_attr <- ifelse(is.na(sub_inf_attr[[laz_var_name]]), 1, 0) #indicator for Y missing
    Y_inf_attr <- sub_inf_attr[[laz_var_name]]
    covariates_inf_attr <- sub_inf_attr[, covariate_list, drop = FALSE]
    severity_inf_attr <- sub_inf_attr[, severity_list, drop = FALSE]
    pathogen_q_inf_attr <- sub_inf_attr[, pathogen_quantity_list, drop = FALSE]
    abx_inf_attr <- sub_inf_attr[, abx_var_name, drop = FALSE]
    
    # Complete data (excluding missing outcome)
    sub_inf_attr_complete <- sub_inf_attr[!is.na(sub_inf_attr[[laz_var_name]]), ]
    
    Y_inf_attr_complete <- sub_inf_attr_complete[[laz_var_name]]
    covariates_inf_attr_complete <- sub_inf_attr_complete[, covariate_list, drop = FALSE]
    severity_inf_attr_complete <- sub_inf_attr_complete[, severity_list, drop = FALSE]
    pathogen_q_inf_attr_complete <- sub_inf_attr_complete[, pathogen_quantity_list, drop = FALSE]
    abx_inf_attr_complete <- sub_inf_attr_complete[, abx_var_name, drop = FALSE]
    
    ## Subset 0b: subset to cases with no etiology
    sub_no_attr <- data[which(data[[infection_var_name]] == 0),]
    
    # if data has var indicating no etiology, use that
    if(!is.na(no_etiology_var_name)){
      sub_no_attr <- sub_no_attr[which(sub_no_attr[[no_etiology_var_name]] == 1),]
    } else {
      # otherwise find which rows have no attr pathogens in pathogen_attributable_list
      sub_no_attr <- sub_no_attr[which(rowSums(sub_no_attr[,pathogen_attributable_list]) == 0),]
    }
    
    I_Y_no_attr <- ifelse(is.na(sub_no_attr[[laz_var_name]]), 1, 0) #indicator for Y missing
    Y_no_attr <- sub_no_attr[[laz_var_name]]
    covariates_no_attr <- sub_no_attr[, covariate_list, drop = FALSE]
    severity_no_attr <- sub_no_attr[, severity_list, drop = FALSE]
    pathogen_q_no_attr <- sub_no_attr[, pathogen_quantity_list, drop = FALSE]
    abx_no_attr <- sub_no_attr[, abx_var_name, drop = FALSE]
    
    sub_no_attr_complete <- sub_no_attr[!is.na(sub_no_attr[[laz_var_name]]), ]
    Y_no_attr_complete <- sub_no_attr_complete[[laz_var_name]]
    covariates_no_attr_complete <- sub_no_attr_complete[, covariate_list, drop = FALSE]
    severity_no_attr_complete <- sub_no_attr_complete[, severity_list, drop = FALSE]
    pathogen_q_no_attr_complete <- sub_no_attr_complete[, pathogen_quantity_list, drop = FALSE]
    abx_no_attr_complete <- sub_no_attr_complete[, abx_var_name, drop = FALSE]
    
    ### STEP 1: Fit outcome models
    
    ## Model 1a: Outcome model in Shigella (or other infection) attributable cases
    # QUESTION should this be outcome type or gaussian?? thesis pseudo-outcome always continuous? but that was dr-learner
    outcome_model_1a <- SuperLearner::SuperLearner(Y = Y_inf_attr_complete, 
                                                   X = data.frame(abx_inf_attr_complete,
                                                                  covariates_inf_attr_complete, 
                                                                  severity_inf_attr_complete, 
                                                                  pathogen_q_inf_attr_complete),
                                                   family = outcome_type, 
                                                   SL.library = sl.library.outcome,
                                                   cvControl = list(V = v_folds))
    
    ## Model 1b: Outcome model in cases with no attribution
    outcome_model_1b <- SuperLearner::SuperLearner(Y = Y_no_attr_complete, 
                                                   X = data.frame(abx_no_attr_complete,
                                                                  covariates_no_attr_complete, 
                                                                  severity_no_attr_complete, 
                                                                  pathogen_q_no_attr_complete),
                                                   family = outcome_type, 
                                                   SL.library = sl.library.outcome,
                                                   cvControl = list(V = v_folds))
    
    ### STEP 2: Fit second stage outcome regressions & predict outcomes
    
    # For each level of abx, predict setting abx = x, infection = 1 & abx = x, infection = 0
    abx_levels <- unique(data[[abx_var_name]])[!is.na(unique(data[[abx_var_name]]))]
    results_df <- data.frame(abx_levels = abx_levels,
                             abx_level_inf_1 = vector("numeric", length = length(abx_levels)),
                             abx_level_inf_0 = vector("numeric", length = length(abx_levels)),
                             effect_inf_abx_level = vector("numeric", length = length(abx_levels)))
    
    outcome_vectors_1a <- data.frame(matrix(ncol = 0, nrow = nrow(data)))
    outcome_vectors_1b <- data.frame(matrix(ncol = 0, nrow = nrow(data)))
    outcome_vectors_2b <- data.frame(matrix(ncol = 0, nrow = nrow(data)))
    
    # Iterate through each abx level
    for(i in 1:length(abx_levels)){
      abx_level <- abx_levels[i]
      
      ## 2a: Predictions for infection = 1 (no 2nd stage model needed)
      #data_level_1 <- data
      #data_level_1[[infection_var_name]] <- 1
      #data_level_1[[abx_var_name]] <- abx_level
    
      # QUESTION should this not be setting to 1?? predictions needed by AIPW supposed to be for everybody? so this is different than what's needed for AIPW?
      #yhat_level_1 <- stats::predict(outcome_model_1a, newdata = data_level_1, type = "response")
      
      data_abx <- data
      data_abx[[abx_var_name]] <- abx_level
      outcome_vectors_1a[,paste0("abx_", abx_level)] <- stats::predict(outcome_model_1a, newdata = data_abx[,c(abx_var_name,
                                                                                                               covariate_list,
                                                                                                               severity_list,
                                                                                                               pathogen_quantity_list)], type = "response")$pred
      
      if(!att){
        exit("Skip non-ATT version for now")
        
        #     # Regress yhat_level_1 on all other non-mediating variables in subset with infection = 1, call this model2
        #     data$yhat_level_1 <- yhat_level_1
        #     sub_inf_1 <- data[data[[infection_var_name]] == 1,]
        #     
        #     # covariate list for the exposed (may or may not include the pathogens)
        #     model2 <- stats::glm(stats::as.formula(paste("yhat_level_1", "~", paste(covariate_list, collapse = "+"))),
        #                          data = sub_inf_1,
        #                          family = outcome_type)
        #     
        #     # Predict from model2 on everyone and average
        #     data$ybar_level_1_preds <- stats::predict(model2, newdata = data, type = "response")
        #     ybar_level_1 <- mean(data$ybar_level_1_preds, na.rm = TRUE)
        
      } else{
        # ybar_level_1 <- mean(yhat_level_1[data[[infection_var_name]] == 1], na.rm = TRUE)
      }
      
        ## 2b: Second stage regression model for no attribution
      
        # Estimate outcome for abx = abx_level, infection = 0
        # data_level_0 <- data
        # data_level_0[[infection_var_name]] <- 0
        # data_level_1[[abx_var_name]] <- abx_level
        # 
        # yhat_level_0 <- stats::predict(outcome_model_1b, newdata = data_level_0, type = "response")
        # data$yhat_level_0 <- yhat_level_0
        
        # QUESTION same as above, sheet says naturally so should this just be data in general? not setting abx level?
        outcome_vectors_1b[,paste0("abx_", abx_level)] <- stats::predict(outcome_model_1b, newdata = data[,c(abx_var_name,
                                                                                                             covariate_list,
                                                                                                             severity_list,
                                                                                                             pathogen_quantity_list)], type = "response")$pred
        
        # Get sub_no_attr_complete with yhat_level_0 predictions using obs_id
        data$set_abx_outcome <- stats::predict(outcome_model_1b, newdata = data_abx[,c(abx_var_name,
                                                                                  covariate_list,
                                                                                  severity_list,
                                                                                  pathogen_quantity_list)], type = "response")$pred
        
        # Take new subset because added set_abx_outcome column
        sub_no_attr_2b <- data[which(data[[infection_var_name]] == 0),]
        
        if(!is.na(no_etiology_var_name)){
          sub_no_attr_2b <- sub_no_attr_2b[which(sub_no_attr_2b[[no_etiology_var_name]] == 1),]
        } else {
          # otherwise find which rows have no attr pathogens in pathogen_attributable_list
          sub_no_attr_2b <- sub_no_attr_2b[which(rowSums(sub_no_attr_2b[,pathogen_attributable_list]) == 0),]
        }
        
        outcome_model_2b <- SuperLearner::SuperLearner(Y = sub_no_attr_2b[['set_abx_outcome']],
                                                       X = sub_no_attr_2b[,covariate_list, drop = FALSE],
                                                       family = outcome_type,
                                                       SL.library = sl.library.outcome,
                                                       cvControl = list(V = v_folds))
        if(!att){
          exit("Skip non-ATT version for now")
          
          # data$ybar_level_0_preds <- stats::predict(model3, newdata = data, type = "response")
          # ybar_level_0 <- mean(data$ybar_level_0_preds, na.rm = TRUE)
        }else{
          # yhat_level_0b <- stats::predict(outcome_model_2b, newdata = sub_inf_attr[,c(covariate_list)], type = "response")
          # ybar_level_0 <- mean(yhat_level_0b, na.rm = TRUE)
          
          # QUESTION same as above not setting abx level?
          outcome_vectors_2b[,paste0("abx_", abx_level)] <- stats::predict(outcome_model_2b, newdata = data[,c(covariate_list)], type = "response")$pred
        }

        # inf_abx_level <- ybar_level_1 - ybar_level_0
        # 
        # # Add results to dataframe
        # results_df[i, "abx_level_inf_1"] <- ybar_level_1
        # results_df[i, "abx_level_inf_0"] <- ybar_level_0
        # results_df[i, "effect_inf_abx_level"] <- inf_abx_level
      
    }
    
    
    ### STEP 3: Propensity models
    
    prop_vectors_1a <- data.frame(matrix(ncol = 0, nrow = nrow(data)))
    prop_vectors_1b <- data.frame(matrix(ncol = 0, nrow = nrow(data)))
    prop_vectors_2a <- data.frame(matrix(ncol = 0, nrow = nrow(data)))
    prop_vectors_3a <- data.frame(matrix(ncol = 0, nrow = nrow(data)))
    prop_vectors_3b <- data.frame(matrix(ncol = 0, nrow = nrow(data)))
    
    ## Part 1: Propensity models for antibiotics
    for(i in 1:length(abx_levels)){
      
      # Abx level to be modeled
      abx_level <- abx_levels[i]
      
      # If not last iteration
      if(i != length(abx_levels)){

        # Leave out previously modeled levels if applicable
        if(i > 1){
          abx_levels_out <- abx_levels[1:(i-1)]
          prop_sub_inf_attr <- sub_inf_attr[-which(sub_inf_attr[[abx_var_name]] %in% abx_levels_out),]
          prop_sub_no_attr <- sub_no_attr[-which(sub_no_attr[[abx_var_name]] %in% abx_levels_out),]
          
        } else{
          prop_sub_inf_attr <- sub_inf_attr
          prop_sub_no_attr <- sub_no_attr
        }
        
        prop_covariates_inf_attr <- prop_sub_inf_attr[, covariate_list, drop = FALSE]
        prop_severity_inf_attr <- prop_sub_inf_attr[, severity_list, drop = FALSE]
        prop_pathogen_inf_attr <- prop_sub_inf_attr[, pathogen_quantity_list, drop = FALSE]
        
        prop_covariates_no_attr <- prop_sub_no_attr[, covariate_list, drop = FALSE]
        prop_severity_no_attr <- prop_sub_no_attr[, severity_list, drop = FALSE]
        prop_pathogen_no_attr <- prop_sub_no_attr[, pathogen_quantity_list, drop = FALSE]
        
        ## 1a. Propensity model for abx shigella attributable
        prop_model_1a <- SuperLearner::SuperLearner(Y = as.numeric(prop_sub_inf_attr[[abx_var_name]] == abx_level),
                                                    X = data.frame(prop_covariates_inf_attr,
                                                                   prop_severity_inf_attr,
                                                                   prop_pathogen_inf_attr),
                                                    newX = data[,c(covariate_list,
                                                                   severity_list,
                                                                   pathogen_quantity_list)], 
                                                    family = stats::binomial(), 
                                                    SL.library = sl.library.treatment,
                                                    cvControl = list(V = v_folds))
        
        # Predictions from full data
        tmp_pred_a <- prop_model_1a$SL.pred
        
        ## 1b. Propensity model for abx no attribution
        prop_model_1b <- SuperLearner::SuperLearner(Y = as.numeric(prop_sub_no_attr[[abx_var_name]] == abx_level),
                                                    X = data.frame(prop_covariates_no_attr,
                                                                   prop_severity_no_attr,
                                                                   prop_pathogen_no_attr),
                                                    newX = data[,c(covariate_list,
                                                                   severity_list,
                                                                   pathogen_quantity_list)], 
                                                    family = stats::binomial(), 
                                                    SL.library = sl.library.treatment,
                                                    cvControl = list(V = v_folds))
        
        # Predictions from full data
        tmp_pred_b <- prop_model_1b$SL.pred
        
        if(i == 1){
          # First prediction
          prop_vectors_1a[,paste0("abx_", abx_level)] <- tmp_pred_a
          prop_vectors_1b[,paste0("abx_", abx_level)] <- tmp_pred_b
        } else{
          # Middle prediction
          prop_vectors_1a[,paste0("abx_", abx_level)] <- tmp_pred_a * (1 - prop_vectors_1a[,ncol(prop_vectors_1a)]) # prediction * previous column 
          prop_vectors_1b[,paste0("abx_", abx_level)] <- tmp_pred_b * (1 - prop_vectors_1a[,ncol(prop_vectors_1b)]) # prediction * previous column
        }

      } else{
        # Last prediction
        prop_vectors_1a[,paste0("abx_", abx_level)] <- 1 - rowSums(prop_vectors_1a)
        prop_vectors_1b[,paste0("abx_", abx_level)] <- 1 - rowSums(prop_vectors_1b)
      }
      
    }
    
    ## Part 2: Propensity models for shigella (or other infection) attribution
    
    # 2a_1 = Shigella Attributable ~ BL Cov
    prop_model_2a_1 <- SuperLearner::SuperLearner(Y = data[[infection_var_name]],
                                                 X = data[, covariate_list, drop = FALSE], 
                                                 family = stats::binomial(),
                                                 SL.library = sl.library.treatment, # QUESTION should this be different than prop model for abx?
                                                 cvControl = list(V = v_folds))
    
    tmp_pred_2a_1 <- prop_model_2a_1$SL.pred
    prop_vectors_2a$inf_attr <- tmp_pred_2a_1
    
    # 2a_2 = No attribution ~ BL Cov | Shigella Attr == 0
    
    sub_no_shig <- data[which(data[[infection_var_name]] == 0),]
    
    # if data has var indicating no etiology, use that. otherwise make a var for it
    if(is.na(no_etiology_var_name)){
      sub_no_shig$no_etiology <- ifelse(rowSums(sub_no_shig[,pathogen_attributable_list]) == 0, 1, 0)
      no_etiology_var_name <- "no_etiology"
    } 
    
    prop_model_2a_2 <- SuperLearner::SuperLearner(Y = sub_no_shig[[no_etiology_var_name]],
                                                  X = sub_no_shig[,covariate_list, drop = FALSE],
                                                  family = stats::binomial(),
                                                  SL.library = sl.library.treatment,
                                                  cvControl = list(V = v_folds))
    
    tmp_pred_2a_2 <- stats::predict(prop_model_2a_2, newdata = data[,covariate_list, drop = FALSE], type = "response")$pred
    prop_vectors_2a$no_attr <- tmp_pred_2a_2 * (1 - prop_vectors_2a$inf_attr)
    
    ## Part 3: Propensity models for missingness
    
    ## Missingness model in infection attributable
    prop_model_3a <- SuperLearner::SuperLearner(Y = I_Y_inf_attr,
                                                X = data.frame(abx_inf_attr,
                                                               covariates_inf_attr,
                                                               severity_inf_attr,
                                                               pathogen_q_inf_attr),
                                                family = stats::binomial(),
                                                SL.library = sl.library.missingness,
                                                cvControl = list(V = v_folds))
    
    ## Missingness model in no etiology 
    prop_model_3b <- SuperLearner::SuperLearner(Y = I_Y_no_attr,
                                                X = data.frame(abx_no_attr,
                                                               covariates_no_attr,
                                                               severity_no_attr,
                                                               pathogen_q_no_attr),
                                                family = stats::binomial(),
                                                SL.library = sl.library.missingness,
                                                cvControl = list(V = v_folds))
    
    # Predict setting each abx level
    for(i in 1:length(abx_levels)){
      abx_level <- abx_levels[i]
      
      pred_data <- data
      pred_data[[abx_var_name]] <- abx_level
      
      prop_vectors_3a[,paste0("abx_", abx_level)] <- stats::predict(prop_model_3a, newdata = pred_data[,c(abx_var_name,
                                                                                                          covariate_list,
                                                                                                          severity_list,
                                                                                                          pathogen_quantity_list)], type = "response")$pred
      prop_vectors_3b[,paste0("abx_", abx_level)] <- stats::predict(prop_model_3b, newdata = pred_data[,c(abx_var_name,
                                                                                                          covariate_list,
                                                                                                          severity_list,
                                                                                                          pathogen_quantity_list)], type = "response")$pred
    }
   
    
    # PUT IT ALL TOGETHER FOR AIPW???
    # Stopped here
    
    # We have dataframes with:
    #outcome_vectors_1a (3vec)
    #outcome_vectors_1b (3vec)
    #outcome_vectors_2b (3vec)
    
    #prop_vectors_1a (3vec)
    #prop_vectors_1b (3vec)
    #prop_vectors_2a (2vec)
    #prop_vectors_3a (3vec)
    #prop_vectors_3b (3vec)
    
    exit("stopped here")
    # Return list with results dataframe (old version, effect est only)
    #return(results_df)
    
  } else {
    
    exit("Skip for now")
    
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

