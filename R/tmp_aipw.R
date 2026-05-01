#' Function for AIPW estimates of effect of infection on growth- two-stage outcome regression
#' 
#' @param data dataframe containing dataset used for gcomp
#' @param laz_var_name name of growth outcome variable
#' @param abx_var_name name of binary antibiotic variable
#' @param infection_var_name name of binary infection variable (for other diarrhea analyses; case_control = FALSE)
#' @param site_var_name name of covariate for site (to exclude from propensity models)
#' @param covariate_list character vector containing names of baseline covariates
#' @param severity_list character vector containing names of severity-related covariates (post-infection). If NULL, use AIPW without second stage regression
#' @param pathogen_quantity_list character vector containing name(s) of pathogen quantity variables in dataset
#' @param pathogen_attributable_list character vector containing name(s) of binary pathogen attributable variables in dataset. Used to identify diarrhea with no etiology if is.null(no_etiology_var_name)
#' @param no_etiology_var_name name of binary variable indicating diarrhea episodes with no known etiology. If null, use pathogen_attributable_list
#' @param outcome_type gaussian or binomial for continuous or binomial outcome
#' @param sl.library.outcome list containing SuperLearner libraries to use for outcome regression 
#' @param sl.library.outcome.2 list containing SuperLearner libraries to use for 2nd stage outcome regression (baseline covariates only)
#' @param sl.library.outcome.case list containing SuperLearner libraries to use for case outcome regression 
#' @param sl.library.outcome.control list containing SuperLearner libraries to use for control outcome regression 
#' @param sl.library.treatment list containing SuperLearner libraries to use for antibiotic prescribing propensity model
#' @param sl.library.infection list containing SuperLearner libraries to use for infection propensity model
#' @param sl.library.missingness list containing SuperLearner libraries to use for outcome missingness model
#' @param v_folds number of cross-validation folds to use in SuperLearner
#' @param return_models boolean return SuperLearner models. Default FALSE.
#' @param first_id_var_name name of variable indicating first_id in data. Used to account for re-enrollment in study. 
#' @param msm boolean indicating use of MSM for effect heterogeneity, default FALSE
#' @param msm_var_name name of variable to use for msm
#' @param msm_formula chatacter vector with formula to use for msm if msm TRUE
#' @param ps_trunc_level numeric value to truncate propensity scores `< ps_trunc_level` or `> 1 - ps_trunc_level`. Default to 0.01
#' @param parsimonious_propensity default TRUE for no etiology (TODO: FILL IN, IDK WHAT THIS MEANS. ALSO CAN WE RENAME IT? THIS IS LONG)
#' @param all_other_diarrhea flag for sensitivity analysis comparison to all_other_diarrhea, default FALSE
#' 
#' @keywords internal
#' 
#' @returns List containing `aipw_other_diarrhea_2` object, models (if return_models = TRUE). `aipw_other_diarrhea_2` object contains dataframe with results, covariance matrix, standard errors.
aipw_other_diarrhea_2 <- function(data,
                                  laz_var_name,
                                  abx_var_name,
                                  infection_var_name,
                                  subinfection_var_name, # = 0 if infection_var_name = 0; = 1, 2, ... otherwise specifying a specific type of infection
                                  site_var_name,
                                  followup_var_names,
                                  covariate_list,
                                  severity_list,
                                  pathogen_quantity_list = NA,
                                  pathogen_attributable_list = NA,
                                  no_etiology_var_name = NA,
                                  outcome_type = "gaussian",
                                  sl.library.outcome = c("SL.glm"),
                                  sl.library.outcome.2 = c("SL.glm"),
                                  sl.library.treatment = c("SL.mean"),
                                  sl.library.infection = c("SL.glm"),
                                  sl.library.missingness = c("SL.glm"),
                                  v_folds = 5,
                                  return_models = FALSE,
                                  first_id_var_name = NULL,
                                  msm = FALSE,
                                  msm_var_name = NULL,
                                  msm_formula = NULL,
                                  ps_trunc_level = 0.01,
                                  parsimonious_propensity = TRUE,
                                  all_other_diarrhea = FALSE){
  
  # ------------------------------------------------------------
  # STEP 0: Create subsets of data for model fitting
  # ------------------------------------------------------------
  
  # Rename abx levels if spaces in them
  abx_levels <- levels(factor(data[[abx_var_name]]))
  abx_levels_new <- ifelse(is.na(abx_levels), NA, gsub("[ /]", "_", abx_levels))
  data[[abx_var_name]] <- factor(data[[abx_var_name]], levels = abx_levels, labels = abx_levels_new)
  
  ## Subset 0a: subset shigella (or other infection variable) attr cases
  inf_attr_idx <- which(data[[infection_var_name]] == 1)
  sub_inf_attr <- data[inf_attr_idx,]
  
  # Full data (including missing outcome)
  I_Y_inf_attr <- ifelse(is.na(sub_inf_attr[[laz_var_name]]), 1, 0) #indicator for Y missing
  Y_inf_attr <- sub_inf_attr[[laz_var_name]]
  
  # remove attributes to avoid error in xgboost
  attributes(Y_inf_attr) <- NULL
  
  covariates_inf_attr <- sub_inf_attr[, covariate_list, drop = FALSE]
  severity_inf_attr <- sub_inf_attr[, severity_list, drop = FALSE]
  pathogen_q_inf_attr <- sub_inf_attr[, pathogen_quantity_list, drop = FALSE]
  abx_inf_attr <- sub_inf_attr[, abx_var_name, drop = FALSE]
  subtype_inf_attr <- sub_inf_attr[, subinfection_var_name, drop = FALSE]
  
  # Complete data (excluding missing outcome)
  sub_inf_attr_complete <- sub_inf_attr[!is.na(sub_inf_attr[[laz_var_name]]), ]
  
  Y_inf_attr_complete <- sub_inf_attr_complete[[laz_var_name]]
  
  # remove attributes to avoid error in xgboost
  attributes(Y_inf_attr_complete) <- NULL
  
  covariates_inf_attr_complete <- sub_inf_attr_complete[, covariate_list, drop = FALSE]
  severity_inf_attr_complete <- sub_inf_attr_complete[, severity_list, drop = FALSE]
  pathogen_q_inf_attr_complete <- sub_inf_attr_complete[, pathogen_quantity_list, drop = FALSE]
  abx_inf_attr_complete <- sub_inf_attr_complete[, abx_var_name, drop = FALSE]
  subtype_inf_attr_complete <- sub_inf_attr_complete[, subinfection_var_name, drop = FALSE]

  ## Subset 0b: subset to cases with no etiology
  if(!is.null(no_etiology_var_name)){
    sub_no_attr <- data[which(data[[no_etiology_var_name]] == 1),]
  } else{
    I_no_attr <- ifelse(data[[infection_var_name]] == 0 & (rowSums(data[,pathogen_attributable_list], na.rm = TRUE) == 0),
                        1, 0)
    data$no_etiology <- I_no_attr
    no_etiology_var_name <- "no_etiology"
    
    sub_no_attr <- data[which(data$no_etiology == 1),]
  }
  
  I_Y_no_attr <- ifelse(is.na(sub_no_attr[[laz_var_name]]), 1, 0) #indicator for Y missing
  Y_no_attr <- sub_no_attr[[laz_var_name]]
  
  # remove attributes to avoid error in xgboost
  attributes(Y_no_attr) <- NULL
  
  covariates_no_attr <- sub_no_attr[, covariate_list, drop = FALSE]
  severity_no_attr <- sub_no_attr[, severity_list, drop = FALSE]
  pathogen_q_no_attr <- sub_no_attr[, pathogen_quantity_list, drop = FALSE]
  abx_no_attr <- sub_no_attr[, abx_var_name, drop = FALSE]
  
  sub_no_attr_complete <- sub_no_attr[!is.na(sub_no_attr[[laz_var_name]]), ]
  Y_no_attr_complete <- sub_no_attr_complete[[laz_var_name]]
  
  # remove attributes to avoid error in xgboost
  attributes(Y_no_attr_complete) <- NULL
  
  covariates_no_attr_complete <- sub_no_attr_complete[, covariate_list, drop = FALSE]
  severity_no_attr_complete <- sub_no_attr_complete[, severity_list, drop = FALSE]
  pathogen_q_no_attr_complete <- sub_no_attr_complete[, pathogen_quantity_list, drop = FALSE]
  abx_no_attr_complete <- sub_no_attr_complete[, abx_var_name, drop = FALSE]
  
  # ------------------------------------------------------------
  # STEP 1: Fit & predict from outcome models
  # ------------------------------------------------------------
  
  ## Model 1a: Outcome model in Shigella (or other infection) attributable cases
  outcome_model_1a <- SuperLearner::SuperLearner(Y = Y_inf_attr_complete, 
                                                 X = data.frame(abx_inf_attr_complete,
                                                                covariates_inf_attr_complete, 
                                                                severity_inf_attr_complete, 
                                                                pathogen_q_inf_attr_complete,
                                                                subtype_inf_attr_complete),
                                                 family = outcome_type, 
                                                 SL.library = sl.library.outcome,
                                                 cvControl = list(V = v_folds),
                                                 control = list(saveCVFitLibrary = TRUE))
  
  ## Model 1b: Outcome model in cases with no attribution
  outcome_model_1b <- SuperLearner::SuperLearner(Y = Y_no_attr_complete, 
                                                 X = data.frame(abx_no_attr_complete,
                                                                covariates_no_attr_complete, 
                                                                severity_no_attr_complete, 
                                                                pathogen_q_no_attr_complete),
                                                 family = outcome_type, 
                                                 SL.library = sl.library.outcome,
                                                 cvControl = list(V = v_folds))
  
  # For each level of abx, predict setting abx = x, infection = 1 & abx = x, infection = 0
  abx_levels <- unique(data[[abx_var_name]])[!is.na(unique(data[[abx_var_name]]))]
  subinfection_var_levels <- unique(data[[subinfection_var_name]])[!is.na(unique(data[[subinfection_var_name]]))]

  outcome_vectors_1b <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))  
  outcome_vectors_2b <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))  
  list_outcome_vectors_1a <- vector(mode = "list", length = length(subinfection_var_levels))
  list_outcome_vectors_2a <- vector(mode = "list", length = length(subinfection_var_levels))

  subinfect_ct <- 0
  for(j in seq_along(subinfection_var_levels)){
    subinfect_ct <- subinfect_ct + 1

    subinfect_var_level <- subinfection_var_levels[j]

    list_outcome_vectors_1a[[subinfect_ct]] <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
    list_outcome_vectors_2a[[subinfect_ct]] <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))

    colnames(outcome_vectors_1b) <- paste0("abx_", abx_levels)
    colnames(outcome_vectors_2b) <- paste0("abx_", abx_levels)
    colnames(list_outcome_vectors_1a[[subinfect_ct]]) <- paste0("abx_", abx_levels)
    colnames(list_outcome_vectors_2a[[subinfect_ct]]) <- paste0("abx_", abx_levels)

    if(return_models){
      stop("models not going to be returned. k thanks.")
      list_outcome_model_1b[[subinfect_ct]] <- vector("list", length = length(abx_levels))
      list_outcome_model_2b[[subinfect_ct]] <- vector("list", length = length(abx_levels))
    }
  
    # Iterate through each abx level
    for (i in 1:length(abx_levels)) {
      abx_level <- abx_levels[i]
      
      data_abx <- data
      data_abx[[abx_var_name]] <- abx_level
      
      ## 2a: Predictions for infection = 1 (no 2nd stage model needed)
      pred_data <- data_abx[, c(abx_var_name, covariate_list, severity_list, pathogen_quantity_list)]
      pred_data[[subinfection_var_name]] <- subinfect_var_level
      list_outcome_vectors_1a[[subinfect_ct]][, i] <- stats::predict(outcome_model_1a, newdata = pred_data, type = "response")$pred
      
      # Get sub_no_attr_complete with yhat_level_0 predictions using obs_id
      data$set_abx_and_subinfect_outcome <-  list_outcome_vectors_1a[[subinfect_ct]][, i]
      
      sub_subinfect_2b <- data[which(data[[subinfection_var_name]] == subinfect_var_level), ]

      if(is.null(sl.library.outcome.2)){
        sl.library.outcome.2 <- sl.library.outcome
      }
      
      ## Model 2b: Second stage regression model for no attribution
      outcome_model_2a <- SuperLearner::SuperLearner(
        Y = sub_subinfect_2b[['set_abx_and_subinfect_outcome']],
        X = sub_subinfect_2b[, covariate_list, drop = FALSE],
        family = outcome_type,
        SL.library = sl.library.outcome.2, 
        cvControl = list(V = v_folds)
      )
      
      list_outcome_vectors_2a[[subinfect_ct]][, i] <- stats::predict(outcome_model_2a, newdata = data[, c(covariate_list)], type = "response")$pred
      
      if(return_models){
        list_outcome_model_2a[[subinfect_ct]][[i]] <- outcome_model_2a
      }

    
      if(subinfect_ct == 1){
        ## 2b: Predictions from 1b + Second stage regression model for no attribution
        
        outcome_vectors_1b[, i] <- stats::predict(outcome_model_1b, newdata = data_abx[, c(abx_var_name,
                                                                                           covariate_list,
                                                                                           severity_list,
                                                                                           pathogen_quantity_list)], type = "response")$pred
        # Get sub_no_attr_complete with yhat_level_0 predictions using obs_id
        data$set_abx_outcome <-  outcome_vectors_1b[, i]
        
        # Take new subset because added set_abx_outcome column
        # sub_no_attr_2b <- data[which(data[[infection_var_name]] == 0), ]
        
        if (!is.na(no_etiology_var_name)) {
          sub_no_attr_2b <- data[which(data[[no_etiology_var_name]] == 1), ]
        } else {
          # otherwise find which rows have no attr pathogens in pathogen_attributable_list
          sub_no_attr_2b <- data[which(rowSums(data[, pathogen_attributable_list]) == 0), ]
        }
        
        if(is.null(sl.library.outcome.2)){
          sl.library.outcome.2 <- sl.library.outcome
        }
        
        ## Model 2b: Second stage regression model for no attribution
        outcome_model_2b <- SuperLearner::SuperLearner(
          Y = sub_no_attr_2b[['set_abx_outcome']],
          X = sub_no_attr_2b[, covariate_list, drop = FALSE],
          family = outcome_type,
          SL.library = sl.library.outcome.2, 
          cvControl = list(V = v_folds)
        )
        
        outcome_vectors_2b[, i] <- stats::predict(outcome_model_2b, newdata = data[, c(covariate_list)], type = "response")$pred
        
        if(return_models){
          list_outcome_model_2b[[subinfect_ct]][[i]] <- outcome_model_2b
        }
      }
    }
  }
  
  if(msm){
    stop("msm not done yet")
    # Subtract predictions
    outcome_msm <- outcome_vectors_1a - outcome_vectors_2b
    
    msm_vectors <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
    colnames(msm_vectors) <- paste0("abx_", abx_levels)
    
    if(return_models){
      msm_model_list <- vector("list", length = length(abx_levels))
      msm_mod_mat_list <- vector("list", length = length(abx_levels))
      cQ_inv_list <- vector("list", length = length(abx_levels))
    }
    
    for(i in 1:length(abx_levels)){
      abx_level <- abx_levels[i]
      abx_level_name <- paste0("abx_", abx_level)
      
      # Fit difference in predictions on msm_formula
      msm_formula_full <- as.formula(paste0(abx_level_name, " ~ ", msm_formula))
      
      # dataframe with difference in outcomes, variable of interest, infection var (for subsetting)
      msm_data <- setNames(
        data.frame(outcome_msm[, i], data[[msm_var_name]], data[[infection_var_name]]), 
        c(abx_level_name, msm_var_name, infection_var_name) 
      )
      
      # compute normalizing matrix needed later for 
      mod_mat <- model.matrix(msm_formula_full, data = msm_data[msm_data[[infection_var_name]] == 1,])
      msm_mod_mat_list[[i]] <- model.matrix(msm_formula_full, data = msm_data)
      
      cQ <- lapply(split(mod_mat, row(mod_mat)), tcrossprod)
      cQ_sum <- Reduce("+", cQ)
      cQ_inv_list[[i]] <- tryCatch(solve(cQ_sum/nrow(mod_mat)), error = function(){ MASS::ginv(cQ_sum/nrow(mod_mat)) })
      
      effect_hetero_msm <- stats::glm(msm_formula_full, 
                                      data = msm_data[msm_data[[infection_var_name]] == 1,],       
                                      family = outcome_type)                                       
      
      msm_vectors[,i] <- stats::predict(effect_hetero_msm, newdata = msm_data, type = 'response')
      
      if(return_models){
        msm_model_list[[i]] <- effect_hetero_msm
      }
      
    }
    
  }
  
  # ------------------------------------------------------------
  # STEP 2: Fit & predict from propensity models
  # ------------------------------------------------------------
  
  # If followup_days related variables, remove from missingness models
  if(!any(is.na(followup_var_names))){
    covariate_list_no_followup <- covariate_list[!(covariate_list %in% followup_var_names)]
  } else{
    covariate_list_no_followup <- covariate_list
  }
  
  # Create matrices to hold predictions from propensity models
  # abx
  prop_vectors_1a <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  prop_vectors_1b <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  # shig attr
  prop_vectors_2a <- data.frame(matrix(ncol = 2, nrow = nrow(data)))
  prop_vectors_2b <- data.frame(matrix(ncol = 2, nrow = nrow(data)))
  # missing
  prop_vectors_3a <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  prop_vectors_3b <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  
  colnames(prop_vectors_1a) <- paste0("abx_", abx_levels)
  colnames(prop_vectors_1b) <- paste0("abx_", abx_levels)
  colnames(prop_vectors_2a) <- c("inf_attr", "no_attr")
  # holds P(inf_attr = 1, subinf = 1 | stuff), P(inf_attr = 1, subinf = 2 | stuff)
  colnames(prop_vectors_2b) <- c("subinf_1", "subinf_2")
  colnames(prop_vectors_3a) <- paste0("abx_", abx_levels)
  colnames(prop_vectors_3b) <- paste0("abx_", abx_levels)
  
  ###############################################
  ## Part 1: Propensity models for antibiotics ##
  ###############################################
  
  if(return_models){
    prop_model_1a_list <- vector("list", length = length(abx_levels))
    prop_model_1b_list <- vector("list", length = length(abx_levels))
  }
  
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
      
      prop_covariates_inf_attr <- prop_sub_inf_attr[, covariate_list , drop = FALSE]
      prop_severity_inf_attr <- prop_sub_inf_attr[, severity_list, drop = FALSE]
      prop_pathogen_inf_attr <- prop_sub_inf_attr[, pathogen_quantity_list, drop = FALSE]
      prop_subinfect_inf_attr <- prop_sub_inf_attr[, subinfection_var_name, drop = FALSE]
      
      prop_covariates_no_attr <- prop_sub_no_attr[, covariate_list , drop = FALSE]
      prop_severity_no_attr <- prop_sub_no_attr[, severity_list, drop = FALSE]
      prop_pathogen_no_attr <- prop_sub_no_attr[, pathogen_quantity_list, drop = FALSE]
    
      
      ## 1a. Propensity model for abx shigella attributable
      prop_model_1a <- SuperLearner::SuperLearner(Y = as.numeric(prop_sub_inf_attr[[abx_var_name]] == abx_level),
                                                  X = data.frame(prop_covariates_inf_attr,
                                                                 prop_severity_inf_attr,
                                                                 prop_pathogen_inf_attr,
                                                                 prop_subinfect_inf_attr),
                                                  newX = data[,c(covariate_list ,
                                                                 severity_list,
                                                                 pathogen_quantity_list,
                                                                 subinfection_var_name)], 
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
                                                  newX = data[,c(covariate_list ,
                                                                 severity_list,
                                                                 pathogen_quantity_list)], 
                                                  family = stats::binomial(), 
                                                  SL.library = sl.library.treatment,
                                                  cvControl = list(V = v_folds))
      
      # Predictions from full data
      tmp_pred_b <- prop_model_1b$SL.pred
      
      # save models in list if returning
      if(return_models){
        prop_model_1a_list[[i]] <- prop_model_1a
        prop_model_1b_list[[i]] <- prop_model_1b
      }
      
      if(i == 1){
        # First prediction
        prop_vectors_1a[,i] <- tmp_pred_a
        prop_vectors_1b[,i] <- tmp_pred_b
      } else{
        # Middle prediction
        # tmp_pred_a * for j in 1:(i-1) (1 - tmp_pred_aj) * TODO
        for(j in 1:(i-1)){
          tmp_pred_a <- tmp_pred_a * (1 - prop_vectors_1a[,j])
          tmp_pred_b <- tmp_pred_b * (1 - prop_vectors_1b[,j])
        }
        prop_vectors_1a[,i] <- tmp_pred_a 
        prop_vectors_1b[,i] <- tmp_pred_b 
      }
      
    } else{
      # Last prediction
      prop_vectors_1a[,i] <- 1 - rowSums(prop_vectors_1a[,1:(ncol(prop_vectors_1a)-1)])
      prop_vectors_1b[,i] <- 1 - rowSums(prop_vectors_1b[,1:(ncol(prop_vectors_1b)-1)])
    }
    
  }
  
  #############################################################################
  ## Part 2: Propensity models for shigella (or other infection) attribution ##
  #############################################################################
  
  # 2a_1 = Shigella Attributable ~ BL Cov
  prop_model_2a_1 <- SuperLearner::SuperLearner(Y = data[[infection_var_name]],
                                                X = data[, covariate_list , drop = FALSE], 
                                                family = stats::binomial(),
                                                SL.library = sl.library.infection, 
                                                cvControl = list(V = v_folds))
  tmp_pred_2a_1 <- prop_model_2a_1$SL.pred
  prop_vectors_2a$inf_attr <- tmp_pred_2a_1

  # 2b_1 = I(subinfection_var_levels == subinfection_var_levels[1]) ~ BL Cov | Shigella Attributable
  prop_model_2b_1 <- SuperLearner::SuperLearner(Y = as.numeric(
                                                      data[[subinfection_var_name]][data[[infection_var_name]] == 1] == subinfection_var_levels[1]
                                                    ),
                                                X = data[data[[infection_var_name]] == 1, covariate_list , drop = FALSE], 
                                                newX = data[, covariate_list, drop = FALSE],
                                                family = stats::binomial(),
                                                SL.library = sl.library.infection, 
                                                cvControl = list(V = v_folds))
  P_subinfect_is_level_1_given_attr <- prop_model_2b_1$SL.pred
  P_subinfect_is_level_1_and_attr <- P_subinfect_is_level_1_given_attr * prop_vectors_2a$inf_attr
  P_subinfect_is_level_2_given_attr <- 1 - P_subinfect_is_level_1_given_attr
  P_subinfect_is_level_2_and_attr <- P_subinfect_is_level_2_given_attr * prop_vectors_2a$inf_attr

  prop_vectors_2b$subinf_1 <- P_subinfect_is_level_1_and_attr
  prop_vectors_2b$subinf_2 <- P_subinfect_is_level_2_and_attr


  # 2a_2 = 1 - 2a_1 OR No attribution ~ BL Cov | Shigella Attr == 0 OR No attribution ~ BL Cov
  if(all_other_diarrhea){
    # comparison to all other diarrhea so prediction is 1 - shigella attr
    prop_model_2a_2 <- NULL
    prop_vectors_2a$no_attr <- 1 - prop_vectors_2a$inf_attr
  } else if(parsimonious_propensity){
    
    # 2a_2 = No attribution ~ BL Cov | Shigella Attr == 0
    sub_no_shig <- data[which(data[[infection_var_name]] == 0),]
    
    # if data has var indicating no etiology, use that. otherwise make a var for it
    if(is.na(no_etiology_var_name)){
      sub_no_shig$no_etiology <- ifelse(rowSums(sub_no_shig[,pathogen_attributable_list]) == 0, 1, 0)
      no_etiology_var_name <- "no_etiology"
    } 
    
    prop_model_2a_2 <- SuperLearner::SuperLearner(Y = sub_no_shig[[no_etiology_var_name]],
                                                  X = sub_no_shig[,covariate_list , drop = FALSE],
                                                  newX = data[,covariate_list , drop = FALSE],
                                                  family = stats::binomial(),
                                                  SL.library = sl.library.infection,
                                                  cvControl = list(V = v_folds))
    
    tmp_pred_2a_2 <- prop_model_2a_2$SL.pred
    prop_vectors_2a$no_attr <- tmp_pred_2a_2 * (1 - prop_vectors_2a$inf_attr)
    
    
  }else{
    # 2a_2 = No attribution ~ BL Cov 
    
    # if data has var indicating no etiology, use that. otherwise make a var for it
    if(is.na(no_etiology_var_name)){
      data$no_etiology <- ifelse(rowSums(data[,pathogen_attributable_list]) == 0, 1, 0)
      no_etiology_var_name <- "no_etiology"
    } 
    
    prop_model_2a_2 <- SuperLearner::SuperLearner(Y = data[[no_etiology_var_name]],
                                                  X = data[,covariate_list , drop = FALSE],
                                                  newX = data[,covariate_list , drop = FALSE],
                                                  family = stats::binomial(),
                                                  SL.library = sl.library.infection,
                                                  cvControl = list(V = v_folds))
    
    prop_vectors_2a$no_attr <- prop_model_2a_2$SL.pred
    
  }
  
  ###############################################
  ## Part 3: Propensity models for missingness ##
  ###############################################
  
  # Data without site
  prop_covariates_inf_attr <- covariates_inf_attr[, colnames(covariates_inf_attr) %in% covariate_list_no_followup , drop = FALSE]
  prop_covariates_no_attr <- covariates_no_attr[, colnames(covariates_no_attr) %in% covariate_list_no_followup , drop = FALSE]
  
  ## Missingness model in infection attributable (if there is missingness)
  if(sum(I_Y_inf_attr) == 0){
    # no missingness; dummy model that returns 0 for all
    prop_model_3a <- list()
    class(prop_model_3a) <- "constant_zero_model"
    
    predict.constant_zero_model <- function(object, newdata, ...) {
      list(pred = rep(0, nrow(newdata)))
    }
    
  } else{
    
    prop_model_3a <- SuperLearner::SuperLearner(Y = I_Y_inf_attr,
                                                X = data.frame(abx_inf_attr,
                                                               prop_covariates_inf_attr,
                                                               severity_inf_attr,
                                                               pathogen_q_inf_attr,
                                                               subtype_inf_attr),
                                                family = stats::binomial(),
                                                SL.library = sl.library.missingness,
                                                cvControl = list(V = v_folds))
  }
  
  ## Missingness model in no etiology 
  
  if(sum(I_Y_no_attr) == 0){
    # no missingness; dummy model that returns 0 for all
    prop_model_ba <- list()
    class(prop_model_3b) <- "constant_zero_model"
    
    predict.constant_zero_model <- function(object, newdata, ...) {
      list(pred = rep(0, nrow(newdata)))
    }
  } else{
    prop_model_3b <- SuperLearner::SuperLearner(Y = I_Y_no_attr,
                                                X = data.frame(abx_no_attr,
                                                               prop_covariates_no_attr,
                                                               severity_no_attr,
                                                               pathogen_q_no_attr),
                                                family = stats::binomial(),
                                                SL.library = sl.library.missingness,
                                                cvControl = list(V = v_folds))
  }
  
  # Predict setting each abx level
  for(i in 1:length(abx_levels)){
    abx_level <- abx_levels[i]
    
    pred_data <- data
    pred_data[[abx_var_name]] <- abx_level
    
    prop_vectors_3a[,i] <- stats::predict(prop_model_3a, newdata = pred_data[,c(abx_var_name,
                                                                                covariate_list_no_followup,
                                                                                severity_list,
                                                                                pathogen_quantity_list,
                                                                                subinfection_var_name)], type = "response")$pred
    prop_vectors_3b[,i] <- stats::predict(prop_model_3b, newdata = pred_data[,c(abx_var_name,
                                                                                covariate_list_no_followup,
                                                                                severity_list,
                                                                                pathogen_quantity_list)], type = "response")$pred
  }
  
  
  # -------------------------------------------------
  # STEP 3: AIPW estimates and  confidence intervals
  # -------------------------------------------------
  
  ## Plug-in estimates
  
  plug_ins_inf_subinfect1 <- colMeans(list_outcome_vectors_2a[[1]][inf_attr_idx, , drop = FALSE])
  plug_ins_inf_subinfect2 <- colMeans(list_outcome_vectors_2a[[2]][inf_attr_idx, , drop = FALSE])
  plug_ins_no_attr <- colMeans(outcome_vectors_2b[inf_attr_idx, , drop = FALSE])
  
  ## Bias corrections
  inf_subinfect1_eifs <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  inf_subinfect2_eifs <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))

  no_attr_eifs <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  
  inf_eifs_msm <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  no_attr_eifs_msm <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  
  colnames(inf_subinfect1_eifs) <- paste0("inf_eif_", abx_levels)
  colnames(inf_subinfect2_eifs) <- paste0("inf_eif_", abx_levels)
  colnames(no_attr_eifs) <- paste0("no_attr_eif_", abx_levels)
  
  colnames(inf_eifs_msm) <- paste0("inf_eif_", abx_levels)
  colnames(no_attr_eifs_msm) <- paste0("no_attr_eif_", abx_levels)
  
  # Truncate any large propensity scores
  if(!is.na(ps_trunc_level)){
    #  any observation that is < ps_trunc_level or > 1-ps_trunc_level should be changed to ps_trunc_level or 1 - ps_trunc_level
    
    truncate_ps <- function(mat, ps_trunc_level){
      mat[mat < ps_trunc_level] <- ps_trunc_level
      mat[mat > 1- ps_trunc_level] <- 1- ps_trunc_level
      return(mat)
    }
    
    prop_vectors_1a <- truncate_ps(prop_vectors_1a, ps_trunc_level)
    prop_vectors_1b <- truncate_ps(prop_vectors_1b, ps_trunc_level)
    prop_vectors_2a <- truncate_ps(prop_vectors_2a, ps_trunc_level)
    prop_vectors_2b <- truncate_ps(prop_vectors_2b, ps_trunc_level)
    prop_vectors_3a <- truncate_ps(prop_vectors_3a, ps_trunc_level)
    prop_vectors_3b <- truncate_ps(prop_vectors_3b, ps_trunc_level)
  }
  
  # 1 - Bias correction for shigella (or other infection) attributable, abx level = a
  I_Inf_subinfect1_1 <- as.numeric(data[[subinfection_var_name]] == subinfection_var_levels[1])
  I_Inf_subinfect2_1 <- as.numeric(data[[subinfection_var_name]] == subinfection_var_levels[2])
  P_Inf_subinfect1_1 <- prop_vectors_2b$subinf_1
  P_Inf_subinfect2_1 <- prop_vectors_2b$subinf_2

  P_Inf_1__Covariates <- prop_vectors_2a[,1]
  I_Inf_1 <- data[[infection_var_name]]
  P_Inf_1 <- mean(prop_vectors_2a[,1])

  for(i in 1:length(abx_levels)){
    abx_level <- abx_levels[i]
    
    I_Abx_a <- as.numeric(data[[abx_var_name]] == abx_level)
    P_Abx_a__Inf_1_Covariates <- prop_vectors_1a[,i]
    
    I_Delta_0 <- as.numeric(!is.na(data[[laz_var_name]])) # Indicator NOT missing
    P_Delta_0__Inf_all <- 1 - prop_vectors_3a[,i]
    
    obs_outcome <- ifelse(is.na(data[[laz_var_name]]), 0, data[[laz_var_name]])  
    Qbar_Inf_1_subinfect1_Abx_a_Covariates <- list_outcome_vectors_1a[[1]][,i]
    Qbar_Inf_1_subinfect2_Abx_a_Covariates <- list_outcome_vectors_1a[[2]][,i]
    
    Qbar_Inf_1_subinfect1_Covariates <- list_outcome_vectors_2a[[1]][,i]
    Qbar_Inf_1_subinfect2_Covariates <- list_outcome_vectors_2a[[2]][,i]
    
    if(!is.null(first_id_var_name)){
      pseudo_n <- mean(P_Delta_0__Inf_all) * P_Inf_1 * length(unique(data[[first_id_var_name]])) # where length(unique(data[[first_id_var_name]])) = number of unique kids in the dataset
    } else{
      pseudo_n <- mean(P_Delta_0__Inf_all) * P_Inf_1 * nrow(data)
    }
    
    # Truncation by 5 / (sqrt(n) ln(n)) (https://academic.oup.com/aje/article/191/9/1640/6580570?login=true)
    truncate_factor <- 5/(sqrt(pseudo_n) * log(pseudo_n))

    full_propensity <- P_Abx_a__Inf_1_Covariates * P_Delta_0__Inf_all
    full_propensity[full_propensity < truncate_factor] <- truncate_factor
    
    eif_vec_inf_subinfect1 <- 
      I_Inf_subinfect1_1 / P_Inf_subinfect1_1 * P_Inf_1__Covariates / P_Inf_1 * 
      (I_Abx_a * I_Delta_0) / full_propensity * (obs_outcome - Qbar_Inf_1_subinfect1_Abx_a_Covariates) +
      I_Inf_subinfect1_1 / P_Inf_subinfect1_1 * P_Inf_1__Covariates / P_Inf_1 * (Qbar_Inf_1_subinfect1_Abx_a_Covariates - Qbar_Inf_1_subinfect1_Covariates) + 
        (I_Inf_1 / P_Inf_1) * (Qbar_Inf_1_subinfect1_Covariates - plug_ins_inf_subinfect1[i])
    
    eif_vec_inf_subinfect2 <- 
      I_Inf_subinfect2_1 / P_Inf_subinfect2_1 * P_Inf_1__Covariates / P_Inf_1 * 
      (I_Abx_a * I_Delta_0) / full_propensity * (obs_outcome - Qbar_Inf_1_subinfect2_Abx_a_Covariates) +
      I_Inf_subinfect2_1 / P_Inf_subinfect2_1 * P_Inf_1__Covariates / P_Inf_1 * (Qbar_Inf_1_subinfect2_Abx_a_Covariates - Qbar_Inf_1_subinfect2_Covariates) + 
        (I_Inf_1 / P_Inf_1) * (Qbar_Inf_1_subinfect2_Covariates - plug_ins_inf_subinfect2[i])
    
    # original
    # eif_vec_inf <- (I_Inf_1 / P_Inf_1) * (I_Abx_a / P_Abx_a__Inf_1_Covariates) * (I_Delta_0 / P_Delta_0__Inf_all) * (obs_outcome - Qbar_Inf_1_Abx_a_Covariates) +
    #   (I_Inf_1 / P_Inf_1) * (Qbar_Inf_1_Abx_a_Covariates - plug_ins_inf[i])
    
    # eif_vec_inf_msm <- (I_Inf_1 / P_Inf_1) * (I_Abx_a * I_Delta_0) / full_propensity * (obs_outcome - Qbar_Inf_1_Abx_a_Covariates)
    
    # 2 - Bias correction for no etiology (some repeats from above for clarity while writing)
    if(is.na(no_etiology_var_name)){
      I_No_attr_1 <- data$no_etiology
    } else{
      I_No_attr_1 <- data[[no_etiology_var_name]]
    }
    
    P_No_attr__Covaritates <- prop_vectors_2a[,2]
    
    P_Inf_1__Covariates <- prop_vectors_2a[,1]
    P_Inf_1 <- mean(prop_vectors_2a[,1])
    
    I_Abx_a <- as.numeric(data[[abx_var_name]] == abx_level)
    P_Abx_a__Covariates <- prop_vectors_1b[,i]
    
    I_Delta_0 <- as.numeric(!is.na(data[[laz_var_name]])) # Indicator NOT missing
    P_Delta_0__No_attr_all <- 1 - prop_vectors_3b[,i]
    
    obs_outcome <- ifelse(is.na(data[[laz_var_name]]), 0, data[[laz_var_name]])  
    Qbar_No_attr_Abx_a_Covariates <- outcome_vectors_1b[,i]
    
    Qbar_No_attr_Covariates <- outcome_vectors_2b[,i]
    
    I_Inf_1 <- data[[infection_var_name]]
    P_Inf_1 <- mean(prop_vectors_2a[,1])
    
    # Truncation by 5 / (sqrt(n) ln(n)) (https://academic.oup.com/aje/article/191/9/1640/6580570?login=true)
    if(!is.null(first_id_var_name)){
      pseudo_n <- mean(P_Delta_0__No_attr_all) * P_Inf_1 * length(unique(data[[first_id_var_name]])) # where length(unique(data[[first_id_var_name]])) = number of unique kids in the dataset
    } else{
      pseudo_n <- mean(P_Delta_0__No_attr_all) * P_Inf_1 * nrow(data)
    }
    truncate_factor <- 5/(sqrt(pseudo_n) * log(pseudo_n))
    full_propensity <- P_Abx_a__Covariates *  P_Delta_0__No_attr_all
    full_propensity[full_propensity < truncate_factor] <- truncate_factor
    # Also truncate these guys by themselves
    P_No_attr__Covaritates[P_No_attr__Covaritates < truncate_factor] <- truncate_factor

    # original
    # eif_vec_no_attr <- (I_No_attr_1 / P_No_attr__Covaritates) * (P_Inf_1__Covariates / P_Inf_1) * (I_Abx_a / P_Abx_a__Covariates) * (I_Delta_0 / P_Delta_0__No_attr_all) * (obs_outcome - Qbar_No_attr_Abx_a_Covariates) +
    #   (I_No_attr_1 / P_No_attr__Covaritates) * (P_Inf_1__Covariates / P_Inf_1) * (Qbar_No_attr_Abx_a_Covariates - Qbar_No_attr_Covariates) + 
    #   (I_Inf_1 / P_Inf_1) * (Qbar_No_attr_Covariates - plug_ins_no_attr[i]) 
    
    eif_vec_no_attr <- (I_No_attr_1 / P_No_attr__Covaritates) * (P_Inf_1__Covariates / P_Inf_1) * (I_Abx_a * I_Delta_0) / full_propensity * (obs_outcome - Qbar_No_attr_Abx_a_Covariates) +
      (I_No_attr_1 / P_No_attr__Covaritates) * (P_Inf_1__Covariates / P_Inf_1) * (Qbar_No_attr_Abx_a_Covariates - Qbar_No_attr_Covariates) + 
      (I_Inf_1 / P_Inf_1) * (Qbar_No_attr_Covariates - plug_ins_no_attr[i]) 
    
    # eif_vec_no_attr_msm <- (I_No_attr_1 / P_No_attr__Covaritates) * (P_Inf_1__Covariates / P_Inf_1) * (I_Abx_a * I_Delta_0) / full_propensity * (obs_outcome - Qbar_No_attr_Abx_a_Covariates) + 
      # (I_No_attr_1 / P_No_attr__Covaritates) * (P_Inf_1__Covariates / P_Inf_1) * (Qbar_No_attr_Abx_a_Covariates - Qbar_No_attr_Covariates) 

    
    inf_subinfect1_eifs[,i] <- eif_vec_inf_subinfect1
    inf_subinfect2_eifs[,i] <- eif_vec_inf_subinfect2
    no_attr_eifs[,i] <- eif_vec_no_attr
    # inf_eifs_msm[,i] <- eif_vec_inf_msm
    # no_attr_eifs_msm[,i] <- eif_vec_no_attr_msm
  }
  
  # Add bias correction
  aipw_inf_subinfect1 <- plug_ins_inf_subinfect1 + colMeans(inf_subinfect1_eifs)
  aipw_inf_subinfect2 <- plug_ins_inf_subinfect2 + colMeans(inf_subinfect2_eifs)
  aipw_no_attr <- plug_ins_no_attr + colMeans(no_attr_eifs)
  
  eif_matrix <- cbind(inf_subinfect1_eifs, inf_subinfect2_eifs, no_attr_eifs)
  cov_matrix <- stats::cov(eif_matrix)
  
  # eif_matrix_msm <- cbind(inf_eifs_msm, no_attr_eifs_msm)
  # eif_effect_msm_list <- vector(mode = "list", length = length(abx_levels))
  
  # # Finish calculating MSM efficient influence function
  # for(i in 1:length(abx_levels)){
    
  #   # Use EIFs to get variance for effect
  #   idx_1 <- i
  #   idx_2 <- length(abx_levels) + i
    
  #   gradient <- rep(0, length(abx_levels)*2)
  #   gradient[idx_1] <- 1
  #   gradient[idx_2] <- -1
    
  #   gradient <- matrix(gradient, ncol = 1)
    
  #   eif_effect_msm_unnormed <- diag(c(as.matrix(eif_matrix_msm) %*% gradient)) %*% msm_mod_mat_list[[i]] + 
  #     + diag((I_Inf_1 / P_Inf_1) * (outcome_msm[,i] - msm_vectors[,i])) %*% msm_mod_mat_list[[i]]
  #   eif_effect_msm_list[[i]] <- eif_effect_msm_unnormed %*% cQ_inv_list[[i]]
  
  # Get id for each participant and recreate EIFs based on this if present
  if(!is.null(first_id_var_name)){
    first_id_eif_matrix <- cbind(data.frame(first_id = data[[first_id_var_name]]), eif_matrix)
    first_id_eif_matrix <- aggregate(. ~ first_id, data = first_id_eif_matrix, FUN = sum)
    
    # restore column names 
    colnames(first_id_eif_matrix)[-c(1)] <- colnames(eif_matrix)
    
    # Get first ids to return at end for piecing things back together afterwards
    eif_first_ids <- first_id_eif_matrix[,1]
    
    scaled_matrix <- first_id_eif_matrix[,-c(1)] * (nrow(first_id_eif_matrix) / nrow(eif_matrix))
    
    # scaled_matrix_msm_list <- lapply(eif_effect_msm_list, function(eif_effect_msm){
    #   first_id_eif_matrix_msm <- cbind(data.frame(first_id = data[[first_id_var_name]]), eif_effect_msm)
    #   first_id_eif_matrix_msm <- aggregate(. ~ first_id, data = first_id_eif_matrix_msm, FUN = sum)
      
    #   return(first_id_eif_matrix_msm[,-c(1)] * (nrow(first_id_eif_matrix_msm) / nrow(eif_effect_msm)))
    # })
    
  } else{
    eif_first_ids <- NULL
    scaled_matrix <- eif_matrix
    # scaled_matrix_msm_list <- eif_effect_msm_list
  }
  
  aipws_effect_subinfect1 <- vector("numeric", length = length(abx_levels))
  aipws_effect_subinfect2 <- vector("numeric", length = length(abx_levels))
  # aipw_msm <- vector("list", length = length(abx_levels))
  eifs_effect_subinfect1 <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(scaled_matrix)))
  eifs_effect_subinfect2 <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(scaled_matrix)))
  
  names(aipws_effect_subinfect1) <- paste0("effect_", abx_levels)
  names(aipws_effect_subinfect2) <- paste0("effect_", abx_levels)
  # names(aipw_msm) <- paste0("msm_", abx_levels)
  
  colnames(eifs_effect_subinfect1) <- paste0("effect_", abx_levels)
  colnames(eifs_effect_subinfect2) <- paste0("effect_", abx_levels)
  
  # Compute effects for all levels of abx
  
  for(i in 1:length(abx_levels)){
    
    # Get AIPW for Effect
    aipws_effect_subinfect1[i] <- aipw_inf_subinfect1[i] - aipw_no_attr[i]
    aipws_effect_subinfect2[i] <- aipw_inf_subinfect2[i] - aipw_no_attr[i]
    
    # Use EIFs to get variance for effect
    idx_1_subinfect1 <- i
    idx_1_subinfect2 <- length(abx_levels) + i
    idx_2 <- 2*length(abx_levels) + i
    
    gradient_subinfect1 <- rep(0, length(abx_levels)*3)
    gradient_subinfect1[idx_1_subinfect1] <- 1
    gradient_subinfect1[idx_2] <- -1
    gradient_subinfect1 <- matrix(gradient_subinfect1, ncol = 1)

    gradient_subinfect2 <- rep(0, length(abx_levels)*3)
    gradient_subinfect2[idx_1_subinfect2] <- 1
    gradient_subinfect2[idx_2] <- -1
    gradient_subinfect2 <- matrix(gradient_subinfect2, ncol = 1)
    
    eif_effect_subinfect1 <- as.numeric(as.matrix(scaled_matrix) %*% gradient_subinfect1)
    eifs_effect_subinfect1[,i] <- eif_effect_subinfect1

    eif_effect_subinfect2 <- as.numeric(as.matrix(scaled_matrix) %*% gradient_subinfect2)
    eifs_effect_subinfect2[,i] <- eif_effect_subinfect2
    
    # aipw_msm[[i]] <- list(
    #   coefficients = coef(msm_model_list[[i]]) + colMeans(scaled_matrix_msm_list[[i]]),
    #   cov = stats::cov(scaled_matrix_msm_list[[i]]) / nrow(scaled_matrix_msm_list[[i]])
    # )
    
  }
  
  # Put all pt ests into dataframe same format as original gcomp analysis
  results_df <- data.frame(abx_levels = abx_levels,
                           abx_level_inf_subinfect1_1 = aipw_inf_subinfect1,
                           abx_level_inf_subinfect2_1 = aipw_inf_subinfect2,
                           abx_level_inf_0 = aipw_no_attr,
                           effect_inf_subinfect1_abx_level = aipws_effect_subinfect1,
                           effect_inf_subinfect2_abx_level = aipws_effect_subinfect2)
  
  # Add EIFs for effects to EIF matrix
  eif_matrix_scaled <- cbind(scaled_matrix, eifs_effect_subinfect1, eifs_effect_subinfect2)
  cov_matrix <- stats::cov(eif_matrix_scaled)
  eif_hat <- sqrt( diag(cov_matrix) / nrow(eif_matrix_scaled) )
  
  # Get marginal effect estimates for MSM (if applicable) and create results object
  if(msm){
    marginal_effect_estimates <- colMeans(msm_vectors[inf_attr_idx, , drop = FALSE])
    
    results_object <- list(results_df = results_df,
                           plug_ins_inf = plug_ins_inf,
                           plug_ins_no_attr = plug_ins_no_attr,
                           eif_matrix = eif_matrix_scaled,
                           eif_first_ids = eif_first_ids,
                           se = eif_hat,
                           aipw_msm = aipw_msm,
                           marginal_effect_estimates = marginal_effect_estimates)
  } else{
    results_object <- list(results_df = results_df,
                           plug_ins_inf_subinfect1 = plug_ins_inf_subinfect1,
                           plug_ins_inf_subinfect2 = plug_ins_inf_subinfect2,
                           plug_ins_no_attr = plug_ins_no_attr,
                           eif_matrix = eif_matrix_scaled,
                           eif_first_ids = eif_first_ids,
                           se = eif_hat)
  }
  
  class(results_object) <- "aipw_other_diarrhea_2"
  
  if(return_models){
    # Make list of models
    if(msm){
      aipw_models <- list(outcome_model_1a = outcome_model_1a,
                          outcome_model_1b = outcome_model_1b,
                          outcome_model_2b_list = outcome_model_2b_list,
                          prop_model_1a_list = prop_model_1a_list,
                          prop_model_1b_list = prop_model_1b_list,
                          prop_model_2a_1 = prop_model_2a_1,
                          prop_model_2a_2 = prop_model_2a_2,
                          prop_model_3a = prop_model_3a,
                          prop_model_3b = prop_model_3b,
                          msm_model_list = msm_model_list)
    } else{
      aipw_models <- list(outcome_model_1a = outcome_model_1a,
                          outcome_model_1b = outcome_model_1b,
                          outcome_model_2b_list = outcome_model_2b_list,
                          prop_model_1a_list = prop_model_1a_list,
                          prop_model_1b_list = prop_model_1b_list,
                          prop_model_2a_1 = prop_model_2a_1,
                          prop_model_2a_2 = prop_model_2a_2,
                          prop_model_3a = prop_model_3a,
                          prop_model_3b = prop_model_3b)
    }
    
    return(list(results_object = results_object,
                aipw_models = aipw_models))
  } else{
    # Return NULL models
    return(list(results_object = results_object,
                aipw_models = NULL))
  }
  
  return(list(results_df = results_df,
              eif_matrix = eif_matrix_scaled,
              eif_first_ids = eif_first_ids,
              se = eif_hat))
}



aipw_case_control <- function(data,
                              laz_var_name,
                              abx_var_name,
                              case_var_name,
                              subinfection_var_name,
                              site_var_name,
                              followup_var_names,
                              covariate_list,
                              severity_list,
                              pathogen_quantity_list = NULL,
                              outcome_type = "gaussian",
                              sl.library.outcome.case = c("SL.glm"),
                              sl.library.outcome.control = c("SL.glm"),
                              sl.library.outcome.2 = c("SL.glm"),
                              sl.library.treatment = c("SL.mean"),
                              sl.library.infection = c("SL.glm"),
                              sl.library.missingness.case = c("SL.glm"),
                              sl.library.missingness.control = c("SL.glm"),
                              v_folds = 5,
                              return_models = FALSE,
                              first_id_var_name = NULL,
                              msm = FALSE,
                              msm_var_name = NULL,
                              msm_formula = NULL,
                              ps_trunc_level = 0.01){

  if(msm){
    stop("msm not done yet for subtype-specific case-control analysis")
  }

  if(return_models){
    stop("return_models not done yet for subtype-specific case-control analysis")
  }

  # ------------------------------------------------------------
  # STEP 0: Create subsets of data for model fitting
  # ------------------------------------------------------------

  # Rename abx levels if spaces in them
  abx_levels <- levels(factor(data[[abx_var_name]]))
  abx_levels_new <- ifelse(is.na(abx_levels), NA, gsub("[ /]", "_", abx_levels))
  data[[abx_var_name]] <- factor(data[[abx_var_name]], levels = abx_levels, labels = abx_levels_new)

  # Get case vs control data
  case_data <- data[data[[case_var_name]] == 1, ]
  control_data <- data[data[[case_var_name]] == 0, ]

  case_data_idx <- which(data[[case_var_name]] == 1)
  control_data_idx <- which(data[[case_var_name]] == 0)

  # Subinfection levels among cases only
  subinfection_var_levels <- unique(case_data[[subinfection_var_name]])
  subinfection_var_levels <- subinfection_var_levels[!is.na(subinfection_var_levels)]

  if(length(subinfection_var_levels) != 2){
    stop("This function currently assumes exactly two non-missing subinfection levels among cases.")
  }

  # Case data prep
  I_Y_case <- ifelse(is.na(case_data[[laz_var_name]]), 1, 0)
  Y_case <- case_data[[laz_var_name]]
  attributes(Y_case) <- NULL

  covariates_case <- case_data[, covariate_list, drop = FALSE]
  severity_case <- case_data[, severity_list, drop = FALSE]
  pathogen_q_case <- case_data[, pathogen_quantity_list, drop = FALSE]
  abx_case <- case_data[, abx_var_name, drop = FALSE]
  subtype_case <- case_data[, subinfection_var_name, drop = FALSE]

  # Complete case data
  case_data_complete <- case_data[!is.na(case_data[[laz_var_name]]), ]

  Y_case_complete <- case_data_complete[[laz_var_name]]
  attributes(Y_case_complete) <- NULL

  covariates_case_complete <- case_data_complete[, covariate_list, drop = FALSE]
  severity_case_complete <- case_data_complete[, severity_list, drop = FALSE]
  pathogen_q_case_complete <- case_data_complete[, pathogen_quantity_list, drop = FALSE]
  abx_case_complete <- case_data_complete[, abx_var_name, drop = FALSE]
  subtype_case_complete <- case_data_complete[, subinfection_var_name, drop = FALSE]

  # Control data prep
  I_Y_control <- ifelse(is.na(control_data[[laz_var_name]]), 1, 0)
  Y_control <- control_data[[laz_var_name]]
  attributes(Y_control) <- NULL

  covariates_control <- control_data[, covariate_list, drop = FALSE]

  # Complete control data
  control_data_complete <- control_data[!is.na(control_data[[laz_var_name]]), ]

  Y_control_complete <- control_data_complete[[laz_var_name]]
  attributes(Y_control_complete) <- NULL

  covariates_control_complete <- control_data_complete[, covariate_list, drop = FALSE]

  # ------------------------------------------------------------
  # STEP 1: Fit & predict from outcome models
  # ------------------------------------------------------------

  ## Model 1a: Outcome model in cases, including subtype
  outcome_model_1a <- SuperLearner::SuperLearner(
    Y = Y_case_complete,
    X = data.frame(
      abx_case_complete,
      covariates_case_complete,
      severity_case_complete,
      pathogen_q_case_complete,
      subtype_case_complete
    ),
    family = outcome_type,
    SL.library = sl.library.outcome.case,
    cvControl = list(V = v_folds)
  )

  ## Model 1b: Outcome model in controls
  outcome_model_1b <- SuperLearner::SuperLearner(
    Y = Y_control_complete,
    X = data.frame(covariates_control_complete),
    family = outcome_type,
    SL.library = sl.library.outcome.control,
    cvControl = list(V = v_folds)
  )

  # Antibiotic levels among cases
  abx_levels <- unique(case_data[[abx_var_name]])
  abx_levels <- abx_levels[!is.na(abx_levels)]

  # Storage
  list_outcome_vectors_1a <- vector("list", length = length(subinfection_var_levels))
  list_outcome_vectors_2a <- vector("list", length = length(subinfection_var_levels))

  outcome_vectors_1b <- data.frame(matrix(ncol = 1, nrow = nrow(data)))
  colnames(outcome_vectors_1b) <- "control"

  for(s in seq_along(subinfection_var_levels)){

    subinfect_level <- subinfection_var_levels[s]

    list_outcome_vectors_1a[[s]] <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
    list_outcome_vectors_2a[[s]] <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))

    colnames(list_outcome_vectors_1a[[s]]) <- paste0("abx_", abx_levels)
    colnames(list_outcome_vectors_2a[[s]]) <- paste0("abx_", abx_levels)

    for(i in seq_along(abx_levels)){

      abx_level <- abx_levels[i]

      data_abx <- data
      data_abx[[abx_var_name]] <- abx_level

      pred_data <- data_abx[, c(
        abx_var_name,
        covariate_list,
        severity_list,
        pathogen_quantity_list
      ), drop = FALSE]

      pred_data[[subinfection_var_name]] <- subinfect_level

      # First-stage subtype-specific case prediction
      list_outcome_vectors_1a[[s]][, i] <- stats::predict(
        outcome_model_1a,
        newdata = pred_data,
        type = "response"
      )$pred

      # Second-stage regression onto baseline covariates among cases of this subtype
      data$set_abx_and_subinfect_outcome <- list_outcome_vectors_1a[[s]][, i]

      sub_subinfect_2a <- data[
        data[[case_var_name]] == 1 & data[[subinfection_var_name]] == subinfect_level,
      ]

      if(is.null(sl.library.outcome.2)){
        sl.library.outcome.2 <- sl.library.outcome.case
      }

      outcome_model_2a <- SuperLearner::SuperLearner(
        Y = sub_subinfect_2a[["set_abx_and_subinfect_outcome"]],
        X = sub_subinfect_2a[, covariate_list, drop = FALSE],
        family = outcome_type,
        SL.library = sl.library.outcome.2,
        cvControl = list(V = v_folds)
      )

      list_outcome_vectors_2a[[s]][, i] <- stats::predict(
        outcome_model_2a,
        newdata = data[, covariate_list, drop = FALSE],
        type = "response"
      )$pred
    }
  }

  # Control outcome predictions, shared across subtype/abx levels
  outcome_vectors_1b[, 1] <- stats::predict(
    outcome_model_1b,
    newdata = data[, covariate_list, drop = FALSE],
    type = "response"
  )$pred

  # ------------------------------------------------------------
  # STEP 2: Fit & predict from propensity models
  # ------------------------------------------------------------

  if(!any(is.na(followup_var_names))){
    covariate_list_no_followup <- covariate_list[!(covariate_list %in% followup_var_names)]
  } else{
    covariate_list_no_followup <- covariate_list
  }

  # Storage
  prop_vectors_1a <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  prop_vectors_2a <- data.frame(matrix(ncol = 1, nrow = nrow(data)))
  prop_vectors_2b <- data.frame(matrix(ncol = 2, nrow = nrow(data)))
  prop_vectors_3a <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  prop_vectors_3b <- data.frame(matrix(ncol = 1, nrow = nrow(data)))

  colnames(prop_vectors_1a) <- paste0("abx_", abx_levels)
  colnames(prop_vectors_2a) <- "case"
  colnames(prop_vectors_2b) <- c("subinf_1", "subinf_2")
  colnames(prop_vectors_3a) <- paste0("abx_", abx_levels)
  colnames(prop_vectors_3b) <- "control"

  ###############################################
  ## Part 1: Propensity models for antibiotics ##
  ###############################################

  for(i in seq_along(abx_levels)){

    abx_level <- abx_levels[i]

    if(i != length(abx_levels)){

      if(i > 1){
        abx_levels_out <- abx_levels[1:(i - 1)]
        prop_sub_case <- case_data[-which(case_data[[abx_var_name]] %in% abx_levels_out), ]
      } else{
        prop_sub_case <- case_data
      }

      prop_covariates_case <- prop_sub_case[, covariate_list, drop = FALSE]
      prop_severity_case <- prop_sub_case[, severity_list, drop = FALSE]
      prop_pathogen_case <- prop_sub_case[, pathogen_quantity_list, drop = FALSE]
      prop_subtype_case <- prop_sub_case[, subinfection_var_name, drop = FALSE]

      ## Antibiotic propensity among cases, conditioning on subtype
      prop_model_1a <- SuperLearner::SuperLearner(
        Y = as.numeric(prop_sub_case[[abx_var_name]] == abx_level),
        X = data.frame(
          prop_covariates_case,
          prop_severity_case,
          prop_pathogen_case,
          prop_subtype_case
        ),
        newX = data[, c(
          covariate_list,
          severity_list,
          pathogen_quantity_list,
          subinfection_var_name
        ), drop = FALSE],
        family = stats::binomial(),
        SL.library = sl.library.treatment,
        cvControl = list(V = v_folds)
      )

      tmp_pred_a <- prop_model_1a$SL.pred

      if(i == 1){
        prop_vectors_1a[, i] <- tmp_pred_a
      } else{
        for(j in 1:(i - 1)){
          tmp_pred_a <- tmp_pred_a * (1 - prop_vectors_1a[, j])
        }
        prop_vectors_1a[, i] <- tmp_pred_a
      }

    } else{
      prop_vectors_1a[, i] <- 1 - rowSums(prop_vectors_1a[, 1:(ncol(prop_vectors_1a) - 1), drop = FALSE])
    }
  }

  ################################################################
  ## Part 2: Propensity models for case and subtype attribution ##
  ################################################################

  # Case ~ baseline covariates
  prop_model_2a <- SuperLearner::SuperLearner(
    Y = data[[case_var_name]],
    X = data[, covariate_list, drop = FALSE],
    family = stats::binomial(),
    SL.library = sl.library.infection,
    cvControl = list(V = v_folds)
  )

  prop_vectors_2a[, 1] <- prop_model_2a$SL.pred

  # Subtype among cases ~ baseline covariates
  prop_model_2b_1 <- SuperLearner::SuperLearner(
    Y = as.numeric(case_data[[subinfection_var_name]] == subinfection_var_levels[1]),
    X = case_data[, covariate_list, drop = FALSE],
    newX = data[, covariate_list, drop = FALSE],
    family = stats::binomial(),
    SL.library = sl.library.infection,
    cvControl = list(V = v_folds)
  )

  P_subinfect_is_level_1_given_case <- prop_model_2b_1$SL.pred
  P_subinfect_is_level_1_and_case <- P_subinfect_is_level_1_given_case * prop_vectors_2a[, 1]

  P_subinfect_is_level_2_given_case <- 1 - P_subinfect_is_level_1_given_case
  P_subinfect_is_level_2_and_case <- P_subinfect_is_level_2_given_case * prop_vectors_2a[, 1]

  prop_vectors_2b$subinf_1 <- P_subinfect_is_level_1_and_case
  prop_vectors_2b$subinf_2 <- P_subinfect_is_level_2_and_case

  ###############################################
  ## Part 3: Propensity models for missingness ##
  ###############################################

  covariates_case_no_followup <- case_data[, covariate_list_no_followup, drop = FALSE]
  covariates_control_no_followup <- control_data[, covariate_list_no_followup, drop = FALSE]

  ## Missingness model in cases
  if(sum(I_Y_case) == 0){

    prop_model_3a <- list()
    class(prop_model_3a) <- "constant_zero_model"

    predict.constant_zero_model <- function(object, newdata, ...) {
      list(pred = rep(0, nrow(newdata)))
    }

  } else{

    prop_model_3a <- SuperLearner::SuperLearner(
      Y = I_Y_case,
      X = data.frame(
        abx_case,
        covariates_case_no_followup,
        severity_case,
        pathogen_q_case,
        subtype_case
      ),
      family = stats::binomial(),
      SL.library = sl.library.missingness.case,
      cvControl = list(V = v_folds)
    )
  }

  ## Missingness model in controls
  if(sum(I_Y_control) == 0){

    prop_model_3b <- list()
    class(prop_model_3b) <- "constant_zero_model"

    predict.constant_zero_model <- function(object, newdata, ...) {
      list(pred = rep(0, nrow(newdata)))
    }

  } else{

    prop_model_3b <- SuperLearner::SuperLearner(
      Y = I_Y_control,
      X = data.frame(covariates_control_no_followup),
      family = stats::binomial(),
      SL.library = sl.library.missingness.control,
      cvControl = list(V = v_folds)
    )
  }

  # Predict case missingness setting each antibiotic level
  for(i in seq_along(abx_levels)){

    abx_level <- abx_levels[i]

    pred_data <- data
    pred_data[[abx_var_name]] <- abx_level

    prop_vectors_3a[, i] <- stats::predict(
      prop_model_3a,
      newdata = pred_data[, c(
        abx_var_name,
        covariate_list_no_followup,
        severity_list,
        pathogen_quantity_list,
        subinfection_var_name
      ), drop = FALSE],
      type = "response"
    )$pred
  }

  # Control missingness predictions
  prop_vectors_3b[, 1] <- stats::predict(
    prop_model_3b,
    newdata = data[, covariate_list_no_followup, drop = FALSE],
    type = "response"
  )$pred

  # -------------------------------------------------
  # STEP 3: AIPW estimates and confidence intervals
  # -------------------------------------------------

  ## Plug-in estimates
  plug_ins_case_subinfect1 <- colMeans(list_outcome_vectors_2a[[1]][case_data_idx, , drop = FALSE])
  plug_ins_case_subinfect2 <- colMeans(list_outcome_vectors_2a[[2]][case_data_idx, , drop = FALSE])
  plug_ins_control <- mean(outcome_vectors_1b[case_data_idx, 1])

  ## EIF storage
  case_subinfect1_eifs <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  case_subinfect2_eifs <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  control_eifs <- data.frame(matrix(ncol = 1, nrow = nrow(data)))

  colnames(case_subinfect1_eifs) <- paste0("case_subinfect1_eif_", abx_levels)
  colnames(case_subinfect2_eifs) <- paste0("case_subinfect2_eif_", abx_levels)
  colnames(control_eifs) <- "control_eif"

  # Truncate propensity scores
  if(!is.na(ps_trunc_level)){

    truncate_ps <- function(mat, ps_trunc_level, upper = TRUE, lower = TRUE){
      if(lower){
        mat[mat < ps_trunc_level] <- ps_trunc_level
      }

      if(upper){
        mat[mat > 1 - ps_trunc_level] <- 1 - ps_trunc_level
      }

      return(mat)
    }

    prop_vectors_1a <- truncate_ps(prop_vectors_1a, ps_trunc_level, upper = FALSE, lower = TRUE)
    prop_vectors_2a <- truncate_ps(prop_vectors_2a, ps_trunc_level, upper = TRUE, lower = TRUE)
    prop_vectors_2b <- truncate_ps(prop_vectors_2b, ps_trunc_level, upper = FALSE, lower = TRUE)
    prop_vectors_3a <- truncate_ps(prop_vectors_3a, ps_trunc_level, upper = TRUE, lower = FALSE)
    prop_vectors_3b <- truncate_ps(prop_vectors_3b, ps_trunc_level, upper = TRUE, lower = FALSE)
  }

  I_Case <- data[[case_var_name]]
  P_Case <- mean(prop_vectors_2a[, 1])
  P_Case__Covariates <- prop_vectors_2a[, 1]

  I_Case_subinfect1 <- as.numeric(data[[case_var_name]] == 1 & data[[subinfection_var_name]] == subinfection_var_levels[1])
  I_Case_subinfect2 <- as.numeric(data[[case_var_name]] == 1 & data[[subinfection_var_name]] == subinfection_var_levels[2])

  P_Case_subinfect1__Covariates <- prop_vectors_2b$subinf_1
  P_Case_subinfect2__Covariates <- prop_vectors_2b$subinf_2

  obs_outcome <- ifelse(is.na(data[[laz_var_name]]), 0, data[[laz_var_name]])
  I_Delta_0 <- as.numeric(!is.na(data[[laz_var_name]]))

  for(i in seq_along(abx_levels)){

    abx_level <- abx_levels[i]

    I_Abx_a <- as.numeric(data[[abx_var_name]] == abx_level)
    I_Abx_a[is.na(I_Abx_a)] <- 0

    P_Abx_a__Case_Covariates <- prop_vectors_1a[, i]

    P_Delta_0__Case_all <- 1 - prop_vectors_3a[, i]

    if(!is.null(first_id_var_name)){
      pseudo_n <- mean(P_Delta_0__Case_all) * P_Case * length(unique(data[[first_id_var_name]]))
    } else{
      pseudo_n <- mean(P_Delta_0__Case_all) * P_Case * nrow(data)
    }

    truncate_factor <- 5 / (sqrt(pseudo_n) * log(pseudo_n))

    full_propensity <- P_Abx_a__Case_Covariates * P_Delta_0__Case_all
    full_propensity[full_propensity < truncate_factor] <- truncate_factor

    P_Case_subinfect1_tmp <- P_Case_subinfect1__Covariates
    P_Case_subinfect2_tmp <- P_Case_subinfect2__Covariates

    P_Case_subinfect1_tmp[P_Case_subinfect1_tmp < truncate_factor] <- truncate_factor
    P_Case_subinfect2_tmp[P_Case_subinfect2_tmp < truncate_factor] <- truncate_factor

    Qbar_Case_subinfect1_Abx_a_Covariates <- list_outcome_vectors_1a[[1]][, i]
    Qbar_Case_subinfect2_Abx_a_Covariates <- list_outcome_vectors_1a[[2]][, i]

    Qbar_Case_subinfect1_Covariates <- list_outcome_vectors_2a[[1]][, i]
    Qbar_Case_subinfect2_Covariates <- list_outcome_vectors_2a[[2]][, i]

    eif_case_subinfect1_vec <-
      (I_Case_subinfect1 / P_Case_subinfect1_tmp) *
      (P_Case__Covariates / P_Case) *
      (I_Abx_a * I_Delta_0) / full_propensity *
      (obs_outcome - Qbar_Case_subinfect1_Abx_a_Covariates) +
      (I_Case_subinfect1 / P_Case_subinfect1_tmp) *
      (P_Case__Covariates / P_Case) *
      (Qbar_Case_subinfect1_Abx_a_Covariates - Qbar_Case_subinfect1_Covariates) +
      (I_Case / P_Case) *
      (Qbar_Case_subinfect1_Covariates - plug_ins_case_subinfect1[i])

    eif_case_subinfect2_vec <-
      (I_Case_subinfect2 / P_Case_subinfect2_tmp) *
      (P_Case__Covariates / P_Case) *
      (I_Abx_a * I_Delta_0) / full_propensity *
      (obs_outcome - Qbar_Case_subinfect2_Abx_a_Covariates) +
      (I_Case_subinfect2 / P_Case_subinfect2_tmp) *
      (P_Case__Covariates / P_Case) *
      (Qbar_Case_subinfect2_Abx_a_Covariates - Qbar_Case_subinfect2_Covariates) +
      (I_Case / P_Case) *
      (Qbar_Case_subinfect2_Covariates - plug_ins_case_subinfect2[i])

    case_subinfect1_eifs[, i] <- ifelse(is.nan(eif_case_subinfect1_vec), 0, eif_case_subinfect1_vec)
    case_subinfect2_eifs[, i] <- ifelse(is.nan(eif_case_subinfect2_vec), 0, eif_case_subinfect2_vec)
  }

  # Control EIF, shared across antibiotic levels
  I_Control <- as.numeric(data[[case_var_name]] == 0)
  P_Control__Covariates <- 1 - prop_vectors_2a[, 1]

  P_Delta_0__Control_all <- 1 - prop_vectors_3b[, 1]

  Qbar_Control_Covariates <- outcome_vectors_1b[, 1]

  eif_control_vec <-
    (I_Control / P_Control__Covariates) *
    (I_Delta_0 / P_Delta_0__Control_all) *
    (P_Case__Covariates / P_Case) *
    (obs_outcome - Qbar_Control_Covariates) +
    (I_Case / P_Case) *
    (Qbar_Control_Covariates - plug_ins_control)

  control_eifs[, 1] <- ifelse(is.nan(eif_control_vec), 0, eif_control_vec)

  # Add bias correction
  aipw_case_subinfect1 <- plug_ins_case_subinfect1 + colMeans(case_subinfect1_eifs)
  aipw_case_subinfect2 <- plug_ins_case_subinfect2 + colMeans(case_subinfect2_eifs)
  aipw_control <- plug_ins_control + mean(control_eifs[, 1])

  eif_matrix <- cbind(case_subinfect1_eifs, case_subinfect2_eifs, control_eifs)

  # Get id for each participant and recreate EIFs based on this if present
  if(!is.null(first_id_var_name)){

    first_id_eif_matrix <- cbind(data.frame(first_id = data[[first_id_var_name]]), eif_matrix)
    first_id_eif_matrix <- aggregate(. ~ first_id, data = first_id_eif_matrix, FUN = sum)

    colnames(first_id_eif_matrix)[-1] <- colnames(eif_matrix)

    eif_first_ids <- first_id_eif_matrix[, 1]

    scaled_matrix <- first_id_eif_matrix[, -1, drop = FALSE] *
      (nrow(first_id_eif_matrix) / nrow(eif_matrix))

  } else{

    eif_first_ids <- NULL
    scaled_matrix <- eif_matrix
  }

  aipws_effect_subinfect1 <- vector("numeric", length = length(abx_levels))
  aipws_effect_subinfect2 <- vector("numeric", length = length(abx_levels))

  eifs_effect_subinfect1 <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(scaled_matrix)))
  eifs_effect_subinfect2 <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(scaled_matrix)))

  names(aipws_effect_subinfect1) <- paste0("effect_", abx_levels)
  names(aipws_effect_subinfect2) <- paste0("effect_", abx_levels)

  colnames(eifs_effect_subinfect1) <- paste0("effect_subinfect1_", abx_levels)
  colnames(eifs_effect_subinfect2) <- paste0("effect_subinfect2_", abx_levels)

  for(i in seq_along(abx_levels)){

    aipws_effect_subinfect1[i] <- aipw_case_subinfect1[i] - aipw_control
    aipws_effect_subinfect2[i] <- aipw_case_subinfect2[i] - aipw_control

    idx_1_subinfect1 <- i
    idx_1_subinfect2 <- length(abx_levels) + i
    idx_2 <- 2 * length(abx_levels) + 1

    gradient_subinfect1 <- rep(0, 2 * length(abx_levels) + 1)
    gradient_subinfect1[idx_1_subinfect1] <- 1
    gradient_subinfect1[idx_2] <- -1
    gradient_subinfect1 <- matrix(gradient_subinfect1, ncol = 1)

    gradient_subinfect2 <- rep(0, 2 * length(abx_levels) + 1)
    gradient_subinfect2[idx_1_subinfect2] <- 1
    gradient_subinfect2[idx_2] <- -1
    gradient_subinfect2 <- matrix(gradient_subinfect2, ncol = 1)

    eif_effect_subinfect1 <- as.numeric(as.matrix(scaled_matrix) %*% gradient_subinfect1)
    eifs_effect_subinfect1[, i] <- eif_effect_subinfect1

    eif_effect_subinfect2 <- as.numeric(as.matrix(scaled_matrix) %*% gradient_subinfect2)
    eifs_effect_subinfect2[, i] <- eif_effect_subinfect2
  }

  results_df <- data.frame(
    abx_levels = abx_levels,
    abx_level_case_subinfect1 = aipw_case_subinfect1,
    abx_level_case_subinfect2 = aipw_case_subinfect2,
    abx_level_control = rep(aipw_control, length(abx_levels)),
    effect_case_subinfect1_abx_level = aipws_effect_subinfect1,
    effect_case_subinfect2_abx_level = aipws_effect_subinfect2
  )

  eif_matrix_scaled <- cbind(
    scaled_matrix,
    eifs_effect_subinfect1,
    eifs_effect_subinfect2
  )

  cov_matrix <- stats::cov(eif_matrix_scaled)
  eif_hat <- sqrt(diag(cov_matrix) / nrow(eif_matrix_scaled))

  results_object <- list(
    results_df = results_df,
    plug_ins_case_subinfect1 = plug_ins_case_subinfect1,
    plug_ins_case_subinfect2 = plug_ins_case_subinfect2,
    plug_ins_control = plug_ins_control,
    eif_matrix = eif_matrix_scaled,
    eif_first_ids = eif_first_ids,
    se = eif_hat
  )

  class(results_object) <- "aipw_case_control"

  return(list(
    results_object = results_object,
    aipw_models = NULL
  ))
}



aipw_other_diarrhea <- function(data,
                                laz_var_name,
                                abx_var_name,
                                infection_var_name,
                                subinfection_var_name,
                                site_var_name,
                                followup_var_names,
                                covariate_list,
                                pathogen_quantity_list = NULL,
                                pathogen_attributable_list = NULL,
                                no_etiology_var_name = NULL,
                                outcome_type = "gaussian",
                                sl.library.outcome = c("SL.glm"),
                                sl.library.treatment = c("SL.mean"),
                                sl.library.infection = c("SL.glm"),
                                sl.library.missingness = c("SL.glm"),
                                v_folds = 5,
                                return_models = FALSE,
                                first_id_var_name = NULL,
                                msm = FALSE,
                                msm_var_name = NULL,
                                msm_formula = NULL,
                                ps_trunc_level = 0.01,
                                all_other_diarrhea = FALSE){

  if(msm){
    stop("msm not done yet for subtype-specific aipw_other_diarrhea")
  }

  if(return_models){
    stop("return_models not done yet for subtype-specific aipw_other_diarrhea")
  }

  # ------------------------------------------------------------
  # STEP 0: Create subsets of data for model fitting
  # ------------------------------------------------------------

  abx_levels <- levels(factor(data[[abx_var_name]]))
  abx_levels_new <- ifelse(is.na(abx_levels), NA, gsub("[ /]", "_", abx_levels))
  data[[abx_var_name]] <- factor(data[[abx_var_name]], levels = abx_levels, labels = abx_levels_new)

  inf_attr_idx <- which(data[[infection_var_name]] == 1)
  sub_inf_attr <- data[inf_attr_idx, ]

  subinfection_var_levels <- unique(sub_inf_attr[[subinfection_var_name]])
  subinfection_var_levels <- subinfection_var_levels[!is.na(subinfection_var_levels)]

  if(length(subinfection_var_levels) != 2){
    stop("This function currently assumes exactly two non-missing subinfection levels among infected observations.")
  }

  I_Y_inf_attr <- ifelse(is.na(sub_inf_attr[[laz_var_name]]), 1, 0)
  Y_inf_attr <- sub_inf_attr[[laz_var_name]]
  attributes(Y_inf_attr) <- NULL

  covariates_inf_attr <- sub_inf_attr[, covariate_list, drop = FALSE]
  abx_inf_attr <- sub_inf_attr[, abx_var_name, drop = FALSE]
  subtype_inf_attr <- sub_inf_attr[, subinfection_var_name, drop = FALSE]

  sub_inf_attr_complete <- sub_inf_attr[!is.na(sub_inf_attr[[laz_var_name]]), ]

  Y_inf_attr_complete <- sub_inf_attr_complete[[laz_var_name]]
  attributes(Y_inf_attr_complete) <- NULL

  covariates_inf_attr_complete <- sub_inf_attr_complete[, covariate_list, drop = FALSE]
  abx_inf_attr_complete <- sub_inf_attr_complete[, abx_var_name, drop = FALSE]
  subtype_inf_attr_complete <- sub_inf_attr_complete[, subinfection_var_name, drop = FALSE]

  if(!is.null(no_etiology_var_name)){
    sub_no_attr <- data[which(data[[no_etiology_var_name]] == 1), ]
  } else{
    I_no_attr <- ifelse(
      data[[infection_var_name]] == 0 &
        rowSums(data[, pathogen_attributable_list, drop = FALSE], na.rm = TRUE) == 0,
      1, 0
    )
    data$no_etiology <- I_no_attr
    no_etiology_var_name <- "no_etiology"

    sub_no_attr <- data[which(data$no_etiology == 1), ]
  }

  I_Y_no_attr <- ifelse(is.na(sub_no_attr[[laz_var_name]]), 1, 0)
  Y_no_attr <- sub_no_attr[[laz_var_name]]
  attributes(Y_no_attr) <- NULL

  covariates_no_attr <- sub_no_attr[, covariate_list, drop = FALSE]
  abx_no_attr <- sub_no_attr[, abx_var_name, drop = FALSE]

  sub_no_attr_complete <- sub_no_attr[!is.na(sub_no_attr[[laz_var_name]]), ]

  Y_no_attr_complete <- sub_no_attr_complete[[laz_var_name]]
  attributes(Y_no_attr_complete) <- NULL

  covariates_no_attr_complete <- sub_no_attr_complete[, covariate_list, drop = FALSE]
  abx_no_attr_complete <- sub_no_attr_complete[, abx_var_name, drop = FALSE]

  # ------------------------------------------------------------
  # STEP 1: Fit & predict from outcome models
  # ------------------------------------------------------------

  outcome_model_1a <- SuperLearner::SuperLearner(
    Y = Y_inf_attr_complete,
    X = data.frame(
      abx_inf_attr_complete,
      covariates_inf_attr_complete,
      subtype_inf_attr_complete
    ),
    family = outcome_type,
    SL.library = sl.library.outcome,
    cvControl = list(V = v_folds)
  )

  outcome_model_1b <- SuperLearner::SuperLearner(
    Y = Y_no_attr_complete,
    X = data.frame(
      abx_no_attr_complete,
      covariates_no_attr_complete
    ),
    family = outcome_type,
    SL.library = sl.library.outcome,
    cvControl = list(V = v_folds)
  )

  abx_levels <- unique(data[[abx_var_name]])
  abx_levels <- abx_levels[!is.na(abx_levels)]

  list_outcome_vectors_1a <- vector("list", length = length(subinfection_var_levels))
  outcome_vectors_1b <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))

  colnames(outcome_vectors_1b) <- paste0("abx_", abx_levels)

  for(s in seq_along(subinfection_var_levels)){

    subinfect_level <- subinfection_var_levels[s]

    list_outcome_vectors_1a[[s]] <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
    colnames(list_outcome_vectors_1a[[s]]) <- paste0("abx_", abx_levels)

    for(i in seq_along(abx_levels)){

      abx_level <- abx_levels[i]

      data_abx <- data
      data_abx[[abx_var_name]] <- abx_level

      pred_data <- data_abx[, c(abx_var_name, covariate_list), drop = FALSE]
      pred_data[[subinfection_var_name]] <- subinfect_level

      list_outcome_vectors_1a[[s]][, i] <- stats::predict(
        outcome_model_1a,
        newdata = pred_data,
        type = "response"
      )$pred
    }
  }

  for(i in seq_along(abx_levels)){

    abx_level <- abx_levels[i]

    data_abx <- data
    data_abx[[abx_var_name]] <- abx_level

    outcome_vectors_1b[, i] <- stats::predict(
      outcome_model_1b,
      newdata = data_abx[, c(abx_var_name, covariate_list), drop = FALSE],
      type = "response"
    )$pred
  }

  # ------------------------------------------------------------
  # STEP 2: Fit & predict from propensity models
  # ------------------------------------------------------------

  if(!any(is.na(followup_var_names))){
    covariate_list_no_followup <- covariate_list[!(covariate_list %in% followup_var_names)]
  } else{
    covariate_list_no_followup <- covariate_list
  }

  prop_vectors_1a <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  prop_vectors_1b <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  prop_vectors_2a <- data.frame(matrix(ncol = 2, nrow = nrow(data)))
  prop_vectors_2b <- data.frame(matrix(ncol = 2, nrow = nrow(data)))
  prop_vectors_3a <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  prop_vectors_3b <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))

  colnames(prop_vectors_1a) <- paste0("abx_", abx_levels)
  colnames(prop_vectors_1b) <- paste0("abx_", abx_levels)
  colnames(prop_vectors_2a) <- c("inf_attr", "no_attr")
  colnames(prop_vectors_2b) <- c("subinf_1", "subinf_2")
  colnames(prop_vectors_3a) <- paste0("abx_", abx_levels)
  colnames(prop_vectors_3b) <- paste0("abx_", abx_levels)

  ###############################################
  ## Part 1: Propensity models for antibiotics ##
  ###############################################

  for(i in seq_along(abx_levels)){

    abx_level <- abx_levels[i]

    if(i != length(abx_levels)){

      if(i > 1){
        abx_levels_out <- abx_levels[1:(i - 1)]
        prop_sub_inf_attr <- sub_inf_attr[-which(sub_inf_attr[[abx_var_name]] %in% abx_levels_out), ]
        prop_sub_no_attr <- sub_no_attr[-which(sub_no_attr[[abx_var_name]] %in% abx_levels_out), ]
      } else{
        prop_sub_inf_attr <- sub_inf_attr
        prop_sub_no_attr <- sub_no_attr
      }

      prop_covariates_inf_attr <- prop_sub_inf_attr[, covariate_list, drop = FALSE]
      prop_subtype_inf_attr <- prop_sub_inf_attr[, subinfection_var_name, drop = FALSE]

      prop_covariates_no_attr <- prop_sub_no_attr[, covariate_list, drop = FALSE]

      prop_model_1a <- SuperLearner::SuperLearner(
        Y = as.numeric(prop_sub_inf_attr[[abx_var_name]] == abx_level),
        X = data.frame(
          prop_covariates_inf_attr,
          prop_subtype_inf_attr
        ),
        newX = data[, c(covariate_list, subinfection_var_name), drop = FALSE],
        family = stats::binomial(),
        SL.library = sl.library.treatment,
        cvControl = list(V = v_folds)
      )

      tmp_pred_a <- prop_model_1a$SL.pred

      prop_model_1b <- SuperLearner::SuperLearner(
        Y = as.numeric(prop_sub_no_attr[[abx_var_name]] == abx_level),
        X = data.frame(prop_covariates_no_attr),
        newX = data[, covariate_list, drop = FALSE],
        family = stats::binomial(),
        SL.library = sl.library.treatment,
        cvControl = list(V = v_folds)
      )

      tmp_pred_b <- prop_model_1b$SL.pred

      if(i == 1){
        prop_vectors_1a[, i] <- tmp_pred_a
        prop_vectors_1b[, i] <- tmp_pred_b
      } else{
        for(j in 1:(i - 1)){
          tmp_pred_a <- tmp_pred_a * (1 - prop_vectors_1a[, j])
          tmp_pred_b <- tmp_pred_b * (1 - prop_vectors_1b[, j])
        }
        prop_vectors_1a[, i] <- tmp_pred_a
        prop_vectors_1b[, i] <- tmp_pred_b
      }

    } else{
      prop_vectors_1a[, i] <- 1 - rowSums(prop_vectors_1a[, 1:(ncol(prop_vectors_1a) - 1), drop = FALSE])
      prop_vectors_1b[, i] <- 1 - rowSums(prop_vectors_1b[, 1:(ncol(prop_vectors_1b) - 1), drop = FALSE])
    }
  }

  ##########################################################
  ## Part 2: Propensity models for infection and subtype ##
  ##########################################################

  prop_model_2a_1 <- SuperLearner::SuperLearner(
    Y = data[[infection_var_name]],
    X = data[, covariate_list, drop = FALSE],
    family = stats::binomial(),
    SL.library = sl.library.infection,
    cvControl = list(V = v_folds)
  )

  prop_vectors_2a$inf_attr <- prop_model_2a_1$SL.pred

  prop_model_2b_1 <- SuperLearner::SuperLearner(
    Y = as.numeric(sub_inf_attr[[subinfection_var_name]] == subinfection_var_levels[1]),
    X = sub_inf_attr[, covariate_list, drop = FALSE],
    newX = data[, covariate_list, drop = FALSE],
    family = stats::binomial(),
    SL.library = sl.library.infection,
    cvControl = list(V = v_folds)
  )

  P_subinfect_is_level_1_given_attr <- prop_model_2b_1$SL.pred
  P_subinfect_is_level_1_and_attr <- P_subinfect_is_level_1_given_attr * prop_vectors_2a$inf_attr

  P_subinfect_is_level_2_given_attr <- 1 - P_subinfect_is_level_1_given_attr
  P_subinfect_is_level_2_and_attr <- P_subinfect_is_level_2_given_attr * prop_vectors_2a$inf_attr

  prop_vectors_2b$subinf_1 <- P_subinfect_is_level_1_and_attr
  prop_vectors_2b$subinf_2 <- P_subinfect_is_level_2_and_attr

  sub_no_shig <- data[which(data[[infection_var_name]] == 0), ]

  if(is.na(no_etiology_var_name)){
    sub_no_shig$no_etiology <- ifelse(
      rowSums(sub_no_shig[, pathogen_attributable_list, drop = FALSE], na.rm = TRUE) == 0,
      1, 0
    )
    no_etiology_var_name <- "no_etiology"
  }

  if(!all_other_diarrhea){

    prop_model_2a_2 <- SuperLearner::SuperLearner(
      Y = sub_no_shig[[no_etiology_var_name]],
      X = sub_no_shig[, covariate_list, drop = FALSE],
      newX = data[, covariate_list, drop = FALSE],
      family = stats::binomial(),
      SL.library = sl.library.infection,
      cvControl = list(V = v_folds)
    )

    prop_vectors_2a$no_attr <- prop_model_2a_2$SL.pred * (1 - prop_vectors_2a$inf_attr)

  } else{

    prop_model_2a_2 <- NULL
    prop_vectors_2a$no_attr <- 1 - prop_vectors_2a$inf_attr
  }

  ###############################################
  ## Part 3: Propensity models for missingness ##
  ###############################################

  prop_covariates_inf_attr <- covariates_inf_attr[
    , colnames(covariates_inf_attr) %in% covariate_list_no_followup,
    drop = FALSE
  ]

  prop_covariates_no_attr <- covariates_no_attr[
    , colnames(covariates_no_attr) %in% covariate_list_no_followup,
    drop = FALSE
  ]

  if(sum(I_Y_inf_attr) == 0){

    prop_model_3a <- list()
    class(prop_model_3a) <- "constant_zero_model"

    predict.constant_zero_model <- function(object, newdata, ...) {
      list(pred = rep(0, nrow(newdata)))
    }

  } else{

    prop_model_3a <- SuperLearner::SuperLearner(
      Y = I_Y_inf_attr,
      X = data.frame(
        abx_inf_attr,
        prop_covariates_inf_attr,
        subtype_inf_attr
      ),
      family = stats::binomial(),
      SL.library = sl.library.missingness,
      cvControl = list(V = v_folds)
    )
  }

  if(sum(I_Y_no_attr) == 0){

    prop_model_3b <- list()
    class(prop_model_3b) <- "constant_zero_model"

    predict.constant_zero_model <- function(object, newdata, ...) {
      list(pred = rep(0, nrow(newdata)))
    }

  } else{

    prop_model_3b <- SuperLearner::SuperLearner(
      Y = I_Y_no_attr,
      X = data.frame(
        abx_no_attr,
        prop_covariates_no_attr
      ),
      family = stats::binomial(),
      SL.library = sl.library.missingness,
      cvControl = list(V = v_folds)
    )
  }

  for(i in seq_along(abx_levels)){

    abx_level <- abx_levels[i]

    pred_data <- data
    pred_data[[abx_var_name]] <- abx_level

    prop_vectors_3a[, i] <- stats::predict(
      prop_model_3a,
      newdata = pred_data[, c(
        abx_var_name,
        covariate_list_no_followup,
        subinfection_var_name
      ), drop = FALSE],
      type = "response"
    )$pred

    prop_vectors_3b[, i] <- stats::predict(
      prop_model_3b,
      newdata = pred_data[, c(
        abx_var_name,
        covariate_list_no_followup
      ), drop = FALSE],
      type = "response"
    )$pred
  }

  # -------------------------------------------------
  # STEP 3: AIPW estimates and confidence intervals
  # -------------------------------------------------

  plug_ins_inf_subinfect1 <- colMeans(list_outcome_vectors_1a[[1]][inf_attr_idx, , drop = FALSE])
  plug_ins_inf_subinfect2 <- colMeans(list_outcome_vectors_1a[[2]][inf_attr_idx, , drop = FALSE])
  plug_ins_no_attr <- colMeans(outcome_vectors_1b[inf_attr_idx, , drop = FALSE])

  inf_subinfect1_eifs <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  inf_subinfect2_eifs <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  no_attr_eifs <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))

  colnames(inf_subinfect1_eifs) <- paste0("inf_subinfect1_eif_", abx_levels)
  colnames(inf_subinfect2_eifs) <- paste0("inf_subinfect2_eif_", abx_levels)
  colnames(no_attr_eifs) <- paste0("no_attr_eif_", abx_levels)

  if(!is.na(ps_trunc_level)){

    truncate_ps <- function(mat, ps_trunc_level){
      mat[mat < ps_trunc_level] <- ps_trunc_level
      mat[mat > 1 - ps_trunc_level] <- 1 - ps_trunc_level
      return(mat)
    }

    prop_vectors_1a <- truncate_ps(prop_vectors_1a, ps_trunc_level)
    prop_vectors_1b <- truncate_ps(prop_vectors_1b, ps_trunc_level)
    prop_vectors_2a <- truncate_ps(prop_vectors_2a, ps_trunc_level)
    prop_vectors_2b <- truncate_ps(prop_vectors_2b, ps_trunc_level)
    prop_vectors_3a <- truncate_ps(prop_vectors_3a, ps_trunc_level)
    prop_vectors_3b <- truncate_ps(prop_vectors_3b, ps_trunc_level)
  }

  I_Inf_1 <- data[[infection_var_name]]
  P_Inf_1 <- mean(prop_vectors_2a[, 1])
  P_Inf_1__Covariates <- prop_vectors_2a[, 1]

  I_Inf_subinfect1_1 <- as.numeric(data[[subinfection_var_name]] == subinfection_var_levels[1])
  I_Inf_subinfect2_1 <- as.numeric(data[[subinfection_var_name]] == subinfection_var_levels[2])

  P_Inf_subinfect1_1 <- prop_vectors_2b$subinf_1
  P_Inf_subinfect2_1 <- prop_vectors_2b$subinf_2

  I_Delta_0 <- as.numeric(!is.na(data[[laz_var_name]]))
  obs_outcome <- ifelse(is.na(data[[laz_var_name]]), 0, data[[laz_var_name]])

  for(i in seq_along(abx_levels)){

    abx_level <- abx_levels[i]

    I_Abx_a <- as.numeric(data[[abx_var_name]] == abx_level)
    P_Abx_a__Inf_1_Covariates <- prop_vectors_1a[, i]

    P_Delta_0__Inf_all <- 1 - prop_vectors_3a[, i]

    if(!is.null(first_id_var_name)){
      pseudo_n <- mean(P_Delta_0__Inf_all) * P_Inf_1 * length(unique(data[[first_id_var_name]]))
    } else{
      pseudo_n <- mean(P_Delta_0__Inf_all) * P_Inf_1 * nrow(data)
    }

    truncate_factor <- 5 / (sqrt(pseudo_n) * log(pseudo_n))

    full_propensity <- P_Abx_a__Inf_1_Covariates * P_Delta_0__Inf_all
    full_propensity[full_propensity < truncate_factor] <- truncate_factor

    P_Inf_subinfect1_tmp <- P_Inf_subinfect1_1
    P_Inf_subinfect2_tmp <- P_Inf_subinfect2_1

    P_Inf_subinfect1_tmp[P_Inf_subinfect1_tmp < truncate_factor] <- truncate_factor
    P_Inf_subinfect2_tmp[P_Inf_subinfect2_tmp < truncate_factor] <- truncate_factor

    Qbar_Inf_1_subinfect1_Abx_a_Covariates <- list_outcome_vectors_1a[[1]][, i]
    Qbar_Inf_1_subinfect2_Abx_a_Covariates <- list_outcome_vectors_1a[[2]][, i]

    eif_vec_inf_subinfect1 <-
      (I_Inf_subinfect1_1 / P_Inf_subinfect1_tmp) *
      (P_Inf_1__Covariates / P_Inf_1) *
      (I_Abx_a * I_Delta_0) / full_propensity *
      (obs_outcome - Qbar_Inf_1_subinfect1_Abx_a_Covariates) +
      (I_Inf_1 / P_Inf_1) *
      (Qbar_Inf_1_subinfect1_Abx_a_Covariates - plug_ins_inf_subinfect1[i])

    eif_vec_inf_subinfect2 <-
      (I_Inf_subinfect2_1 / P_Inf_subinfect2_tmp) *
      (P_Inf_1__Covariates / P_Inf_1) *
      (I_Abx_a * I_Delta_0) / full_propensity *
      (obs_outcome - Qbar_Inf_1_subinfect2_Abx_a_Covariates) +
      (I_Inf_1 / P_Inf_1) *
      (Qbar_Inf_1_subinfect2_Abx_a_Covariates - plug_ins_inf_subinfect2[i])

    if(is.na(no_etiology_var_name)){
      I_No_attr_1 <- data$no_etiology
    } else{
      I_No_attr_1 <- data[[no_etiology_var_name]]
    }

    P_No_attr__Covariates <- prop_vectors_2a[, 2]
    P_Abx_a__Covariates <- prop_vectors_1b[, i]
    P_Delta_0__No_attr_all <- 1 - prop_vectors_3b[, i]

    Qbar_No_attr_Abx_a_Covariates <- outcome_vectors_1b[, i]

    if(!is.null(first_id_var_name)){
      pseudo_n <- mean(P_Delta_0__No_attr_all) * P_Inf_1 * length(unique(data[[first_id_var_name]]))
    } else{
      pseudo_n <- mean(P_Delta_0__No_attr_all) * P_Inf_1 * nrow(data)
    }

    truncate_factor <- 5 / (sqrt(pseudo_n) * log(pseudo_n))

    full_propensity <- P_Abx_a__Covariates * P_Delta_0__No_attr_all
    full_propensity[full_propensity < truncate_factor] <- truncate_factor

    P_No_attr_tmp <- P_No_attr__Covariates
    P_No_attr_tmp[P_No_attr_tmp < truncate_factor] <- truncate_factor

    eif_vec_no_attr <-
      (I_No_attr_1 / P_No_attr_tmp) *
      (P_Inf_1__Covariates / P_Inf_1) *
      (I_Abx_a * I_Delta_0) / full_propensity *
      (obs_outcome - Qbar_No_attr_Abx_a_Covariates) +
      (I_Inf_1 / P_Inf_1) *
      (Qbar_No_attr_Abx_a_Covariates - plug_ins_no_attr[i])

    inf_subinfect1_eifs[, i] <- eif_vec_inf_subinfect1
    inf_subinfect2_eifs[, i] <- eif_vec_inf_subinfect2
    no_attr_eifs[, i] <- eif_vec_no_attr
  }

  aipw_inf_subinfect1 <- plug_ins_inf_subinfect1 + colMeans(inf_subinfect1_eifs)
  aipw_inf_subinfect2 <- plug_ins_inf_subinfect2 + colMeans(inf_subinfect2_eifs)
  aipw_no_attr <- plug_ins_no_attr + colMeans(no_attr_eifs)

  eif_matrix <- cbind(inf_subinfect1_eifs, inf_subinfect2_eifs, no_attr_eifs)

  if(!is.null(first_id_var_name)){

    first_id_eif_matrix <- cbind(data.frame(first_id = data[[first_id_var_name]]), eif_matrix)
    first_id_eif_matrix <- aggregate(. ~ first_id, data = first_id_eif_matrix, FUN = sum)

    colnames(first_id_eif_matrix)[-1] <- colnames(eif_matrix)

    eif_first_ids <- first_id_eif_matrix[, 1]

    scaled_matrix <- first_id_eif_matrix[, -1, drop = FALSE] *
      (nrow(first_id_eif_matrix) / nrow(eif_matrix))

  } else{

    eif_first_ids <- NULL
    scaled_matrix <- eif_matrix
  }

  aipws_effect_subinfect1 <- vector("numeric", length = length(abx_levels))
  aipws_effect_subinfect2 <- vector("numeric", length = length(abx_levels))

  eifs_effect_subinfect1 <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(scaled_matrix)))
  eifs_effect_subinfect2 <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(scaled_matrix)))

  names(aipws_effect_subinfect1) <- paste0("effect_", abx_levels)
  names(aipws_effect_subinfect2) <- paste0("effect_", abx_levels)

  colnames(eifs_effect_subinfect1) <- paste0("effect_subinfect1_", abx_levels)
  colnames(eifs_effect_subinfect2) <- paste0("effect_subinfect2_", abx_levels)

  for(i in seq_along(abx_levels)){

    aipws_effect_subinfect1[i] <- aipw_inf_subinfect1[i] - aipw_no_attr[i]
    aipws_effect_subinfect2[i] <- aipw_inf_subinfect2[i] - aipw_no_attr[i]

    idx_1_subinfect1 <- i
    idx_1_subinfect2 <- length(abx_levels) + i
    idx_2 <- 2 * length(abx_levels) + i

    gradient_subinfect1 <- rep(0, 3 * length(abx_levels))
    gradient_subinfect1[idx_1_subinfect1] <- 1
    gradient_subinfect1[idx_2] <- -1
    gradient_subinfect1 <- matrix(gradient_subinfect1, ncol = 1)

    gradient_subinfect2 <- rep(0, 3 * length(abx_levels))
    gradient_subinfect2[idx_1_subinfect2] <- 1
    gradient_subinfect2[idx_2] <- -1
    gradient_subinfect2 <- matrix(gradient_subinfect2, ncol = 1)

    eif_effect_subinfect1 <- as.numeric(as.matrix(scaled_matrix) %*% gradient_subinfect1)
    eifs_effect_subinfect1[, i] <- eif_effect_subinfect1

    eif_effect_subinfect2 <- as.numeric(as.matrix(scaled_matrix) %*% gradient_subinfect2)
    eifs_effect_subinfect2[, i] <- eif_effect_subinfect2
  }

  results_df <- data.frame(
    abx_levels = abx_levels,
    abx_level_inf_subinfect1_1 = aipw_inf_subinfect1,
    abx_level_inf_subinfect2_1 = aipw_inf_subinfect2,
    abx_level_inf_0 = aipw_no_attr,
    effect_inf_subinfect1_abx_level = aipws_effect_subinfect1,
    effect_inf_subinfect2_abx_level = aipws_effect_subinfect2
  )

  eif_matrix_scaled <- cbind(
    scaled_matrix,
    eifs_effect_subinfect1,
    eifs_effect_subinfect2
  )

  cov_matrix <- stats::cov(eif_matrix_scaled)
  eif_hat <- sqrt(diag(cov_matrix) / nrow(eif_matrix_scaled))

  results_object <- list(
    results_df = results_df,
    plug_ins_inf_subinfect1 = plug_ins_inf_subinfect1,
    plug_ins_inf_subinfect2 = plug_ins_inf_subinfect2,
    plug_ins_no_attr = plug_ins_no_attr,
    eif_matrix = eif_matrix_scaled,
    eif_first_ids = eif_first_ids,
    se = eif_hat
  )

  class(results_object) <- "aipw_other_diarrhea"

  return(list(
    results_object = results_object,
    aipw_models = NULL
  ))
}