#' Function for AIPW estimates of effect of infection on growth- single outcome regression
#' 
#' @param data dataframe containing dataset used for gcomp
#' @param laz_var_name name of growth outcome variable
#' @param abx_var_name name of binary antibiotic variable
#' @param infection_var_name name of binary infection variable (for other diarrhea analyses; case_control = FALSE)
#' @param site_var_name name of covariate for site (to exclude from propensity models)
#' @param covariate_list character vector containing names of baseline covariates
#' @param pathogen_quantity_list character vector containing name(s) of pathogen quantity variables in dataset
#' @param pathogen_attributable_list character vector containing name(s) of binary pathogen attributable variables in dataset. Used to identify diarrhea with no etiology if is.null(no_etiology_var_name)
#' @param no_etiology_var_name name of binary variable indicating diarrhea episodes with no known etiology. If null, use pathogen_attributable_list
#' @param outcome_type gaussian or binomial for continuous or binomial outcome
#' @param sl.library.outcome list containing SuperLearner libraries to use for outcome regression 
#' @param sl.library.treatment list containing SuperLearner libraries to use for antibiotic prescribing propensity model
#' @param sl.library.infection list containing SuperLearner libraries to use for infection propensity model
#' @param sl.library.missingness list containing SuperLearner libraries to use for outcome missingness model
#' @param v_folds number of cross-validation folds to use in SuperLearner
#' @param return_models boolean return SuperLearner models. Default FALSE.
#' @param first_id_var_name name of variable indicating first_id in data. Used to account for re-enrollment in study. 
#' @param msm boolean indicating use of MSM for effect heterogeneity, default FALSE
#' @param msm_var_name name of variable to use for msm
#' @param msm_formula chatacter vector with formula to use for msm if msm TRUE
#' 
#' @keywords internal
#' 
#' @returns List containing `aipw_other_diarrhea` object, models (if return_models = TRUE). `aipw_other_diarrhea` object contains dataframe with results, covariance matrix, standard errors.
aipw_other_diarrhea <- function(data,
                                laz_var_name,
                                abx_var_name,
                                infection_var_name,
                                site_var_name,
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
                                msm_formula = NULL){
  
  # ------------------------------------------------------------
  # STEP 0: Create subsets of data for model fitting
  # ------------------------------------------------------------
  
  # Rename abx levels if spaces in them
  abx_levels <- levels(factor(data[[abx_var_name]]))
  abx_levels_new <- ifelse(is.na(abx_levels), NA, gsub("[ /]", "_", abx_levels))
  data[[abx_var_name]] <- factor(data[[abx_var_name]], levels = abx_levels, labels = abx_levels_new)
  
  ## Subset 0a: subset shigella (or other infection variable) cases
  inf_attr_idx <- which(data[[infection_var_name]] == 1)
  sub_inf_attr <- data[inf_attr_idx,]
  
  # Full data (including missing outcome)
  I_Y_inf_attr <- ifelse(is.na(sub_inf_attr[[laz_var_name]]), 1, 0) #indicator for Y missing
  Y_inf_attr <- sub_inf_attr[[laz_var_name]]
  covariates_inf_attr <- sub_inf_attr[, covariate_list, drop = FALSE]
  abx_inf_attr <- sub_inf_attr[, abx_var_name, drop = FALSE]
  
  # Complete data (excluding missing outcome)
  sub_inf_attr_complete <- sub_inf_attr[!is.na(sub_inf_attr[[laz_var_name]]), ]
  
  Y_inf_attr_complete <- sub_inf_attr_complete[[laz_var_name]]
  covariates_inf_attr_complete <- sub_inf_attr_complete[, covariate_list, drop = FALSE]
  abx_inf_attr_complete <- sub_inf_attr_complete[, abx_var_name, drop = FALSE]
  
  ## Subset 0b: subset to cases with no etiology
  if(!is.null(no_etiology_var_name)){
    sub_no_attr <- data[which(data[[no_etiology_var_name]] == 1),]
  } else{
    I_no_attr <- ifelse(data[[infection_var_name]] == 0 & (rowSums(data[,pathogen_attributable_list]) == 0),
                        1, 0)
    data$no_etiology <- I_no_attr
    no_etiology_var_name <- "no_etiology"
    
    sub_no_attr <- data[which(data$no_etiology == 1),]
  }
  
  I_Y_no_attr <- ifelse(is.na(sub_no_attr[[laz_var_name]]), 1, 0) #indicator for Y missing
  Y_no_attr <- sub_no_attr[[laz_var_name]]
  covariates_no_attr <- sub_no_attr[, covariate_list, drop = FALSE]
  abx_no_attr <- sub_no_attr[, abx_var_name, drop = FALSE]
  
  sub_no_attr_complete <- sub_no_attr[!is.na(sub_no_attr[[laz_var_name]]), ]
  Y_no_attr_complete <- sub_no_attr_complete[[laz_var_name]]
  covariates_no_attr_complete <- sub_no_attr_complete[, covariate_list, drop = FALSE]
  abx_no_attr_complete <- sub_no_attr_complete[, abx_var_name, drop = FALSE]
  
  # ------------------------------------------------------------
  # STEP 1: Fit & predict from outcome models
  # ------------------------------------------------------------
  
  ## Model 1a: Outcome model in Shigella (or other infection) attributable cases
  outcome_model_1a <- SuperLearner::SuperLearner(Y = Y_inf_attr_complete, 
                                                 X = data.frame(abx_inf_attr_complete,
                                                                covariates_inf_attr_complete),
                                                 family = outcome_type, 
                                                 SL.library = sl.library.outcome,
                                                 cvControl = list(V = v_folds))
  
  ## Model 1b: Outcome model in cases with no attribution
  outcome_model_1b <- SuperLearner::SuperLearner(Y = Y_no_attr_complete, 
                                                 X = data.frame(abx_no_attr_complete,
                                                                covariates_no_attr_complete),
                                                 family = outcome_type, 
                                                 SL.library = sl.library.outcome,
                                                 cvControl = list(V = v_folds))
  
  
  ## Predict outcomes
  
  # For each level of abx, predict setting abx = x, infection = 1 & abx = x, infection = 0
  abx_levels <- unique(data[[abx_var_name]])[!is.na(unique(data[[abx_var_name]]))]
  
  outcome_vectors_1a <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  outcome_vectors_1b <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  
  colnames(outcome_vectors_1a) <- paste0("abx_", abx_levels)
  colnames(outcome_vectors_1b) <- paste0("abx_", abx_levels)
  
  # If looking at effect heterogeneity, create matrix to hold difference in outcome vectors
  if(msm){
    outcome_msm <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
    colnames(outcome_msm) <- paste0("abx_", abx_levels)
  }
  
  # Iterate through each abx level
  for (i in 1:length(abx_levels)) {
    abx_level <- abx_levels[i]
    
    data_abx <- data
    data_abx[[abx_var_name]] <- abx_level
    
    ## 2a: Predictions from 1a (cases)
    
    outcome_vectors_1a[, i] <- stats::predict(outcome_model_1a, newdata = data_abx[, c(abx_var_name,
                                                                                       covariate_list)], type = "response")$pred
    
    
    ## 2b: Predictions from 1b (controls)
    
    outcome_vectors_1b[, i] <- stats::predict(outcome_model_1b, newdata = data_abx[, c(abx_var_name,
                                                                                       covariate_list)], type = "response")$pred
  }
  
  if(msm){
    # Subtract predictions
    outcome_msm <- outcome_vectors_1a - outcome_vectors_1b
    
    msm_vectors <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
    colnames(msm_vectors) <- paste0("abx_", abx_levels)
    
    if(return_models){
      msm_model_list <- vector("list", length = length(abx_levels))
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
      
      effect_hetero_msm <- stats::glm(msm_formula_full, 
                                      data = msm_data[msm_data[[infection_var_name]] == 1,],       # QUESTION subset to shigella people here? then predict on everyone?
                                      family = outcome_type)                                       # QUESTION in general does this only work for continuous? because difference in outcome vectors?
      
      msm_vectors[,i] <- stats::predict(effect_hetero_msm, newdata = msm_data, type = 'response')
      
      if(return_models){
        msm_model_list[[i]] <- effect_hetero_msm
      }
      
    }
    
  }
  
  # ------------------------------------------------------------
  # STEP 2: Fit & predict from propensity models
  # ------------------------------------------------------------
  
  # If site in covariate_list, remove
  if(!any(is.na(site_var_name))){
    covariate_list_no_site <- covariate_list[!(covariate_list %in% site_var_name)]
  }
  
  prop_vectors_1a <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  prop_vectors_1b <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  prop_vectors_2a <- data.frame(matrix(ncol = 2, nrow = nrow(data)))
  prop_vectors_3a <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  prop_vectors_3b <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  
  colnames(prop_vectors_1a) <- paste0("abx_", abx_levels)
  colnames(prop_vectors_1b) <- paste0("abx_", abx_levels)
  colnames(prop_vectors_2a) <- c("inf_attr", "no_attr")
  colnames(prop_vectors_3a) <- paste0("abx_", abx_levels)
  colnames(prop_vectors_3b) <- paste0("abx_", abx_levels)
  
  if(return_models){
    prop_model_1a_list <- vector("list", length = length(abx_levels))
    prop_model_1b_list <- vector("list", length = length(abx_levels))
  }
  
  ## Part 1: Propensity models for antibiotics - RCT when using this function with no 2nd stage regression so pass in hard coded GLM
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
      
      prop_covariates_inf_attr <- prop_sub_inf_attr[, covariate_list_no_site , drop = FALSE]
      prop_pathogen_inf_attr <- prop_sub_inf_attr[, pathogen_quantity_list, drop = FALSE]
      
      prop_covariates_no_attr <- prop_sub_no_attr[, covariate_list_no_site, drop = FALSE]
      prop_pathogen_no_attr <- prop_sub_no_attr[, pathogen_quantity_list, drop = FALSE]
      
      ## 1a. Propensity model for abx shigella attributable
      prop_model_1a <- SuperLearner::SuperLearner(Y = as.numeric(prop_sub_inf_attr[[abx_var_name]] == abx_level),
                                                  X = data.frame(prop_covariates_inf_attr,
                                                                 prop_pathogen_inf_attr),
                                                  newX = data[,c(covariate_list_no_site ,
                                                                 pathogen_quantity_list)], 
                                                  family = stats::binomial(), 
                                                  SL.library = sl.library.treatment,
                                                  cvControl = list(V = v_folds))
      
      # Predictions from full data
      tmp_pred_a <- prop_model_1a$SL.pred
      
      ## 1b. Propensity model for abx no attribution
      prop_model_1b <- SuperLearner::SuperLearner(Y = as.numeric(prop_sub_no_attr[[abx_var_name]] == abx_level),
                                                  X = data.frame(prop_covariates_no_attr,
                                                                 prop_pathogen_no_attr),
                                                  newX = data[,c(covariate_list_no_site ,
                                                                 pathogen_quantity_list)], 
                                                  family = stats::binomial(), 
                                                  SL.library = sl.library.treatment,
                                                  cvControl = list(V = v_folds))
      
      # Predictions from full data
      tmp_pred_b <- prop_model_1b$SL.pred
      
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
  
  
  ## Part 2: Propensity models for shigella (or other infection) attribution
  
  # 2a_1 = Shigella Attributable ~ BL Cov
  prop_model_2a_1 <- SuperLearner::SuperLearner(Y = data[[infection_var_name]],
                                                X = data[, covariate_list_no_site , drop = FALSE], 
                                                family = stats::binomial(),
                                                SL.library = sl.library.infection, 
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
                                                X = sub_no_shig[,covariate_list_no_site , drop = FALSE],
                                                newX = data[,covariate_list_no_site , drop = FALSE],
                                                family = stats::binomial(),
                                                SL.library = sl.library.infection,
                                                cvControl = list(V = v_folds))
  
  tmp_pred_2a_2 <- prop_model_2a_2$SL.pred
  prop_vectors_2a$no_attr <- tmp_pred_2a_2 * (1 - prop_vectors_2a$inf_attr)
  
  ## Part 3: Propensity models for missingness
  
  # Data without site
  prop_covariates_inf_attr <- covariates_inf_attr[, colnames(covariates_inf_attr) %in% covariate_list_no_site , drop = FALSE]
  prop_covariates_no_attr <- covariates_no_attr[, colnames(covariates_no_attr) %in% covariate_list_no_site , drop = FALSE]
  
  ## Missingness model in infection attributable
  prop_model_3a <- SuperLearner::SuperLearner(Y = I_Y_inf_attr,
                                              X = data.frame(abx_inf_attr,
                                                             prop_covariates_inf_attr),
                                              family = stats::binomial(),
                                              SL.library = sl.library.missingness,
                                              cvControl = list(V = v_folds))
  
  ## Missingness model in no etiology 
  prop_model_3b <- SuperLearner::SuperLearner(Y = I_Y_no_attr,
                                              X = data.frame(abx_no_attr,
                                                             prop_covariates_no_attr),
                                              family = stats::binomial(),
                                              SL.library = sl.library.missingness,
                                              cvControl = list(V = v_folds))
  
  # Predict setting each abx level
  for(i in 1:length(abx_levels)){
    abx_level <- abx_levels[i]
    
    pred_data <- data
    pred_data[[abx_var_name]] <- abx_level
    
    prop_vectors_3a[,i] <- stats::predict(prop_model_3a, newdata = pred_data[,c(abx_var_name,
                                                                                covariate_list_no_site)], type = "response")$pred
    prop_vectors_3b[,i] <- stats::predict(prop_model_3b, newdata = pred_data[,c(abx_var_name,
                                                                                covariate_list_no_site)], type = "response")$pred
  }
  
  # -------------------------------------------------
  # STEP 3: AIPW estimates and confidence intervals
  # -------------------------------------------------
  
  ## Plug-in estimates
  
  plug_ins_inf <- colMeans(outcome_vectors_1a[inf_attr_idx, , drop = FALSE])
  plug_ins_no_attr <- colMeans(outcome_vectors_1b[inf_attr_idx, , drop = FALSE])
  
  ## Bias corrections
  inf_eifs <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  no_attr_eifs <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  
  colnames(inf_eifs) <- paste0("inf_eif_", abx_levels)
  colnames(no_attr_eifs) <- paste0("no_attr_eif_", abx_levels)
  
  for(i in 1:length(abx_levels)){
    abx_level <- abx_levels[i]
    
    # 1 - Bias correction for shigella (or other infection) attributable, abx level = a
    I_Inf_1 <- data[[infection_var_name]]
    P_Inf_1 <- mean(prop_vectors_2a[,1])
    
    I_Abx_a <- as.numeric(data[[abx_var_name]] == abx_level)
    P_Abx_a__Inf_1_Covariates <- prop_vectors_1a[,i]
    
    I_Delta_0 <- as.numeric(!is.na(data[[infection_var_name]])) # Indicator NOT missing
    P_Delta_0__Inf_all <- 1 - prop_vectors_3a[,i]
    
    obs_outcome <- ifelse(is.na(data[[laz_var_name]]), 0, data[[laz_var_name]])  
    Qbar_Inf_1_Abx_a_Covariates <- outcome_vectors_1a[,i]
    
    eif_vec_inf <- (I_Inf_1 / P_Inf_1) * (I_Abx_a / P_Abx_a__Inf_1_Covariates) * (I_Delta_0 / P_Delta_0__Inf_all) * (obs_outcome - Qbar_Inf_1_Abx_a_Covariates) +
      (I_Inf_1 / P_Inf_1) * (Qbar_Inf_1_Abx_a_Covariates - plug_ins_inf[i])
    
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
    
    I_Delta_0 <- as.numeric(!is.na(data[[infection_var_name]])) # Indicator NOT missing
    P_Delta_0__No_attr_all <- 1 - prop_vectors_3b[,i]
    
    obs_outcome <- ifelse(is.na(data[[laz_var_name]]), 0, data[[laz_var_name]])  
    Qbar_No_attr_Abx_a_Covariates <- outcome_vectors_1b[,i]
    
    Qbar_No_attr_Covariates <- outcome_vectors_1b[,i]
    
    I_Inf_1 <- data[[infection_var_name]]
    P_Inf_1 <- mean(prop_vectors_2a[,1])
    
    eif_vec_no_attr <- (I_No_attr_1 / P_No_attr__Covaritates) * (P_Inf_1__Covariates / P_Inf_1) * (I_Abx_a / P_Abx_a__Covariates) * (I_Delta_0 / P_Delta_0__No_attr_all) * (obs_outcome - Qbar_No_attr_Abx_a_Covariates) +
      (I_Inf_1 / P_Inf_1) * (Qbar_No_attr_Abx_a_Covariates - plug_ins_no_attr[i]) 
    
    inf_eifs[,i] <- eif_vec_inf
    no_attr_eifs[,i] <- eif_vec_no_attr
    
  }
  
  aipw_inf <- plug_ins_inf + colMeans(inf_eifs)
  aipw_no_attr <- plug_ins_no_attr + colMeans(no_attr_eifs)
  
  eif_matrix <- cbind(inf_eifs, no_attr_eifs)
  cov_matrix <- stats::cov(eif_matrix)
  
  # Compute effects for all levels of abx
  
  # Get id for each participant and recreate EIFs based on this if present
  if(!is.null(first_id_var_name)){
    first_id_eif_matrix <- cbind(data.frame(first_id = data[[first_id_var_name]]), eif_matrix)
    first_id_eif_matrix <- aggregate(. ~ first_id, data = first_id_eif_matrix, FUN = sum)
    
    scaled_matrix <- first_id_eif_matrix[,-c(1)] * (nrow(first_id_eif_matrix) / nrow(eif_matrix))
  }else{
    scaled_matrix <- eif_matrix
  }
  
  aipws_effect <- vector("numeric", length = length(abx_levels))
  eifs_effect <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(scaled_matrix)))
  
  names(aipws_effect) <- paste0("effect_", abx_levels)
  colnames(eifs_effect) <- paste0("effect_", abx_levels)
  
  for(i in 1:length(abx_levels)){
    
    # Get AIPW for Effect
    aipw_effect <- aipw_inf[i] - aipw_no_attr[i]
    aipws_effect[i] <- aipw_effect
    
    # Use EIFs to get variance for effect
    idx_1 <- i
    idx_2 <- length(abx_levels) + i
    
    gradient <- rep(0, length(abx_levels)*2)
    gradient[idx_1] <- 1
    gradient[idx_2] <- -1
    
    gradient <- matrix(gradient, ncol = 1)
    
    eif_effect <- as.numeric(as.matrix(scaled_matrix) %*% gradient)
    eifs_effect[,i] <- eif_effect
    
  }
  
  # Put all pt ests into dataframe same format as original gcomp analysis
  results_df <- data.frame(abx_levels = abx_levels,
                           abx_level_inf_1 = aipw_inf,
                           abx_level_inf_0 = aipw_no_attr,
                           effect_inf_abx_level = aipws_effect)
  
  # Add EIFs for effects to EIF matrix
  eif_matrix_scaled <- cbind(scaled_matrix, eifs_effect)
  cov_matrix <- stats::cov(eif_matrix_scaled)
  eif_hat <- sqrt( diag(cov_matrix) / nrow(eif_matrix_scaled) )
  
  # Get marginal effect estimates for MSM (if applicable) and create results object
  if(msm){
    marginal_effect_estimates <- colMeans(msm_vectors[inf_attr_idx, , drop = FALSE])
    
    results_object <- list(results_df = results_df,
                           eif_matrix = eif_matrix_scaled,
                           se = eif_hat,
                           marginal_effect_estimates = marginal_effect_estimates)
  } else{
    results_object <- list(results_df = results_df,
                           eif_matrix = eif_matrix_scaled,
                           se = eif_hat)
  }
  
  class(results_object) <- "aipw_other_diarrhea"

  if(return_models){
    # Make list of models
    if(msm){
      aipw_models <- list(outcome_model_1a = outcome_model_1a,
                          outcome_model_1b = outcome_model_1b,
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
  
}

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
#' 
#' @keywords internal
#' 
#' @returns List containing `aipw_other_diarrhea_2` object, models (if return_models = TRUE). `aipw_other_diarrhea_2` object contains dataframe with results, covariance matrix, standard errors.
aipw_other_diarrhea_2 <- function(data,
                                  laz_var_name,
                                  abx_var_name,
                                  infection_var_name,
                                  site_var_name,
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
                                  msm_formula = NULL){
  
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
  if(!is.null(no_etiology_var_name)){
    sub_no_attr <- data[which(data[[no_etiology_var_name]] == 1),]
  } else{
    I_no_attr <- ifelse(data[[infection_var_name]] == 0 & (rowSums(data[,pathogen_attributable_list]) == 0),
                        1, 0)
    data$no_etiology <- I_no_attr
    no_etiology_var_name <- "no_etiology"
    
    sub_no_attr <- data[which(data$no_etiology == 1),]
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
  
  # ------------------------------------------------------------
  # STEP 1: Fit & predict from outcome models
  # ------------------------------------------------------------
  
  ## Model 1a: Outcome model in Shigella (or other infection) attributable cases
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
  
  # For each level of abx, predict setting abx = x, infection = 1 & abx = x, infection = 0
  abx_levels <- unique(data[[abx_var_name]])[!is.na(unique(data[[abx_var_name]]))]
  
  outcome_vectors_1a <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  outcome_vectors_1b <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  outcome_vectors_2b <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  
  colnames(outcome_vectors_1a) <- paste0("abx_", abx_levels)
  colnames(outcome_vectors_1b) <- paste0("abx_", abx_levels)
  colnames(outcome_vectors_2b) <- paste0("abx_", abx_levels)
  
  if(return_models){
    outcome_model_2b_list <- vector("list", length = length(abx_levels))
  }
  
  # Iterate through each abx level
  for (i in 1:length(abx_levels)) {
    abx_level <- abx_levels[i]
    
    data_abx <- data
    data_abx[[abx_var_name]] <- abx_level
    
    ## 2a: Predictions for infection = 1 (no 2nd stage model needed)
    
    outcome_vectors_1a[, i] <- stats::predict(outcome_model_1a, newdata = data_abx[, c(abx_var_name,
                                                                                       covariate_list,
                                                                                       severity_list,
                                                                                       pathogen_quantity_list)], type = "response")$pred
    
    
    ## 2b: Predictions from 1b + Second stage regression model for no attribution
    
    outcome_vectors_1b[, i] <- stats::predict(outcome_model_1b, newdata = data_abx[, c(abx_var_name,
                                                                                       covariate_list,
                                                                                       severity_list,
                                                                                       pathogen_quantity_list)], type = "response")$pred
    
    # Get sub_no_attr_complete with yhat_level_0 predictions using obs_id
    data$set_abx_outcome <-  outcome_vectors_1b[, i]
    
    # Take new subset because added set_abx_outcome column
    sub_no_attr_2b <- data[which(data[[infection_var_name]] == 0), ]
    
    if (!is.na(no_etiology_var_name)) {
      sub_no_attr_2b <- sub_no_attr_2b[which(sub_no_attr_2b[[no_etiology_var_name]] == 1), ]
    } else {
      # otherwise find which rows have no attr pathogens in pathogen_attributable_list
      sub_no_attr_2b <- sub_no_attr_2b[which(rowSums(sub_no_attr_2b[, pathogen_attributable_list]) == 0), ]
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
      outcome_model_2b_list[[i]] <- outcome_model_2b
    }
  }
  
  if(msm){
    # Subtract predictions
    outcome_msm <- outcome_vectors_1a - outcome_vectors_2b
    
    msm_vectors <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
    colnames(msm_vectors) <- paste0("abx_", abx_levels)
    
    if(return_models){
      msm_model_list <- vector("list", length = length(abx_levels))
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
      
      effect_hetero_msm <- stats::glm(msm_formula_full, 
                                      data = msm_data[msm_data[[infection_var_name]] == 1,],       # QUESTION subset to shigella people here? then predict on everyone?
                                      family = outcome_type)                                       # QUESTION in general does this only work for continuous? because difference in outcome vectors?
      
      msm_vectors[,i] <- stats::predict(effect_hetero_msm, newdata = msm_data, type = 'response')
      
      if(return_models){
        msm_model_list[[i]] <- effect_hetero_msm
      }
      
    }
    
  }
  
  # ------------------------------------------------------------
  # STEP 2: Fit & predict from propensity models
  # ------------------------------------------------------------
  
  # If site in covariate_list, remove (any for if site one hot encoded)
  if(!any(is.na(site_var_name))){
    covariate_list_no_site <- covariate_list[!(covariate_list %in% site_var_name)]
  }
  
  # Create matrices to hold predictions from propensity models
  prop_vectors_1a <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  prop_vectors_1b <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  prop_vectors_2a <- data.frame(matrix(ncol = 2, nrow = nrow(data)))
  prop_vectors_3a <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  prop_vectors_3b <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  
  colnames(prop_vectors_1a) <- paste0("abx_", abx_levels)
  colnames(prop_vectors_1b) <- paste0("abx_", abx_levels)
  colnames(prop_vectors_2a) <- c("inf_attr", "no_attr")
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
      
      prop_covariates_inf_attr <- prop_sub_inf_attr[, covariate_list_no_site , drop = FALSE]
      prop_severity_inf_attr <- prop_sub_inf_attr[, severity_list, drop = FALSE]
      prop_pathogen_inf_attr <- prop_sub_inf_attr[, pathogen_quantity_list, drop = FALSE]
      
      prop_covariates_no_attr <- prop_sub_no_attr[, covariate_list_no_site , drop = FALSE]
      prop_severity_no_attr <- prop_sub_no_attr[, severity_list, drop = FALSE]
      prop_pathogen_no_attr <- prop_sub_no_attr[, pathogen_quantity_list, drop = FALSE]
      
      ## 1a. Propensity model for abx shigella attributable
      prop_model_1a <- SuperLearner::SuperLearner(Y = as.numeric(prop_sub_inf_attr[[abx_var_name]] == abx_level),
                                                  X = data.frame(prop_covariates_inf_attr,
                                                                 prop_severity_inf_attr,
                                                                 prop_pathogen_inf_attr),
                                                  newX = data[,c(covariate_list_no_site ,
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
                                                  newX = data[,c(covariate_list_no_site ,
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
                                                X = data[, covariate_list_no_site , drop = FALSE], 
                                                family = stats::binomial(),
                                                SL.library = sl.library.infection, 
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
                                                X = sub_no_shig[,covariate_list_no_site , drop = FALSE],
                                                newX = data[,covariate_list_no_site , drop = FALSE],
                                                family = stats::binomial(),
                                                SL.library = sl.library.infection,
                                                cvControl = list(V = v_folds))
  
  tmp_pred_2a_2 <- prop_model_2a_2$SL.pred
  prop_vectors_2a$no_attr <- tmp_pred_2a_2 * (1 - prop_vectors_2a$inf_attr)
  
  ###############################################
  ## Part 3: Propensity models for missingness ##
  ###############################################
  
  # Data without site
  prop_covariates_inf_attr <- covariates_inf_attr[, colnames(covariates_inf_attr) %in% covariate_list_no_site , drop = FALSE]
  prop_covariates_no_attr <- covariates_no_attr[, colnames(covariates_no_attr) %in% covariate_list_no_site , drop = FALSE]
  
  ## Missingness model in infection attributable
  prop_model_3a <- SuperLearner::SuperLearner(Y = I_Y_inf_attr,
                                              X = data.frame(abx_inf_attr,
                                                             prop_covariates_inf_attr,
                                                             severity_inf_attr,
                                                             pathogen_q_inf_attr),
                                              family = stats::binomial(),
                                              SL.library = sl.library.missingness,
                                              cvControl = list(V = v_folds))
  
  ## Missingness model in no etiology 
  prop_model_3b <- SuperLearner::SuperLearner(Y = I_Y_no_attr,
                                              X = data.frame(abx_no_attr,
                                                             prop_covariates_no_attr,
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
    
    prop_vectors_3a[,i] <- stats::predict(prop_model_3a, newdata = pred_data[,c(abx_var_name,
                                                                                covariate_list_no_site,
                                                                                severity_list,
                                                                                pathogen_quantity_list)], type = "response")$pred
    prop_vectors_3b[,i] <- stats::predict(prop_model_3b, newdata = pred_data[,c(abx_var_name,
                                                                                covariate_list_no_site,
                                                                                severity_list,
                                                                                pathogen_quantity_list)], type = "response")$pred
  }
  
  
  # -------------------------------------------------
  # STEP 3: AIPW estimates and  confidence intervals
  # -------------------------------------------------
  
  ## Plug-in estimates
  
  plug_ins_inf <- colMeans(outcome_vectors_1a[inf_attr_idx, , drop = FALSE])
  plug_ins_no_attr <- colMeans(outcome_vectors_2b[inf_attr_idx, , drop = FALSE])
  
  ## Bias corrections
  inf_eifs <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  no_attr_eifs <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  
  colnames(inf_eifs) <- paste0("inf_eif_", abx_levels)
  colnames(no_attr_eifs) <- paste0("no_attr_eif_", abx_levels)
  
  for(i in 1:length(abx_levels)){
    abx_level <- abx_levels[i]
    
    # 1 - Bias correction for shigella (or other infection) attributable, abx level = a
    I_Inf_1 <- data[[infection_var_name]]
    P_Inf_1 <- mean(prop_vectors_2a[,1])
    
    I_Abx_a <- as.numeric(data[[abx_var_name]] == abx_level)
    P_Abx_a__Inf_1_Covariates <- prop_vectors_1a[,i]
    
    I_Delta_0 <- as.numeric(!is.na(data[[infection_var_name]])) # Indicator NOT missing
    P_Delta_0__Inf_all <- 1 - prop_vectors_3a[,i]
    
    obs_outcome <- ifelse(is.na(data[[laz_var_name]]), 0, data[[laz_var_name]])  
    Qbar_Inf_1_Abx_a_Covariates <- outcome_vectors_1a[,i]
    
    eif_vec_inf <- (I_Inf_1 / P_Inf_1) * (I_Abx_a / P_Abx_a__Inf_1_Covariates) * (I_Delta_0 / P_Delta_0__Inf_all) * (obs_outcome - Qbar_Inf_1_Abx_a_Covariates) +
      (I_Inf_1 / P_Inf_1) * (Qbar_Inf_1_Abx_a_Covariates - plug_ins_inf[i])
    
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
    
    I_Delta_0 <- as.numeric(!is.na(data[[infection_var_name]])) # Indicator NOT missing
    P_Delta_0__No_attr_all <- 1 - prop_vectors_3b[,i]
    
    obs_outcome <- ifelse(is.na(data[[laz_var_name]]), 0, data[[laz_var_name]])  
    Qbar_No_attr_Abx_a_Covariates <- outcome_vectors_1b[,i]
    
    Qbar_No_attr_Covariates <- outcome_vectors_2b[,i]
    
    I_Inf_1 <- data[[infection_var_name]]
    P_Inf_1 <- mean(prop_vectors_2a[,1])
    
    eif_vec_no_attr <- (I_No_attr_1 / P_No_attr__Covaritates) * (P_Inf_1__Covariates / P_Inf_1) * (I_Abx_a / P_Abx_a__Covariates) * (I_Delta_0 / P_Delta_0__No_attr_all) * (obs_outcome - Qbar_No_attr_Abx_a_Covariates) +
      (I_No_attr_1 / P_No_attr__Covaritates) * (P_Inf_1__Covariates / P_Inf_1) * (Qbar_No_attr_Abx_a_Covariates - Qbar_No_attr_Covariates) + 
      (I_Inf_1 / P_Inf_1) * (Qbar_No_attr_Covariates - plug_ins_no_attr[i]) 
    
    inf_eifs[,i] <- eif_vec_inf
    no_attr_eifs[,i] <- eif_vec_no_attr
    
  }
  
  # Add bias correction
  aipw_inf <- plug_ins_inf + colMeans(inf_eifs)
  aipw_no_attr <- plug_ins_no_attr + colMeans(no_attr_eifs)
  
  eif_matrix <- cbind(inf_eifs, no_attr_eifs)
  cov_matrix <- stats::cov(eif_matrix)
  
  # Get id for each participant and recreate EIFs based on this if present
  if(!is.null(first_id_var_name)){
    first_id_eif_matrix <- cbind(data.frame(first_id = data[[first_id_var_name]]), eif_matrix)
    first_id_eif_matrix <- aggregate(. ~ first_id, data = first_id_eif_matrix, FUN = sum)
    
    scaled_matrix <- first_id_eif_matrix[,-c(1)] * (nrow(first_id_eif_matrix) / nrow(eif_matrix))
  } else{
    scaled_matrix <- eif_matrix
  }
  
  aipws_effect <- vector("numeric", length = length(abx_levels))
  eifs_effect <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(scaled_matrix)))
  
  names(aipws_effect) <- paste0("effect_", abx_levels)
  colnames(eifs_effect) <- paste0("effect_", abx_levels)
  
  # Compute effects for all levels of abx
  
  for(i in 1:length(abx_levels)){
    
    # Get AIPW for Effect
    aipw_effect <- aipw_inf[i] - aipw_no_attr[i]
    aipws_effect[i] <- aipw_effect
    
    # Use EIFs to get variance for effect
    idx_1 <- i
    idx_2 <- length(abx_levels) + i
    
    gradient <- rep(0, length(abx_levels)*2)
    gradient[idx_1] <- 1
    gradient[idx_2] <- -1
    
    gradient <- matrix(gradient, ncol = 1)
    
    eif_effect <- as.numeric(as.matrix(scaled_matrix) %*% gradient)
    eifs_effect[,i] <- eif_effect
    
  }
  
  # Put all pt ests into dataframe same format as original gcomp analysis
  results_df <- data.frame(abx_levels = abx_levels,
                           abx_level_inf_1 = aipw_inf,
                           abx_level_inf_0 = aipw_no_attr,
                           effect_inf_abx_level = aipws_effect)
  
  # Add EIFs for effects to EIF matrix
  eif_matrix_scaled <- cbind(scaled_matrix, eifs_effect)
  cov_matrix <- stats::cov(eif_matrix_scaled)
  eif_hat <- sqrt( diag(cov_matrix) / nrow(eif_matrix_scaled) )
  
  # Get marginal effect estimates for MSM (if applicable) and create results object
  if(msm){
    marginal_effect_estimates <- colMeans(msm_vectors[inf_attr_idx, , drop = FALSE])
    
    results_object <- list(results_df = results_df,
                           eif_matrix = eif_matrix_scaled,
                           se = eif_hat,
                           marginal_effect_estimates = marginal_effect_estimates)
  } else{
    results_object <- list(results_df = results_df,
                           eif_matrix = eif_matrix_scaled,
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
              se = eif_hat))
}


#' Function for AIPW estimates of effect of infection on growth- case-control analysis
#' 
#' @param data dataframe containing dataset used for gcomp
#' @param laz_var_name name of growth outcome variable
#' @param abx_var_name name of binary antibiotic variable
#' @param case_var_name name of variable indicating case (for case control analysis; case_control = TRUE)
#' @param site_var_name name of covariate for site (to exclude from propensity models)
#' @param covariate_list character vector containing names of baseline covariates
#' @param severity_list character vector containing names of severity-related covariates (post-infection). If NULL, use AIPW without second stage regression
#' @param pathogen_quantity_list character vector containing name(s) of pathogen quantity variables in dataset
#' @param outcome_type gaussian or binomial for continuous or binomial outcome
#' @param sl.library.outcome.case list containing SuperLearner libraries to use for case outcome regression 
#' @param sl.library.outcome.control list containing SuperLearner libraries to use for control outcome regression 
#' @param sl.library.treatment list containing SuperLearner libraries to use for antibiotic prescribing propensity model
#' @param sl.library.infection list containing SuperLearner libraries to use for infection propensity model
#' @param sl.library.missingness.case list containing SuperLearner libraries to use for outcome missingness model (cases)
#' @param sl.library.missingness.control list containing SuperLearner libraries to use for outcome missingness model (controls)
#' @param v_folds number of cross-validation folds to use in SuperLearner
#' @param return_models boolean return SuperLearner models. Default FALSE.
#' @param first_id_var_name name of variable indicating first_id in data. Used to account for re-enrollment in study. 
#' @param msm boolean indicating use of MSM for effect heterogeneity, default FALSE
#' @param msm_var_name name of variable to use for msm
#' @param msm_formula chatacter vector with formula to use for msm if msm TRUE
#' 
#' @keywords internal
#' 
#' @returns List containing `aipw_case_control` object, models (if return_models = TRUE). `aipw_case_control` object contains dataframe with results, covariance matrix, standard errors.
aipw_case_control <- function(data,
                              laz_var_name,
                              abx_var_name,
                              case_var_name,
                              site_var_name,
                              covariate_list,
                              severity_list,
                              pathogen_quantity_list = NULL,
                              outcome_type = "gaussian",
                              sl.library.outcome.case = c("SL.glm"),
                              sl.library.outcome.control = c("SL.glm"),
                              sl.library.treatment = c("SL.mean"),
                              sl.library.infection = c("SL.glm"),
                              sl.library.missingness.case = c("SL.glm"),
                              sl.library.missingness.control = c("SL.glm"),
                              v_folds = 5,
                              return_models = FALSE,
                              first_id_var_name = NULL,
                              msm = FALSE,
                              msm_var_name = NULL,
                              msm_formula = NULL){
  
  # ------------------------------------------------------------
  # STEP 0: Create subsets of data for model fitting
  # ------------------------------------------------------------
  
  # Rename abx levels if spaces in them
  abx_levels <- levels(factor(data[[abx_var_name]]))
  abx_levels_new <- ifelse(is.na(abx_levels), NA, gsub("[ /]", "_", abx_levels))
  data[[abx_var_name]] <- factor(data[[abx_var_name]], levels = abx_levels, labels = abx_levels_new)
  
  # get case vs control data
  case_data <- data[data[[case_var_name]] == 1,]
  control_data <- data[data[[case_var_name]] == 0,]
  
  case_data_idx <- which(data[[case_var_name]] == 1)
  control_data_idx <- which(data[[case_var_name]] == 0)
  
  # Case data prep
  I_Y_case <- ifelse(is.na(case_data[[laz_var_name]]), 1, 0)
  Y_case <- case_data[[laz_var_name]]
  covariates_case <- case_data[, covariate_list, drop = FALSE]
  severity_case <- case_data[, severity_list, drop = FALSE]
  pathogen_q_case <- case_data[, pathogen_quantity_list, drop = FALSE]
  abx_case <- case_data[, abx_var_name, drop = FALSE]
  
  # Complete case data (excluding missing outcome)
  case_data_complete <- case_data[!is.na(case_data[[laz_var_name]]), ]
  
  Y_case_complete <- case_data_complete[[laz_var_name]]
  covariates_case_complete <- case_data_complete[, covariate_list, drop = FALSE]
  severity_case_complete <- case_data_complete[, severity_list, drop = FALSE]
  pathogen_q_case_complete <- case_data_complete[, pathogen_quantity_list, drop = FALSE]
  abx_case_complete <- case_data_complete[, abx_var_name, drop = FALSE]
  
  # Control data prep
  I_Y_control <- ifelse(is.na(control_data[[laz_var_name]]), 1, 0)
  Y_control <- control_data[[laz_var_name]]
  if(is.null(covariate_list)){
    covariate_list <- covariate_list
  }
  covariates_control <- control_data[, covariate_list, drop = FALSE]
  
  # Complete control data (excluding missing outcome)
  control_data_complete <- control_data[!is.na(control_data[[laz_var_name]]), ]
  
  Y_control_complete <- control_data_complete[[laz_var_name]]
  covariates_control_complete <- control_data_complete[, covariate_list, drop = FALSE]
  
  # ------------------------------------------------------------
  # STEP 1: Fit & predict from outcome models
  # ------------------------------------------------------------
  
  ## Model 1a: Outcome model in cases
  outcome_model_1a <- SuperLearner::SuperLearner(Y = Y_case_complete, 
                                                 X = data.frame(abx_case_complete,
                                                                covariates_case_complete, 
                                                                severity_case_complete, 
                                                                pathogen_q_case_complete),
                                                 family = outcome_type, 
                                                 SL.library = sl.library.outcome.case,
                                                 cvControl = list(V = v_folds))
  
  ## Model 1b: Outcome model in controls
  outcome_model_1b <- SuperLearner::SuperLearner(Y = Y_control_complete, 
                                                 X = data.frame(covariates_control_complete),
                                                 family = outcome_type, 
                                                 SL.library = sl.library.outcome.control,
                                                 cvControl = list(V = v_folds))
  
  ## Predict outcomes
  
  # For each level of abx, predict setting abx = x, infection = 1 & abx = x, infection = 0
  abx_levels <- unique(case_data[[abx_var_name]])[!is.na(unique(case_data[[abx_var_name]]))]
  
  outcome_vectors_1a <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  outcome_vectors_1b <- data.frame(matrix(ncol = 1, nrow = nrow(data)))
  
  colnames(outcome_vectors_1a) <- paste0("abx_", abx_levels)
  colnames(outcome_vectors_1b) <- paste0("control")
  
  # Iterate through each abx level
  for (i in 1:length(abx_levels)) {
    abx_level <- abx_levels[i]
    
    data_abx <- case_data
    data_abx[[abx_var_name]] <- abx_level
    
    outcome_vectors_1a[case_data_idx, i] <- stats::predict(outcome_model_1a, newdata = data_abx[, c(abx_var_name,
                                                                                                    covariate_list,
                                                                                                    severity_list,
                                                                                                    pathogen_quantity_list)], type = "response")$pred
    
  }
  
  # Replace NA controls with 0
  outcome_vectors_1a[is.na(outcome_vectors_1a)] <- 0
  
  outcome_vectors_1b[, 1] <- stats::predict(outcome_model_1b, newdata = data[,covariate_list], type = "response")$pred
  
  if(msm){
    
    # Replicate outcome_vectors_1b across columns to match the dimensions of outcome_vectors_1a
    outcome_matrix_1b <- data.frame(matrix(outcome_vectors_1b[, 1], nrow = nrow(data), ncol = length(abx_levels), byrow = FALSE))
    
    # Subtract predictions
    outcome_msm <- outcome_vectors_1a - outcome_matrix_1b
    
    msm_vectors <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
    colnames(msm_vectors) <- paste0("abx_", abx_levels)
    
    if(return_models){
      msm_model_list <- vector("list", length = length(abx_levels))
    }
    
    for(i in 1:length(abx_levels)){
      abx_level <- abx_levels[i]
      abx_level_name <- paste0("abx_", abx_level)
      
      # Fit difference in predictions on msm_formula
      msm_formula_full <- as.formula(paste0(abx_level_name, " ~ ", msm_formula))
      
      # dataframe with difference in outcomes, variable of interest, infection var (for subsetting)
      msm_data <- setNames(
        data.frame(outcome_msm[, i], data[[msm_var_name]], data[[case_var_name]]), 
        c(abx_level_name, msm_var_name, case_var_name) 
      )
      
      effect_hetero_msm <- stats::glm(msm_formula_full, 
                                      data = msm_data[msm_data[[case_var_name]] == 1,],       # QUESTION subset to shigella people here? then predict on everyone?
                                      family = outcome_type)                                       # QUESTION in general does this only work for continuous? because difference in outcome vectors?
      
      msm_vectors[,i] <- stats::predict(effect_hetero_msm, newdata = msm_data, type = 'response')
      
      if(return_models){
        msm_model_list[[i]] <- effect_hetero_msm
      }
      
    }
    
  }
  
  # ------------------------------------------------------------
  # STEP 2: Fit & predict from propensity models
  # ------------------------------------------------------------
  
  # If site in covariate_list, remove
  if(!any(is.na(site_var_name))){
    covariate_list_no_site <- covariate_list[!(covariate_list %in% site_var_name)]
  }
  
  prop_vectors_1a <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  prop_vectors_2a <- data.frame(matrix(ncol = 1, nrow = nrow(data)))
  prop_vectors_3a <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  prop_vectors_3b <- data.frame(matrix(ncol = 1, nrow = nrow(data)))
  
  colnames(prop_vectors_1a) <- paste0("abx_", abx_levels)
  colnames(prop_vectors_2a) <- c("case")
  colnames(prop_vectors_3a) <- paste0("abx_", abx_levels)
  colnames(prop_vectors_3b) <- c("control")
  
  if(return_models){
    prop_model_1a_list <- vector("list", length = length(abx_levels))
  }
  
  ## Part 1: Propensity model for antibiotics 
  for(i in 1:length(abx_levels)){
    
    # Abx level to be modeled
    abx_level <- abx_levels[i]
    
    # If not last iteration
    if(i != length(abx_levels)){
      
      # Leave out previously modeled levels if applicable
      if(i > 1){
        abx_levels_out <- abx_levels[1:(i-1)]
        prop_sub_case <- case_data[-which(case_data[[abx_var_name]] %in% abx_levels_out),]
        
      } else{
        prop_sub_case <- case_data
      }
      
      prop_covariates_case <- prop_sub_case[, covariate_list_no_site, drop = FALSE]
      prop_severity_case <- prop_sub_case[, severity_list, drop = FALSE]
      prop_pathogen_case <- prop_sub_case[, pathogen_quantity_list, drop = FALSE]
      
      ## 1a. Propensity model for abx cases
      prop_model_1a <- SuperLearner::SuperLearner(Y = as.numeric(prop_sub_case[[abx_var_name]] == abx_level),
                                                  X = data.frame(prop_covariates_case,
                                                                 prop_severity_case,
                                                                 prop_pathogen_case),
                                                  family = stats::binomial(), 
                                                  SL.library = sl.library.treatment,
                                                  cvControl = list(V = v_folds))
      
      tmp_pred_a <- stats::predict(prop_model_1a, newdata = case_data[,c(covariate_list_no_site,
                                                                         severity_list,
                                                                         pathogen_quantity_list)], type = "response")$pred
      
      if(return_models){
        prop_model_1a_list[[i]] <- prop_model_1a
      }
      
      if(i == 1){
        # First prediction
        prop_vectors_1a[case_data_idx,i] <- tmp_pred_a
      } else{
        # Middle prediction
        # tmp_pred_a * for j in 1:(i-1) (1 - tmp_pred_aj) 
        for(j in 1:(i-1)){
          tmp_pred_a <- tmp_pred_a * (1 - prop_vectors_1a[case_data_idx,j])
        }
        prop_vectors_1a[case_data_idx,i] <- tmp_pred_a 
      }
      
    } else{
      # Last prediction
      prop_vectors_1a[case_data_idx,i] <- 1 - rowSums(prop_vectors_1a[case_data_idx,1:(ncol(prop_vectors_1a)-1)])
    }
    
  }
  
  # Fill in controls with 0
  prop_vectors_1a[is.na(prop_vectors_1a)] <- 0
  
  ## Part 2: Propensity models for shigella (or other infection) attribution
  
  # 2a_1 = Shigella Attributable ~ BL Cov
  prop_model_2a <- SuperLearner::SuperLearner(Y = data[[case_var_name]],
                                              X = data[, covariate_list_no_site, drop = FALSE],
                                              family = stats::binomial(),
                                              SL.library = sl.library.infection, 
                                              cvControl = list(V = v_folds))
  tmp_pred_2a <- prop_model_2a$SL.pred
  prop_vectors_2a[,1] <- tmp_pred_2a
  
  ## Part 3: Propensity models for missingness
  
  covariates_case_no_site <- case_data[, covariate_list_no_site, drop = FALSE]
  covariates_control_no_site <- control_data[, covariate_list_no_site, drop = FALSE]
  
  ## Missingness model in cases
  prop_model_3a <- SuperLearner::SuperLearner(Y = I_Y_case,
                                              X = data.frame(abx_case,
                                                             covariates_case_no_site,
                                                             severity_case,
                                                             pathogen_q_case),
                                              family = stats::binomial(),
                                              SL.library = sl.library.missingness.case,
                                              cvControl = list(V = v_folds))
  
  ## Missingness model in controls
  prop_model_3b <- SuperLearner::SuperLearner(Y = I_Y_control,
                                              X = data.frame(covariates_control_no_site),
                                              family = stats::binomial(),
                                              SL.library = sl.library.missingness.control,
                                              cvControl = list(V = v_folds))
  
  # Predict setting each abx level
  for(i in 1:length(abx_levels)){
    abx_level <- abx_levels[i]
    
    pred_data <- case_data_no_site
    pred_data[[abx_var_name]] <- abx_level
    
    prop_vectors_3a[case_data_idx,i] <- stats::predict(prop_model_3a, newdata = pred_data[,c(abx_var_name,
                                                                                             covariate_list_no_site,
                                                                                             pathogen_quantity_list,
                                                                                             severity_list)], type = "response")$pred
  }
  
  # Fill in NAs with 0
  prop_vectors_3a[is.na(prop_vectors_3a)] <- 0
  
  prop_vectors_3b[control_data_idx,1] <- prop_model_3b$SL.pred
  prop_vectors_3b[is.na(prop_vectors_3b)] <- 0
  
  # -------------------------------------------------
  # STEP 3: AIPW estimates and confidence intervals
  # -------------------------------------------------
  
  ## Plug-in estimates
  plug_ins_case <- colMeans(outcome_vectors_1a[case_data_idx,])
  plug_ins_control <- colMeans(outcome_vectors_1b[case_data_idx,])
  
  ## EIF for bias corrections
  case_eifs <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
  control_eifs <- data.frame(matrix(ncol = 1, nrow = nrow(data)))
  
  colnames(case_eifs) <- paste0("case_eif_", abx_levels)
  colnames(control_eifs) <- paste0("control_eif")
  
  for(i in 1:length(abx_levels)){
    abx_level <- abx_levels[i]
    
    # 1 - Bias correction for case, abx level = a
    
    I_Case <- data[[case_var_name]]
    P_Case <- mean(prop_vectors_2a[,1])
    
    I_Abx_a <- as.numeric(data[[abx_var_name]] == abx_level)
    # Replace NA with 0 (so doesn't error, will 0 out in bias term anyways)
    I_Abx_a[is.na(I_Abx_a)] <- 0
    
    P_Abx_a__Case_Covariates <- prop_vectors_1a[,i]
    
    I_Delta_0 <- as.numeric(!is.na(data[[laz_var_name]])) # Indicator NOT missing
    P_Delta_0__Case_all <- 1 - prop_vectors_3a[,i]
    
    obs_outcome <- ifelse(is.na(data[[laz_var_name]]), 0, data[[laz_var_name]])  
    Qbar_Case_Abx_a_Covariates <- outcome_vectors_1a[,i]
    
    eif_case_vec <- (I_Case / P_Case) * (I_Abx_a / P_Abx_a__Case_Covariates) * (I_Delta_0 / P_Delta_0__Case_all) * (obs_outcome - Qbar_Case_Abx_a_Covariates) +
      (I_Case / P_Case) * (Qbar_Case_Abx_a_Covariates - plug_ins_case[i])
    
    # correct for / 0 with P_Abx_a__Case_Covariates, should be 0'd out
    eif_case_vec <- ifelse(is.nan(eif_case_vec), 0, eif_case_vec)
    
    case_eifs[,i] <- eif_case_vec
    
  }
  
  # 2 - Bias correction for control
  
  I_control <- ifelse(data[[case_var_name]] == 1, 0, 1)
  P_control__Covariates <- 1 - prop_vectors_2a[,1]
  
  I_Delta_0 <- as.numeric(!is.na(data[[laz_var_name]])) # Indicator NOT missing
  I_Delta_0__Control_all <- 1 - prop_vectors_3b[,1]
  
  P_Case__all <- prop_vectors_2a[,1]
  
  Qbar_Control_Covariates <- outcome_vectors_1b[,1]
  
  eif_control_vec <- (I_control / P_control__Covariates) * (I_Delta_0 / I_Delta_0__Control_all) * (P_Case__all / P_Case) * (obs_outcome - Qbar_Control_Covariates) +
    (I_Case / P_Case) * (Qbar_Control_Covariates - plug_ins_control)
  
  control_eif <- eif_control_vec
  
  # Get AIPWs
  aipw_case <- plug_ins_case + colMeans(case_eifs)
  aipw_control <- plug_ins_control + colMeans(control_eif)  
  
  eif_matrix <- cbind(case_eifs, control_eif)
  
  # Compute effects for all levels of abx
  
  # Get id for each participant and recreate EIFs based on this if present
  if(!is.null(first_id_var_name)){
    first_id_eif_matrix <- cbind(data.frame(first_id = data[[first_id_var_name]]), eif_matrix)
    first_id_eif_matrix <- aggregate(. ~ first_id, data = first_id_eif_matrix, FUN = sum)
    scaled_matrix <- first_id_eif_matrix[,-c(1)] * (nrow(first_id_eif_matrix) / nrow(eif_matrix))
  }else{
    scaled_matrix <- eif_matrix
  }
  
  aipws_effect <- vector("numeric", length = length(abx_levels))
  eifs_effect <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(scaled_matrix)))
  
  names(aipws_effect) <- paste0("effect_", abx_levels)
  colnames(eifs_effect) <- paste0("effect_", abx_levels)
  
  for(i in 1:length(abx_levels)){
    aipw_effect <- aipw_case[i] - aipw_control[1]
    aipws_effect[i] <- aipw_effect
    
    idx_1 <- i
    idx_2 <- length(abx_levels) + 1
    
    gradient <- rep(0, length(abx_levels) + 1)
    
    gradient[idx_1] <- 1
    gradient[idx_2] <- -1
    
    gradient <- matrix(gradient, ncol = 1)
    eif_effect <- as.numeric(as.matrix(scaled_matrix) %*% gradient)
    eifs_effect[,i] <- eif_effect
    
  }
  
  results_df <- data.frame(abx_levels = abx_levels,
                           abx_level_case = aipw_case,
                           abx_level_control = rep(aipw_control, length(abx_levels)),
                           effect_inf_abx_level = aipws_effect)
  
  eif_matrix_scaled <- cbind(scaled_matrix, eifs_effect)
  cov_matrix <- stats::cov(eif_matrix_scaled)
  eif_hat <- sqrt( diag(cov_matrix) / nrow(eif_matrix_scaled) )
  
  # Get marginal effect estimates for MSM (if applicable) and create results object
  if(msm){
    marginal_effect_estimates <- colMeans(msm_vectors[case_data_idx, , drop = FALSE])
    
    results_object <- list(results_df = results_df,
                           eif_matrix = eif_matrix_scaled,
                           se = eif_hat,
                           marginal_effect_estimates = marginal_effect_estimates)
  } else{
    results_object <- list(results_df = results_df,
                           eif_matrix = eif_matrix_scaled,
                           se = eif_hat)
  }
  
  class(results_object) <- "aipw_case_control"
  
  if(return_models){
    # Make list of models
    if(msm){
      aipw_models <- list(outcome_model_1a = outcome_model_1a,
                          outcome_model_1b = outcome_model_1b,
                          prop_model_1a_list = prop_model_1a_list,
                          prop_model_2a = prop_model_2a,
                          prop_model_3a = prop_model_3a,
                          prop_model_3b = prop_model_3b,
                          msm_model_list = msm_model_list)
    } else{
      aipw_models <- list(outcome_model_1a = outcome_model_1a,
                          outcome_model_1b = outcome_model_1b,
                          prop_model_1a_list = prop_model_1a_list,
                          prop_model_2a = prop_model_2a,
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
}
