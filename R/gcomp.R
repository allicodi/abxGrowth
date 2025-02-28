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
                             sl.library.outcome.2 = NULL, 
                             sl.library.treatment = NULL,
                             sl.library.infection = NULL,
                             sl.library.missingness = NULL,
                             v_folds = 5,
                             att = TRUE,
                             return_models = FALSE){
  
  if(!is.null(severity_list)){
    
    # -----------------------------------
    # Longitudinal G-Computation 
    # -----------------------------------
    
    ### STEP 0: Create subsets of data for model fitting
    
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
    if(!is.na(no_etiology_var_name)){
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
    
    ### STEP 1: Fit outcome models
    
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
    
    ### STEP 2: Fit second stage outcome regressions & predict outcomes
    
    # For each level of abx, predict setting abx = x, infection = 1 & abx = x, infection = 0
    abx_levels <- unique(data[[abx_var_name]])[!is.na(unique(data[[abx_var_name]]))]
    
    outcome_vectors_1a <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
    outcome_vectors_1b <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
    outcome_vectors_2b <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
    
    colnames(outcome_vectors_1a) <- paste0("abx_", abx_levels)
    colnames(outcome_vectors_1b) <- paste0("abx_", abx_levels)
    colnames(outcome_vectors_2b) <- paste0("abx_", abx_levels)
    
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
      
      outcome_model_2b <- SuperLearner::SuperLearner(
        Y = sub_no_attr_2b[['set_abx_outcome']],
        X = sub_no_attr_2b[, covariate_list, drop = FALSE],
        family = outcome_type,
        SL.library = sl.library.outcome.2, 
        cvControl = list(V = v_folds)
      )
      
      outcome_vectors_2b[, i] <- stats::predict(outcome_model_2b, newdata = data[, c(covariate_list)], type = "response")$pred
      
    }
    
    
    ### STEP 3: Propensity models
    
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
                                                 X = data[, covariate_list, drop = FALSE], 
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
                                                  X = sub_no_shig[,covariate_list, drop = FALSE],
                                                  newX = data[,covariate_list, drop = FALSE],
                                                  family = stats::binomial(),
                                                  SL.library = sl.library.infection,
                                                  cvControl = list(V = v_folds))
    
    tmp_pred_2a_2 <- prop_model_2a_2$SL.pred
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
      
      prop_vectors_3a[,i] <- stats::predict(prop_model_3a, newdata = pred_data[,c(abx_var_name,
                                                                                  covariate_list,
                                                                                  severity_list,
                                                                                  pathogen_quantity_list)], type = "response")$pred
      prop_vectors_3b[,i] <- stats::predict(prop_model_3b, newdata = pred_data[,c(abx_var_name,
                                                                                  covariate_list,
                                                                                  severity_list,
                                                                                  pathogen_quantity_list)], type = "response")$pred
    }
   
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
      
      bias_correct_inf <- (I_Inf_1 / P_Inf_1) * (I_Abx_a / P_Abx_a__Inf_1_Covariates) * (I_Delta_0 / P_Delta_0__Inf_all) * (obs_outcome - Qbar_Inf_1_Abx_a_Covariates) +
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
      
      bias_correction_no_attr <- (I_No_attr_1 / P_No_attr__Covaritates) * (P_Inf_1__Covariates / P_Inf_1) * (I_Abx_a / P_Abx_a__Covariates) * (I_Delta_0 / P_Delta_0__No_attr_all) * (obs_outcome - Qbar_No_attr_Abx_a_Covariates) +
        (I_No_attr_1 / P_No_attr__Covaritates) * (P_Inf_1__Covariates / P_Inf_1) * (Qbar_No_attr_Abx_a_Covariates - Qbar_No_attr_Covariates) +
        (I_Inf_1 / P_Inf_1) * (Qbar_No_attr_Covariates - plug_ins_no_attr[i])
    
      inf_eifs[,i] <- bias_correct_inf
      no_attr_eifs[,i] <- bias_correction_no_attr
      
    }
    
    aipw_inf <- plug_ins_inf + colMeans(inf_eifs)
    aipw_no_attr <- plug_ins_no_attr + colMeans(no_attr_eifs)
    
    eif_matrix <- cbind(inf_eifs, no_attr_eifs)
    cov_matrix <- stats::cov(eif_matrix)
    
    # Compute effects for all levels of abx
    
    aipws_effect <- vector("numeric", length = length(abx_levels))
    vars_effect <- vector("numeric", length = length(abx_levels))
    eifs_effect <- data.frame(matrix(ncol = length(abx_levels), nrow = nrow(data)))
    
    names(aipws_effect) <- paste0("effect_", abx_levels)
    names(vars_effect) <- paste0("effect_", abx_levels)
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
      
      # QUESTION i don't think we need variance if using EIF at end?? 
      # was just trying to copy from my thesis
      var_effect <- t(gradient) %*% cov_matrix %*% gradient
      eif_effect <- as.numeric(as.matrix(eif_matrix) %*% gradient)
      eifs_effect[,i] <- eif_effect
      vars_effect[i] <- var_effect
      
    }
    
    # Put all pt ests into dataframe same format as original gcomp analysis
    results_df <- data.frame(abx_levels = abx_levels,
                             abx_level_inf_1 = aipw_inf,
                             abx_level_inf_0 = aipw_no_attr,
                             effect_inf_abx_level = aipws_effect)
    
    # Add EIFs for effects to EIF matrix
    eif_matrix <- cbind(eif_matrix, eifs_effect)
    cov_matrix <- stats::cov(eif_matrix)
    eif_hat <- sqrt( diag(cov_matrix) / nrow(data) )
    
    return(list(results_df = results_df,
                eif_matrix = eif_matrix,
                se = eif_hat))
    
    # OLD FORMAT RESULTS
    # (try to put in this form at end so don't have to redo printing fn again?? )
    
    # results_df <- data.frame(abx_levels = abx_levels,
    #                          abx_level_inf_1 = vector("numeric", length = length(abx_levels)),
    #                          abx_level_inf_0 = vector("numeric", length = length(abx_levels)),
    #                          effect_inf_abx_level = vector("numeric", length = length(abx_levels)))
    
    
    # Return list with results dataframe (old version, effect est only)
    #return(results_df)
    
  } else {
    
    stop("Skip for now")
    
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

