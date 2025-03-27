#' Function to do one bootstrap replicate
#' 
#' @param data dataframe containing dataset used for gcomp
#' @param laz_var_name name of growth outcome variable
#' @param abx_var_name name of binary antibiotic variable
#' @param infection_var_name name of binary infection variable
#' @param site_var_name name of site covariate in dataset (if applicable, else null)
#' @param site_interaction TRUE or FALSE indicating interaction between site and antibiotics 
#' @param age_var_name name of age covariate in dataset (if applicable to add spline, else NULL)
#' @param severity_list character vector containing names of severity-related covariates (post-infection). 
#' @param covariate_list character vector containing names of baseline covariates
#' @param case_control TRUE if case control analysis, FALSE otherwise
#' @param case_var_name name of variable indicating case (only needed if case_control = TRUE)
#' @param covariate_list_control covariate list for controls in case-control analysis, if not case control or covariate lists same then NULL
#' @param outcome_type gaussian or binomial for continuous or binomial outcome
#' @param att boolean if effect should be estimated among people who would naturally get infection, default TRUE
#' 
#' @keywords internal
#' 
#' @returns 
#' \describe{
#'  List of containing the following for the bootstrapped dataset:
#'  \item{\code{effect_inf_no_abx}}{numeric effect of infection on growth in subgroup that did *not* receive antibiotics}
#'  \item{\code{effect_inf_abx}}{numeric effect of infection on growth in subgroup that received antibiotics}
#'  \item{\code{abx_0_inf_1}}{expected growth outcome in infected subgroup who did not receive abx}
#'  \item{\code{abx_0_inf_0}}{expected growth outcome in uninfected subgroup who did not receive abx}
#'  \item{\code{abx_1_inf_1}}{expected growth outcome in infected subgroup who received abx}
#'  \item{\code{abx_1_inf_0}}{expected growth outcome in uninfected subgroup who receieved abx}
#'  }
one_boot <- function(data, 
                     laz_var_name,
                     abx_var_name,
                     infection_var_name,
                     site_var_name,
                     site_interaction,
                     age_var_name,
                     covariate_list, 
                     severity_list, 
                     case_var_name = NULL,
                     case_control = FALSE,
                     covariate_list_control = NULL, 
                     outcome_type = "gaussian",
                     att = TRUE){
  
  if(case_control == FALSE){
    ### Create Bootstrap Dataset ###
    
    # if dataset has 'first_id', there may be some children who are repeated in dataset- sample based on this
    # it's okay if boot_data slightly different size
    if("first_id" %in% colnames(data)){
      boot_child_id <- sample(unique(data$first_id), replace = TRUE)
      boot_child_rows_list <- sapply(boot_child_id, function(id){
        which(data$child_id == id)
      })
      boot_child_rows_vec <- Reduce(c, boot_child_rows_list)
      boot_data <- data[boot_child_rows_vec, , drop = FALSE]
    } else {
      boot_rows <- sample(1:nrow(data), replace = TRUE)
      boot_data <- data[boot_rows, , drop = FALSE]
    }
    
    ### Do g-computation using bootstrap data ###
    
    boot_res <- abx_growth_gcomp(boot_data, 
                                 laz_var_name = laz_var_name,
                                 abx_var_name = abx_var_name,
                                 infection_var_name = infection_var_name,
                                 site_var_name = site_var_name,
                                 site_interaction = site_interaction,
                                 covariate_list = covariate_list,
                                 severity_list = severity_list,
                                 age_var_name = age_var_name,
                                 att = att)
    return(boot_res)
    
  } else{
    ### Create Bootstrap Dataset ###
    
    # make sure to include the matched case/control for each sampled
    # make sure to include all instances of that child if sampled
    
    # GEMS:
    # first_id = associated with the child
    # case_id = associated with the case
    # child_id = associated with the episode 
    
    # make VIDA have same structure except first_id = child_id for now
    
    # (child_id = case_id if case, child_id = first_id if first time included in study)
    
    # TODO make more generalizable, go back and have these as inputs
    
    unique_first_id <- unique(data$first_id[data$case == 1])        # unique children who are cases in the study to sample from
    n_children <- length(unique_first_id)                           # number of unique cases in the study to sample
    i <- 0                                                          # counter for number of CASES to sample
    
    # ---------------------------------------------------------------
    # Helper function for child-level bootstrap sampling
    # ---------------------------------------------------------------
    get_matching <- function(sampled, data){
      # look for anywhere else with matching first_id
      match_first <- data[data$first_id %in% sampled,]
      
      # get case ids (episodes) each first id falls into
      episode_data <- data[data$case_id %in% unique(match_first$case_id),]
      
      # If new first_ids are found, repeat recursively
      if (!setequal(sampled, unique(episode_data$first_id))) {
        # Merge the new first_ids into sampled 
        sampled <- unique(c(sampled, episode_data$first_id))
        # Repeat function with new ids
        return(get_matching(sampled, data))
      } 
      
      # return rows of dataframe corresponding to the first ids needed
      return(data[data$first_id %in% sampled,])
      
    }
    
    
    # dataframe to hold bootstrap dataset
    boot_data <- data.frame()
    
    while(i < n_children){
      #Sample a child
      sampled <- sample(unique_first_id, 1)
      
      # Get rows corresponding to sampled child, 
      # all cases/controls matched to sample child, 
      # and any cases/controls matched to subsequent matched children
      sampled_rows <- get_matching(sampled, data)
      
      # Add to boot data
      boot_data <- rbind(boot_data, sampled_rows)
      
      # Add number of sampled children to children counter
      # i <- i + length(unique(sampled_rows$first_id))
      
      # Add number of sampled CASES to CASES counter
      i <- i + nrow(sampled_rows[sampled_rows[[case_var_name]] == 1,])
    }
    
    boot_res <- abx_growth_gcomp_case_control(data = boot_data,
                                              laz_var_name = laz_var_name,
                                              abx_var_name = abx_var_name,
                                              case_var_name = case_var_name,
                                              severity_list = severity_list,
                                              covariate_list = covariate_list,
                                              covariate_list_control = covariate_list_control,
                                              site_var_name = site_var_name,
                                              site_interaction = site_interaction, 
                                              age_var_name = age_var_name, 
                                              outcome_type = outcome_type,
                                              att = att)
    return(boot_res)
    
  }
}

#' Function to do n_boot bootstrap replicates and get bootstrap standard error
#' 
#' @param data dataframe containing dataset used for gcomp
#' @param laz_var_name name of growth outcome variable
#' @param abx_var_name name of binary antibiotic variable
#' @param infection_var_name name of binary infection variable
#' @param site_var_name name of site covariate in dataset (if applicable, else null)
#' @param site_interaction TRUE or FALSE indicating interaction between site and antibiotics 
#' @param age_var_name name of age covariate in dataset (if applicable to add spline, else NULL)
#' @param severity_list character vector containing names of severity-related covariates (post-infection). If NULL, perform traditional gcomputation. Else, perform longitudinal gcomputation
#' @param covariate_list character vector containing names of baseline covariates
#' @param n_boot number of bootstrap replicates to repeat
#' @param case_control TRUE if case control analysis, FALSE otherwise
#' @param case_var_name name of variable indicating case (only needed if case_control = TRUE)
#' @param covariate_list_control covariate list for controls in case-control analysis, if not case control or covariate lists same then NULL
#' @param outcome_type gaussian or binomial for continuous or binomial outcome
#' @param att boolean if effect should be estimated among people who would naturally get infection, default TRUE
#' 
#' @keywords internal
#' 
#' @returns 
#' \describe{
#'  List of containing bootstrap standard error and confidence intervals for each effect of interest:
#'  \item{\code{se_abx}}{boostrap standard error for subgroup that received antibiotics}
#'  \item{\code{lower_se_abx}}{lower bound of 95% confidence interval for the subgroup that received antibiotics}
#'  \item{\code{upper_se_abx}}{upper bound of 95% confidence interval for the subgroup that received antibiotics}
#'  \item{\code{se_no_abx}}{boostrap standard error for subgroup that did *not* receive antibiotics}
#'  \item{\code{lower_se_no_abx}}{lower bound of 95% confidence interval for the subgroup that did *not* receive antibiotics}
#'  \item{\code{upper_se_no_abx}}{upper bound of 95% confidence interval for the subgroup that that did *not* receive antibiotics}
#'  \item{\code{se_abx_0_inf_1}}{boostrap standard error for infected subgroup who did *not* receive abx and was infected}
#'  \item{\code{lower_ci_abx_0_inf_1}}{lower bound of 95% confidence interval for the subgroup that did *not* receive antibiotics and was infected}
#'  \item{\code{upper_ci_abx_0_inf_1}}{upper bound of 95% confidence interval for the subgroup that did *not* receive antibiotics and was infected}
#'  \item{\code{se_abx_0_inf_0}}{boostrap standard error for uninfected subgroup who did not receive abx}
#'  \item{\code{lower_ci_abx_0_inf_0}}{lower bound of 95% confidence interval for the subgroup that did *not* receive antibiotics and was not infected}
#'  \item{\code{upper_ci_abx_0_inf_0}}{upper bound of 95% confidence interval for the subgroup that did *not* receive antibiotics and was not infected}
#'  \item{\code{se_abx_1_inf_1}}{boostrap standard error for infected subgroup who received abx}
#'  \item{\code{lower_ci_abx_1_inf_1}}{lower bound of 95% confidence interval for the subgroup that did receive antibiotics and was infected}
#'  \item{\code{upper_ci_abx_1_inf_1}}{upper bound of 95% confidence interval for the subgroup that did receive antibiotics and was infected}
#'  \item{\code{se_abx_1_inf_0}}{boostrap standard error for uninfected subgroup who receieved abx - only in comparison to other diarrhea}
#'  \item{\code{lower_ci_abx_1_inf_0}}{lower bound of 95% confidence interval for the subgroup that did receive antibiotics and was not infected - only in comparison to other diarrhea}
#'  \item{\code{upper_ci_abx_1_inf_0}}{upper bound of 95% confidence interval for the subgroup that did receive antibiotics and was not infected - only in comparison to other diarrhea}
#'  }
bootstrap_estimates <- function(data, 
                                laz_var_name,
                                abx_var_name,
                                site_var_name,
                                site_interaction,
                                age_var_name,
                                covariate_list, 
                                severity_list,
                                n_boot, 
                                infection_var_name = NULL,
                                case_var_name = NULL,
                                case_control = FALSE,
                                covariate_list_control = NULL, 
                                outcome_type = "gaussian",
                                att = TRUE){
  
  # Replicate one_boot function n_boot times
  boot_estimates <- replicate(n_boot, one_boot(data = data, 
                                               laz_var_name = laz_var_name,
                                               abx_var_name = abx_var_name,
                                               infection_var_name = infection_var_name,
                                               site_var_name = site_var_name,
                                               site_interaction = site_interaction,
                                               covariate_list = covariate_list,
                                               severity_list = severity_list,
                                               age_var_name = age_var_name,
                                               case_var_name = case_var_name,
                                               case_control = case_control,
                                               covariate_list_control = covariate_list_control, 
                                               outcome_type = outcome_type,
                                               att = att), simplify = FALSE) 
  
  # Extract the coefs
  # coef_vecs <- lapply(boot_estimates, function(x) x[['outcome_coefs']])
  # ses_vecs <- lapply(boot_estimates, function(x) x[['ses_coefs']])
  # 
  # # Combine the numeric vectors into a single dataframe
  # coefs_df <- data.frame(do.call(rbind, coef_vecs))
  # ses_df <- data.frame(do.call(rbind, ses_vecs))
  # 
  boot_res <- data.frame(do.call(rbind, boot_estimates))
  
  abx_levels <- unique(boot_res$abx_levels)
  
  # Make list for output with entry for each abx level
  out <- vector("list", length = length(abx_levels))
  names(out) <- paste0("Abx = ", abx_levels)
  
  for(i in 1:length(abx_levels)){
    abx_level <- abx_levels[i]
    
    sub_boot_res <- boot_res[boot_res$abx_levels == abx_level,]
    
    out[[i]]$abx_level <- abx_level
    
    if(case_control == FALSE){
      out[[i]]$se_abx_level_inf_1 <- stats::sd(sub_boot_res$abx_level_inf_1)
      out[[i]]$lower_ci_abx_level_inf_1 <- stats::quantile(sub_boot_res$abx_level_inf_1, p = 0.025, names = FALSE)
      out[[i]]$upper_ci_abx_level_inf_1 <- stats::quantile(sub_boot_res$abx_level_inf_1, p = 0.975, names = FALSE)
      
      out[[i]]$se_abx_level_inf_0 <- stats::sd(sub_boot_res$abx_level_inf_0)
      out[[i]]$lower_ci_abx_level_inf_0 <- stats::quantile(sub_boot_res$abx_level_inf_0, p = 0.025, names = FALSE)
      out[[i]]$upper_ci_abx_level_inf_0 <- stats::quantile(sub_boot_res$abx_level_inf_0, p = 0.975, names = FALSE)
      
    } else{
      
      out[[i]]$se_abx_level_case <- stats::sd(sub_boot_res$abx_level_case)
      out[[i]]$lower_ci_abx_level_case <- stats::quantile(sub_boot_res$abx_level_case, p = 0.025, names = FALSE)
      out[[i]]$upper_ci_abx_level_case <- stats::quantile(sub_boot_res$abx_level_case, p = 0.975, names = FALSE)
      
      out[[i]]$se_abx_level_control <- stats::sd(sub_boot_res$abx_level_control)
      out[[i]]$lower_ci_abx_level_control <- stats::quantile(sub_boot_res$abx_level_control, p = 0.025, names = FALSE)
      out[[i]]$upper_ci_abx_level_control <- stats::quantile(sub_boot_res$abx_level_control, p = 0.975, names = FALSE)

    }
    
    out[[i]]$se_effect_inf_abx_level <- stats::sd(sub_boot_res$effect_inf_abx_level)
    out[[i]]$lower_ci_effect_inf_abx_level<- stats::quantile(sub_boot_res$effect_inf_abx_level, p = 0.025, names = FALSE)
    out[[i]]$upper_ci_effect_inf_abx_level <- stats::quantile(sub_boot_res$effect_inf_abx_level, p = 0.975, names = FALSE)
    
  }
  
  return(out)
}
