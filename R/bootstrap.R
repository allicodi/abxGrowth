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
                     outcome_type = "gaussian"){
  
  if(case_control == FALSE){
    ### Create Bootstrap Dataset ###
    
    # if dataset has 'child_id', there may be some children who are repeated in dataset- sample based on this
    # it's okay if boot_data slightly different size
    if("child_id" %in% colnames(data)){
      boot_child_id <- sample(unique(data$child_id), replace = TRUE)
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
                                 age_var_name = age_var_name)
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
    
    # NOT POSITIVE THIS IS WORKING RIGHT
   
    # get unique pids for episodes
    boot_child_id <- sample(unique(data$child_id), replace = TRUE)

    i <- 1
    n <- 0
    final_idxs <- c()

    # ----------------------------------------------------------------------------
    # Recursive helper function to get all indexes that are associated with
    # cases + controls for the given episode's case and control idxs
    # ----------------------------------------------------------------------------
    get_all_matching <- function(idxs){
      # Find any other instances of the child in the dataset and the matches that go along with that
      first_ids <- data$first_id[idxs]
      first_id_idxs <- which(data$first_id %in% first_ids)

      # Add any new observations to idxs
      idxs <- unique(c(idxs, first_id_idxs))

      # Make sure all cases/controls associated with those indexes are included
      case_ids <- data$case_id[idxs]
      child_ids <- data$child_id[idxs]
      new_idxs <- unique(c(which(data$case_id %in% case_ids), which(data$child_id %in% child_ids)))

      if (!identical(unique(new_idxs), unique(idxs))) {
        #print("this is where recursion would be helpful but brain blah")
        get_all_matching(new_idxs)
      } else{
        return(new_idxs)
      }
    }

    while(n < nrow(data)){
      child_idx <- which(data$child_id == boot_child_id[i])

      if(data[child_idx, "case"] == 1){
        # Sampled Case; find the control(s) that match
        idxs <- which(data$case_id == data$case_id[child_idx])
        final_idxs <- c(final_idxs, get_all_matching(idxs))

      } else{
        # Sampled Control; find the case and any other controls that match
        idxs <- which(data$case_id == data[child_idx, "case_id"])
        final_idxs <- c(final_idxs, get_all_matching(idxs))

      }

      n <- length(final_idxs)
      i <- i + 1
    }
    
    boot_data <- data[final_idxs, , drop = FALSE]
    
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
                                              outcome_type = outcome_type)
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
                                outcome_type = "gaussian"){
  
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
                                               outcome_type = outcome_type), simplify = FALSE) 
  
  boot_res <- data.frame(do.call(rbind, boot_estimates))
  boot_res$effect_inf_no_abx <- unlist(boot_res$effect_inf_no_abx)
  boot_res$effect_inf_abx <- unlist(boot_res$effect_inf_abx)
  boot_res$abx_0_inf_1 <- unlist(boot_res$abx_0_inf_1)
  boot_res$abx_0_inf_0 <- unlist(boot_res$abx_0_inf_0)
  boot_res$abx_1_inf_1 <- unlist(boot_res$abx_1_inf_1)
  
  if(case_control == FALSE){
    boot_res$abx_1_inf_0 <- unlist(boot_res$abx_1_inf_0)
  }
  
  out <- list()
  
  # Get standard error and confidence intervals for each effect 
  
  out$se_abx <- stats::sd(boot_res$effect_inf_abx)
  out$lower_ci_abx <- stats::quantile(boot_res$effect_inf_abx, p = 0.025, names = FALSE)
  out$upper_ci_abx <- stats::quantile(boot_res$effect_inf_abx, p = 0.975, names = FALSE)
  
  out$se_no_abx <- stats::sd(boot_res$effect_inf_no_abx)
  out$lower_ci_no_abx <- stats::quantile(boot_res$effect_inf_no_abx, p = 0.025, names = FALSE)
  out$upper_ci_no_abx <- stats::quantile(boot_res$effect_inf_no_abx, p = 0.975, names = FALSE)
  
  if(case_control == FALSE){
    out$se_abx_1_inf_0 <- stats::sd(boot_res$abx_1_inf_0)
    out$lower_ci_abx_1_inf_0 <- stats::quantile(boot_res$abx_1_inf_0, p = 0.025, names = FALSE)
    out$upper_ci_abx_1_inf_0 <- stats::quantile(boot_res$abx_1_inf_0, p = 0.975, names = FALSE)
  }
  
  out$se_abx_0_inf_1 <- stats::sd(boot_res$abx_0_inf_1)
  out$lower_ci_abx_0_inf_1 <- stats::quantile(boot_res$abx_0_inf_1, p = 0.025, names = FALSE)
  out$upper_ci_abx_0_inf_1 <- stats::quantile(boot_res$abx_0_inf_1, p = 0.975, names = FALSE)
  
  out$se_abx_0_inf_0 <- stats::sd(boot_res$abx_0_inf_0)
  out$lower_ci_abx_0_inf_0 <- stats::quantile(boot_res$abx_0_inf_0, p = 0.025, names = FALSE)
  out$upper_ci_abx_0_inf_0 <- stats::quantile(boot_res$abx_0_inf_0, p = 0.975, names = FALSE)
  
  out$se_abx_1_inf_1 <- stats::sd(boot_res$abx_1_inf_1)
  out$lower_ci_abx_1_inf_1 <- stats::quantile(boot_res$abx_1_inf_1, p = 0.025, names = FALSE)
  out$upper_ci_abx_1_inf_1 <- stats::quantile(boot_res$abx_1_inf_1, p = 0.975, names = FALSE)
  
  return(out)
}
