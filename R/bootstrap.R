#' Function to do one bootstrap replicate
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
                     severity_list){
  
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
#'  \item{\code{se_abx_1_inf_0}}{boostrap standard error for uninfected subgroup who receieved abx}
#'  \item{\code{lower_ci_abx_1_inf_0}}{lower bound of 95% confidence interval for the subgroup that did receive antibiotics and was not infected}
#'  \item{\code{upper_ci_abx_1_inf_0}}{upper bound of 95% confidence interval for the subgroup that did receive antibiotics and was not infected}
#'  }
bootstrap_estimates <- function(data, 
                                laz_var_name,
                                abx_var_name,
                                infection_var_name,
                                site_var_name,
                                site_interaction,
                                age_var_name,
                                covariate_list, 
                                severity_list,
                                n_boot){
  
  # Replicate one_boot function n_boot times
  boot_estimates <- replicate(n_boot, one_boot(data = data, 
                                               laz_var_name = laz_var_name,
                                               abx_var_name = abx_var_name,
                                               infection_var_name = infection_var_name,
                                               site_var_name = site_var_name,
                                               site_interaction = site_interaction,
                                               covariate_list = covariate_list,
                                               severity_list = severity_list,
                                               age_var_name = age_var_name), simplify = FALSE) 
  
  boot_res <- data.frame(do.call(rbind, boot_estimates))
  boot_res$effect_inf_no_abx <- unlist(boot_res$effect_inf_no_abx)
  boot_res$effect_inf_abx <- unlist(boot_res$effect_inf_abx)
  boot_res$abx_0_inf_1 <- unlist(boot_res$abx_0_inf_1)
  boot_res$abx_0_inf_0 <- unlist(boot_res$abx_0_inf_0)
  boot_res$abx_1_inf_1 <- unlist(boot_res$abx_1_inf_1)
  boot_res$abx_1_inf_0 <- unlist(boot_res$abx_1_inf_0)
  
  out <- list()
  
  # Get standard error and confidence intervals for each effect 
  
  out$se_abx <- stats::sd(boot_res$effect_inf_abx)
  out$lower_ci_abx <- stats::quantile(boot_res$effect_inf_abx, p = 0.025, names = FALSE)
  out$upper_ci_abx <- stats::quantile(boot_res$effect_inf_abx, p = 0.975, names = FALSE)
  
  out$se_no_abx <- stats::sd(boot_res$effect_inf_no_abx)
  out$lower_ci_no_abx <- stats::quantile(boot_res$effect_inf_no_abx, p = 0.025, names = FALSE)
  out$upper_ci_no_abx <- stats::quantile(boot_res$effect_inf_no_abx, p = 0.975, names = FALSE)
  
  out$se_abx_0_inf_1 <- stats::sd(boot_res$abx_0_inf_1)
  out$lower_ci_abx_0_inf_1 <- stats::quantile(boot_res$abx_0_inf_1, p = 0.025, names = FALSE)
  out$upper_ci_abx_0_inf_1 <- stats::quantile(boot_res$abx_0_inf_1, p = 0.975, names = FALSE)
  
  out$se_abx_0_inf_0 <- stats::sd(boot_res$abx_0_inf_0)
  out$lower_ci_abx_0_inf_0 <- stats::quantile(boot_res$abx_0_inf_0, p = 0.025, names = FALSE)
  out$upper_ci_abx_0_inf_0 <- stats::quantile(boot_res$abx_0_inf_0, p = 0.975, names = FALSE)
  
  out$se_abx_1_inf_1 <- stats::sd(boot_res$abx_1_inf_1)
  out$lower_ci_abx_1_inf_1 <- stats::quantile(boot_res$abx_1_inf_1, p = 0.025, names = FALSE)
  out$upper_ci_abx_1_inf_1 <- stats::quantile(boot_res$abx_1_inf_1, p = 0.975, names = FALSE)
  
  out$se_abx_1_inf_0 <- stats::sd(boot_res$abx_1_inf_0)
  out$lower_ci_abx_1_inf_0 <- stats::quantile(boot_res$abx_1_inf_0, p = 0.025, names = FALSE)
  out$upper_ci_abx_1_inf_0 <- stats::quantile(boot_res$abx_1_inf_0, p = 0.975, names = FALSE)
  
  return(out)
}

#' Main function to get point estimates and bootstrap confidence intervals for given data and parameters
#' 
#' @param data dataframe containing dataset used for gcomp
#' @param laz_var_name name of growth outcome variable
#' @param abx_var_name name of binary antibiotic variable
#' @param infection_var_name name of binary infection variable
#' @param site_var_name name of site variable in dataset
#' @param site_interaction TRUE or FALSE indicating interaction between site and antibiotics 
#' @param age_var_name name of age covariate in dataset (if applicable to add spline, else NULL)
#' @param severity_list character vector containing names of severity-related covariates (post-infection). If NULL, perform traditional gcomputation. Else, perform longitudinal gcomputation
#' @param covariate_list character vector containing names of baseline covariates
#' @param n_boot number of bootstrap replicates to repeat
#' @param seed seed to set for bootstrap 
#' @param case_control TRUE if case control analysis, FALSE otherwise
#' @param case_var_name name of variable indicating case (only needed if case_control = TRUE)
#' 
#' @export
#' 
#' @returns List of class `aggcomp_res` containing point estimates, standard error, and 95% confidence intervals for each effect of interest
aggcomp <- function(data,
                    laz_var_name,
                    abx_var_name,
                    infection_var_name = "hazdiff",
                    site_var_name,
                    site_interaction,
                    age_var_name,
                    covariate_list, 
                    severity_list,
                    n_boot = 1000,
                    seed = 12345,
                    case_control = FALSE,
                    case_var_name = "case"){
                      
  # set seed
  set.seed(seed)
  
  # Get point estimates for effect_inf_no_abx, effect_inf_abx, and the estimated growth that makes up each
  pt_est <- abx_growth_gcomp(data = data, 
                             laz_var_name = laz_var_name,
                             abx_var_name = abx_var_name,
                             infection_var_name = infection_var_name,
                             site_var_name = site_var_name,
                             site_interaction = site_interaction,
                             age_var_name = age_var_name, 
                             covariate_list = covariate_list,
                             severity_list = severity_list,
                             case_control = case_control,
                             case_var_name = case_var_name)
  
  # Get standard error and confidence intervals for those point estimates
  bootstrap_results <- bootstrap_estimates(data = data, 
                                           laz_var_name = laz_var_name,
                                           abx_var_name = abx_var_name,
                                           infection_var_name = infection_var_name,
                                           site_var_name = site_var_name,
                                           site_interaction = site_interaction,
                                           age_var_name = age_var_name, 
                                           covariate_list = covariate_list,
                                           severity_list = severity_list,
                                           n_boot = n_boot)
  
  # Return list of all point estimates and standard errors
  # Class `aggcomp_res`
  results <- c(unlist(pt_est), 
               unlist(bootstrap_results))
  
  parameters <- list(laz_var_name = laz_var_name,
                     abx_var_name = abx_var_name,
                     infection_var_name = infection_var_name,
                     site_var_name = site_var_name,
                     site_interaction = site_interaction,
                     age_var_name = age_var_name, 
                     covariate_list = covariate_list,
                     severity_list = severity_list)
  
  aggcomp_results <- list(results = results,
                          parameters = parameters)
  
  class(aggcomp_results) <- "aggcomp_res"
  
  return(aggcomp_results)
  
}



