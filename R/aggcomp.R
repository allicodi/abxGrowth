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
#' @param covariate_list_control covariate list for controls in case-control analysis, if not case control or covariate lists same then NULL
#' @param outcome_type gaussian or binomial for continuous or binomial outcome
#' 
#' @export
#' 
#' @returns List of class `aggcomp_res` containing point estimates, standard error, and 95% confidence intervals for each effect of interest
aggcomp <- function(data,
                    laz_var_name,
                    abx_var_name,
                    site_var_name,
                    site_interaction,
                    age_var_name,
                    covariate_list, 
                    severity_list,
                    n_boot = 1000,
                    seed = 12345,
                    case_control = FALSE,
                    infection_var_name = NULL,
                    case_var_name = "case",
                    covariate_list_control = NULL,
                    outcome_type = "gaussian"){
  
  # set seed
  set.seed(seed)
  
  if(case_control == FALSE){
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
                               outcome_type = outcome_type)
    
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
    
    class(results) <- 'diarrhea_comp_res'
    
    parameters <- list(laz_var_name = laz_var_name,
                       abx_var_name = abx_var_name,
                       infection_var_name = infection_var_name,
                       site_var_name = site_var_name,
                       site_interaction = site_interaction,
                       age_var_name = age_var_name, 
                       covariate_list = covariate_list,
                       severity_list = severity_list)
    
  } else if (case_control == TRUE){
    pt_est <- abx_growth_gcomp_case_control(data = data,
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
    
    bootstrap_results <- bootstrap_estimates(data, 
                                             laz_var_name = laz_var_name,
                                             abx_var_name = abx_var_name,
                                             site_var_name = site_var_name,
                                             site_interaction = site_interaction,
                                             age_var_name = age_var_name,
                                             covariate_list = covariate_list, 
                                             severity_list = severity_list,
                                             n_boot = n_boot, 
                                             case_var_name = case_var_name,
                                             case_control = case_control,
                                             covariate_list_control = covariate_list_control, 
                                             outcome_type = outcome_type)
    
    # Return list of all point estimates and standard errors
    # Class `aggcomp_res`
    results <- c(unlist(pt_est), 
                 unlist(bootstrap_results))
    
    class(results) <- 'case_control_res'
    
    parameters <- list(laz_var_name = laz_var_name,
                       abx_var_name = abx_var_name,
                       case_var_name = case_var_name,
                       site_var_name = site_var_name,
                       site_interaction = site_interaction,
                       age_var_name = age_var_name, 
                       covariate_list = covariate_list,
                       covariate_list_control = covariate_list_control,
                       severity_list = severity_list)
    
  }
  
  aggcomp_results <- list(results = results,
                          parameters = parameters)
  
  class(aggcomp_results) <- "aggcomp_res"
  
  return(aggcomp_results)
  
}



