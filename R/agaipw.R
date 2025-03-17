#' Main function to get point estimates and confidence intervals using AIPW
#' 
#' Note parameters include all options across other diarrhea and case control analyses
#' 
#' @param data dataframe containing dataset used for gcomp
#' @param laz_var_name name of growth outcome variable
#' @param abx_var_name name of binary antibiotic variable
#' @param infection_var_name name of binary infection variable (for other diarrhea analyses; case_control = FALSE)
#' @param case_var_name name of variable indicating case (for case control analysis; case_control = TRUE)
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
#' @param sl.library.missingness.case list containing SuperLearner libraries to use for outcome missingness model (cases)
#' @param sl.library.missingness.control list containing SuperLearner libraries to use for outcome missingness model (controls)
#' @param case_control TRUE if case control analysis, FALSE for other diarrhea analysis. Default to FALSE
#' @param seed seed for reproducibility in SuperLearner cross-validation
#' @param v_folds number of cross-validation folds to use in SuperLearner
#' @param return_models boolean return SuperLearner models. Default FALSE.
#' @param child_id_var_name name of variable indicating child_id in data. Used to account for re-enrollment in study. 
#' 
#' @export
#' 
#' @returns List of class `agaipw_res` containing objects for AIPW result, AIPW models.
agaipw <- function(data,
                    laz_var_name,
                    abx_var_name,
                    infection_var_name = NA,
                    case_var_name = NA,
                    covariate_list, 
                    severity_list = NULL,
                    pathogen_quantity_list = NULL,
                    pathogen_attributable_list = NULL,
                    no_etiology_var_name = NULL,
                    outcome_type = "gaussian",
                    sl.library.outcome = c("SL.glm"),
                    sl.library.outcome.2 = c("SL.glm"),
                    sl.library.outcome.case = c("SL.glm"),
                    sl.library.outcome.control = c("SL.glm"),
                    sl.library.treatment = c("SL.mean"),
                    sl.library.infection = c("SL.glm"),
                    sl.library.missingness = c("SL.glm"),
                    sl.library.missingness.case = c("SL.glm"),
                    sl.library.missingness.control = c("SL.glm"),
                    case_control = FALSE,
                    seed = 12345,
                    v_folds = 5,
                    return_models = FALSE,
                    child_id_var_name = NULL){
  
  # Set seed for reproducibility
  set.seed(seed)
  
  if(case_control == FALSE){
    
    ## Other diarrhea analysis
    if(!is.null(severity_list)){
      # Include second stage outcome regression
      aipw_est <- aipw_other_diarrhea_2(data = data,
                                        laz_var_name = laz_var_name,
                                        abx_var_name = abx_var_name,
                                        infection_var_name = infection_var_name,
                                        covariate_list = covariate_list,
                                        severity_list = severity_list,
                                        pathogen_quantity_list = pathogen_quantity_list,
                                        pathogen_attributable_list = pathogen_attributable_list,
                                        no_etiology_var_name = no_etiology_var_name,
                                        outcome_type = outcome_type,
                                        sl.library.outcome = sl.library.outcome,
                                        sl.library.outcome.2 = sl.library.outcome.2,
                                        sl.library.treatment = sl.library.treatment,
                                        sl.library.infection = sl.library.infection,
                                        sl.library.missingness = sl.library.missingness,
                                        v_folds = v_folds,
                                        return_models = return_models,
                                        child_id_var_name = child_id_var_name)
    } else {
      # Do not include second stage outcome regression
      aipw_est <- aipw_other_diarrhea(data = data,
                                      laz_var_name = laz_var_name,
                                      abx_var_name = abx_var_name,
                                      infection_var_name = infection_var_name,
                                      covariate_list = covariate_list,
                                      pathogen_quantity_list = pathogen_quantity_list,
                                      pathogen_attributable_list = pathogen_attributable_list,
                                      no_etiology_var_name = no_etiology_var_name,
                                      outcome_type = outcome_type,
                                      sl.library.outcome = sl.library.outcome,
                                      sl.library.treatment = sl.library.treatment,
                                      sl.library.infection = sl.library.infection,
                                      sl.library.missingness = sl.library.missingness,
                                      v_folds = v_folds,
                                      return_models = return_models,
                                      child_id_var_name = child_id_var_name)
    }
    
  } else{
    
    ## Case-control analysis
    aipw_est <- aipw_case_control(data = data,
                                  laz_var_name = laz_var_name,
                                  abx_var_name = abx_var_name,
                                  case_var_name = case_var_name,
                                  covariate_list = covariate_list,
                                  severity_list = severity_list,
                                  pathogen_quantity_list = pathogen_quantity_list,
                                  outcome_type = outcome_type,
                                  sl.library.outcome.case = sl.library.outcome.case,
                                  sl.library.outcome.control = sl.library.outcome.control,
                                  sl.library.treatment = sl.library.treatment,
                                  sl.library.infection = sl.library.infection,
                                  sl.library.missingness.case = sl.library.missingness.case,
                                  sl.library.missingness.control = sl.library.missingness.control,
                                  v_folds = v_folds,
                                  return_models = return_models,
                                  child_id_var_name = child_id_var_name)
    
  }
  
  class(aipw_est) <- "agaipw_res"
  
  return(aipw_est)
  
}

