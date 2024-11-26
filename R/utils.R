
#' Print the output of a \code{"aggcomp_res"} object
#' 
#' @param x An \code{"aggcomp_res"} object.
#' @param ... Other arguments (not used)
#'
#' @method print aggcomp_res
#' @export
print.aggcomp_res <- function(x, ...){
  
  res <- x$results
  
  row_names <- c("E[Growth | Infection = 1, Abx = 0]",
                 "E[Growth | Infection = 0, Abx = 0]",
                 "Effect of infection in no antibiotics subgroup",
                 "E[Growth | Infection = 1, Abx = 1]",
                 "E[Growth | Infection = 0, Abx = 1]", 
                 "Effect of infection in antibiotics subgroup")
  
  col_names <- c("Estimate", "Standard Error", "95% CI: Lower", "95% CI: Upper")
  
  tmp <- data.frame(
    c(res["abx_0_inf_1"],
      res["abx_0_inf_0"],
      res["effect_inf_no_abx"],
      res["abx_1_inf_1"],
      res["abx_1_inf_0"],
      res["effect_inf_abx"]),
    c(res["se_abx_0_inf_1"],
      res["se_abx_0_inf_0"],
      res["se_no_abx"],
      res["se_abx_1_inf_1"],
      res["se_abx_1_inf_0"],
      res["se_abx"]),
    c(res["lower_ci_abx_0_inf_1"],
      res["lower_ci_abx_0_inf_0"],
      res["lower_ci_no_abx"],
      res["lower_ci_abx_1_inf_1"],
      res["lower_ci_abx_1_inf_0"],
      res["lower_ci_abx"]),
    c(res["upper_ci_abx_0_inf_1"],
      res["upper_ci_abx_0_inf_0"],
      res["upper_ci_no_abx"],
      res["upper_ci_abx_1_inf_1"],
      res["upper_ci_abx_1_inf_0"],
      res["upper_ci_abx"])
  )
  
  colnames(tmp) <- col_names
  rownames(tmp) <- row_names
  
  # Print header with dashed line
  cat(paste("                           Effect of ", x$parameters$infection_var_name, "on ", x$parameters$laz_var_name, "in ", x$parameters$abx_var_name, "subgroups \n"))
  cat(paste(rep("-", 135), collapse = ""), "\n")
  cat(sprintf("%-60s%-20s%-20s%-20s%-20s\n", "", col_names[1], col_names[2], col_names[3], col_names[4]))
  cat(paste(rep("-", 135), collapse = ""), "\n")
  
  for(i in 1:nrow(tmp)){
    row_to_print <- tmp[i, ]
    
    # Adjust the widths as needed
    formatted_row <- sprintf("%-60s%-20s%-20s%-20s%-20s\n",
                             row.names(row_to_print),
                             round(row_to_print[1],4),
                             round(row_to_print[2],4),
                             round(row_to_print[3],4),
                             round(row_to_print[4],4))
    
    # Print the formatted row
    cat(paste(formatted_row))
    
  }
  
  cat(paste("\nNon-mediating covariates: ", paste(x$parameters$covariate_list, collapse = ", ")))
  cat(paste("\nSeverity related covariates: ", paste(x$parameters$severity_list, collapse = ", "), "\n"))
  
  invisible(tmp)
  
}

#' Helper function to combine multiple "aggcomp_res" objects into a CSV file of results
#' 
#' @param res_list list of "aggcomp_res" objects
combine_res <- function(res_list) {
  
  # Process each object in the list
  res_df <- lapply(res_list, function(row) {
    main_res <- row$results
    params <- row$parameters
    
    # Convert covariate_list and severity_list to character vectors
    params$covariate_list <- paste(params$covariate_list, collapse = ", ")
    params$severity_list <- paste(params$severity_list, collapse = ", ")
    
    # Combine results and parameters into a single row
    c(main_res, params)
  })
  
  # Combine all rows into a dataframe
  res_df <- do.call(rbind, res_df)
  
  # Convert to dataframe with proper column names
  res_df <- as.data.frame(res_df, stringsAsFactors = FALSE)
  
  return(res_df)
}


#' Helper function for one-hot encoding dataset
#' 
#' @param data dataframe containing dataset used for gcomp
#' @param laz_var_name name of growth outcome variable
#' @param abx_var_name name of binary antibiotic variable
#' @param infection_var_name name of binary infection variable
#' @param site_var_name name of site variable in dataset
#' @param site_interaction TRUE or FALSE indicating interaction between site and antibiotics 
#' @param age_name name of age covariate in dataset (if applicable to add spline, else NULL)
#' @param severity_list character vector containing names of severity-related covariates (post-infection). If NULL, perform traditional gcomputation. Else, perform longitudinal gcomputation
#' @param covariate_list character vector containing names of baseline covariates
#' 
#' @returns dataset with factors one-hot encoded
one_hot_encode <- function(data, 
                           laz_var_name,
                           abx_var_name,
                           infection_var_name,
                           site_var_name,
                           site_interaction,
                           age_name,
                           covariate_list, 
                           severity_list){
  
  ### "One hot encode" factors in dataset
  
  # outcome var has missing data, remove that and add back on after to avoid losing rows in data
  outcome <- data[[laz_var_name]]
  tmp_data <- data[,colnames(data) != laz_var_name]
  
  # set contrasts.arg for factors so all levels are included vs losing one to reference (causes issues if trying to stratify by reference variable)
  contrasts.arg.all <- lapply(tmp_data[, sapply(tmp_data, is.factor), drop = FALSE], stats::contrasts, contrasts = FALSE)
  
  # get one hot encoded full data
  one_hot_data <- data.frame(stats::model.matrix(~ . -1, data = tmp_data, contrasts.arg = contrasts.arg.all))
  
  # add outcome var back
  one_hot_data[[laz_var_name]] <- outcome
  
  ### Get new column names for covariates (covariate_list will have changed if any are factors)
  
  covariate_data <- data[,covariate_list, drop = FALSE]
  if(any(sapply(covariate_data, is.factor) == TRUE)){
    one_hot_covariate_colnames <- colnames(data.frame(stats::model.matrix(~ . - 1, data = covariate_data, contrasts.arg = contrasts.arg.all)))
  } else {
    one_hot_covariate_colnames <- colnames(covariate_data)
  }
  
  ### Get new column names for severity variables
  severity_data <- data[,severity_list, drop = FALSE]
  if(any(sapply(severity_data, is.factor) == TRUE)){
    one_hot_severity_colnames <- colnames(data.frame(stats::model.matrix(~ . - 1, data = severity_data, contrasts.arg = contrasts.arg.all)))
  } else {
    one_hot_severity_colnames <- colnames(severity_data)
  }
  
  # If there is an interaction with site, get site names in one-hot encoded data to be used later to add interaction term
  if(site_interaction == "TRUE"){
    one_hot_enroll_site_colnames <- colnames(data.frame(stats::model.matrix(~ . - 1, data = data[,site_var_name, drop = FALSE])))
  } else {
    one_hot_enroll_site_colnames <- NULL
  }
  
  return(list(data = one_hot_data,
              covariate_list = one_hot_covariate_colnames,
              severity_list = one_hot_severity_colnames,
              site_var_names = one_hot_enroll_site_colnames))
  
}
