
#' Print the output of a \code{"aggcomp_res"} object
#' 
#' @param x An \code{"aggcomp_res"} object.
#' @param ... Other arguments (not used)
#'
#' @method print aggcomp_res
#' @export
print.aggcomp_res <- function(x, ...){
  
  res <- x$results
  
  if(class(res) == "case_control_res"){
    row_names <- c("E[Growth | Infection = 1, Abx = 0]",
                   "E[Growth | Infection = 1, Abx = 1]",
                   "E[Growth | Healthy Control]",
                   "Effect of case not treated with antibiotics vs healthy control",
                   "Effect of case treated with antibiotics vs healthy control")
    
    col_names <- c("Estimate", "Standard Error", "95% CI: Lower", "95% CI: Upper")
    
    tmp <- data.frame(
      c(res["abx_0_inf_1"],
        res["abx_1_inf_1"],
        res["abx_0_inf_0"],
        res["effect_inf_no_abx"],
        res["effect_inf_abx"]),
      c(res["se_abx_0_inf_1"],
        res["se_abx_1_inf_1"],
        res["se_abx_0_inf_0"],
        res["se_no_abx"],
        res["se_abx"]),
      c(res["lower_ci_abx_0_inf_1"],
        res["lower_ci_abx_1_inf_1"],
        res["lower_ci_abx_0_inf_0"],
        res["lower_ci_no_abx"],
        res["lower_ci_abx"]),
      c(res["upper_ci_abx_0_inf_1"],
        res["upper_ci_abx_1_inf_1"],
        res["upper_ci_abx_0_inf_0"],
        res["upper_ci_no_abx"],
        res["upper_ci_abx"])
    )
    
    colnames(tmp) <- col_names
    rownames(tmp) <- row_names
    
    # Print header with dashed line
    cat(paste("                                       Effect of ", x$parameters$case_var_name, "on ", x$parameters$laz_var_name, "in ", x$parameters$abx_var_name, "subgroups: Case-Control Study \n"))
    cat(paste(rep("-", 145), collapse = ""), "\n")
    cat(sprintf("%-70s%-20s%-20s%-20s%-20s\n", "", col_names[1], col_names[2], col_names[3], col_names[4]))
    cat(paste(rep("-", 145), collapse = ""), "\n")
    
    for(i in 1:nrow(tmp)){
      row_to_print <- tmp[i, ]
      
      # Adjust the widths as needed
      formatted_row <- sprintf("%-70s%-20s%-20s%-20s%-20s\n",
                               row.names(row_to_print),
                               round(row_to_print[1],4),
                               round(row_to_print[2],4),
                               round(row_to_print[3],4),
                               round(row_to_print[4],4))
      
      # Print the formatted row
      cat(paste(formatted_row))
      
    }
    
    cat(paste("\nNon-mediating covariates in cases: ", paste(x$parameters$covariate_list, collapse = ", ")))
    cat(paste("\nSeverity related covariates in cases: ", paste(x$parameters$severity_list, collapse = ", "), "\n"))
  
  } else{
  
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
    cat(sprintf("%-70s%-20s%-20s%-20s%-20s\n", "", col_names[1], col_names[2], col_names[3], col_names[4]))
    cat(paste(rep("-", 135), collapse = ""), "\n")
    
    for(i in 1:nrow(tmp)){
      row_to_print <- tmp[i, ]
      
      # Adjust the widths as needed
      formatted_row <- sprintf("%-70s%-20s%-20s%-20s%-20s\n",
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
    
  }
  
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
#' @param age_var_name name of age covariate in dataset (if applicable to add spline, else NULL)
#' @param severity_list character vector containing names of severity-related covariates (post-infection). If NULL, perform traditional gcomputation. Else, perform longitudinal gcomputation
#' @param covariate_list character vector containing names of baseline covariates
#' @param covariate_list_control character vector containing names of baseline covariates for controls (if applicable, else NULL)
#' @param id_var_name name of ID variable in dataset
#' 
#' @import fastDummies
#' 
#' @returns dataset with factors one-hot encoded
one_hot_encode <- function(data, 
                           laz_var_name,
                           abx_var_name,
                           infection_var_name,
                           site_var_name,
                           site_interaction,
                           age_var_name,
                           covariate_list, 
                           severity_list,
                           covariate_list_control = NULL,
                           id_var_name = "pid"){
  
  # Apply one-hot encoding to the relevant columns
  one_hot_data <- fastDummies::dummy_cols(data,
                                          select_columns = colnames(data)[
                                            !(colnames(data) %in% c(id_var_name, "child_id", 
                                                                    "sid", "first_id", "case_pid", 
                                                                    "case_id", "case_sid")) &
                                              (sapply(data, is.character) | sapply(data, is.factor))
                                          ],
                                          remove_first_dummy = FALSE,
                                          remove_selected_columns = TRUE,
                                          ignore_na = TRUE)
  
  # Remove spaces from column names in the one-hot encoded data
  colnames(one_hot_data) <- gsub(" ", "_", colnames(one_hot_data))
  
  ### Get new column names for covariates (covariate_list will have changed if any are factors)
  covariate_data <- data[,covariate_list, drop = FALSE]
  if(any(sapply(covariate_data, is.factor) == TRUE)){
    one_hot_covariate_data <- fastDummies::dummy_cols(covariate_data,
                                                      remove_first_dummy = FALSE,
                                                      remove_selected_columns = TRUE,
                                                      ignore_na = TRUE)
    colnames(one_hot_covariate_data) <- gsub(" ", "_", colnames(one_hot_covariate_data)) # remove spaces here as well
    one_hot_covariate_colnames <- colnames(one_hot_covariate_data)
  } else {
    one_hot_covariate_colnames <- colnames(covariate_data)
  }
  
  if(!is.null(covariate_list_control)){
    ### Get new column names for covariates (covariate_list will have changed if any are factors)
    covariate_data_control <- data[,covariate_list_control, drop = FALSE]
    if(any(sapply(covariate_data_control, is.factor) == TRUE)){
      one_hot_covariate_data_control <- fastDummies::dummy_cols(covariate_data_control,
                                                        remove_first_dummy = FALSE,
                                                        remove_selected_columns = TRUE,
                                                        ignore_na = TRUE)
      colnames(one_hot_covariate_data_control) <- gsub(" ", "_", colnames(one_hot_covariate_data_control)) # remove spaces here as well
      one_hot_covariate_colnames_control <- colnames(one_hot_covariate_data_control)
    } else {
      one_hot_covariate_colnames_control <- colnames(covariate_data_control)
    }
  } else {
    one_hot_covariate_colnames_control <- NULL
  }
  
  ### Get new column names for severity variables
  severity_data <- data[,severity_list, drop = FALSE]
  if(any(sapply(severity_data, is.factor) == TRUE)){
    one_hot_severity_data <- fastDummies::dummy_cols(severity_data,
                                                     remove_first_dummy = FALSE,
                                                     remove_selected_columns = TRUE,
                                                     ignore_na = TRUE)
    colnames(one_hot_severity_data) <- gsub(" ", "_", colnames(one_hot_severity_data)) # remove spaces here as well
    one_hot_severity_colnames <- colnames(one_hot_severity_data)
  } else {
    one_hot_severity_colnames <- colnames(severity_data)
  }
  
  # If there is an interaction with site, get site names in one-hot encoded data to be used later to add interaction term
  if(site_interaction == "TRUE"){
    site_data <- data[,site_var_name, drop = FALSE]
    one_hot_site_data <- fastDummies::dummy_cols(site_data,
                                                 remove_first_dummy = FALSE,
                                                 remove_selected_columns = TRUE,
                                                 ignore_na = TRUE)
    colnames(one_hot_site_data) <- gsub(" ", "_", colnames(one_hot_site_data)) # remove spaces here as well
    one_hot_enroll_site_colnames <- colnames(one_hot_site_data)
  } else {
    one_hot_enroll_site_colnames <- NULL
  }
  
  return(list(data = one_hot_data,
              covariate_list = one_hot_covariate_colnames,
              covariate_list_control = one_hot_covariate_colnames_control,
              severity_list = one_hot_severity_colnames,
              site_var_names = one_hot_enroll_site_colnames))
}


#' Function to make and print descriptive tables to aid in interpreting gcomp
#' 
#' @param data dataset used in analysis
#' @param laz_var_name name of growth outcome variable
#' @param abx_var_name name of binary antibiotic variable
#' @param infection_var_name name of binary infection variable
#' @param site_var_name name of site covariate in dataset (if applicable, else null)
#' @param site_interaction TRUE or FALSE indicating interaction between site and antibiotics 
#' @param age_var_name name of age covariate in dataset (if applicable, else NULL)
#' @param severity_list character vector containing names of severity-related covariates (post-infection). If NULL, perform traditional gcomputation. Else, perform longitudinal gcomputation
#' @param covariate_list character vector containing names of baseline covariates
#' 
#' @import gtsummary
#' @import table1
#' @import dplyr
#' @import labelled
#' 
#' @returns list of descriptive tables
make_tables <- function(data,
                        laz_var_name,
                        abx_var_name,
                        infection_var_name,
                        site_var_name,
                        site_interaction,
                        age_var_name,
                        severity_list, 
                        covariate_list){
  
  # if age in covariate_list, change to spline with 3 knots to increase flexibility
  if(age_var_name %in% covariate_list){
    covariate_list <- covariate_list[covariate_list != age_var_name]
    covariate_list <- c(covariate_list, paste0("splines::ns(", age_var_name, ", df = 4)"))
  }
  
  if(!is.null(severity_list)){
    
    # -----------------------------------
    # Longitudinal G-Computation 
    # -----------------------------------
    
    # Get model formula 
    if(site_interaction == "TRUE"){
      
      # if site in covariate list, remove for model1 (only want interaction terms)
      if(site_var_name %in% covariate_list){
        covariate_list <- covariate_list[covariate_list != site_var_name]
      }
      
      model_1_formula <- stats::as.formula(paste(laz_var_name, "~", 
                                                 abx_var_name, "+", infection_var_name, "+", 
                                                 paste(covariate_list, collapse = "+"), "+",
                                                 paste(severity_list, collapse = "+"), "+",
                                                 paste(site_var_name, "*", abx_var_name, collapse = "+"), "+",
                                                 paste0(infection_var_name, "*", abx_var_name)))
      
      # put site back in covariate list
      if(!is.null(site_var_name)){
        covariate_list <- c(covariate_list, site_var_name)
      }
      
    } else {
      
      # if site was listed but is not already in covariate list, add
      if(!is.null(site_var_name) & (!(site_var_name %in% covariate_list))){
        covariate_list <- c(covariate_list, site_var_name)
      }
      
      model_1_formula <- stats::as.formula(paste(laz_var_name, "~", 
                                                 abx_var_name, "+", infection_var_name, "+", 
                                                 paste(covariate_list, collapse = "+"), "+",
                                                 paste(severity_list, collapse = "+"), "+",
                                                 paste0(infection_var_name, "*", abx_var_name)))
    }
    
    # Step 1: regress LAZ on abx, infection, all severity / non-mediating variables
    model_1 <- stats::glm(model_1_formula,
                          data = data, 
                          family = "gaussian")
    
    caption <- paste0("Summary of model for ", labelled::var_label(model_1$y), " ~ ", labelled::var_label(data[abx_var_name]), 
                      " + \n\n ", labelled::var_label(data[infection_var_name]), " + Severity Covariates + Baseline Covariates")
    
    
    # Summary of logistic regression of abx ~ shigella + severity + BL covariates
    abx_logistic <- stats::glm(stats::as.formula(paste(abx_var_name, "~", 
                                         infection_var_name, "+",
                                         paste(covariate_list, collapse = "+"), "+", 
                                         paste(severity_list, collapse = "+"))),
                        data = data,
                        family = "binomial")
    
    caption2 <- paste0("Summary of model for ", labelled::var_label(abx_logistic$y), " ~ ", labelled::var_label(data[infection_var_name]), 
                       " + Severity Covariates + Baseline Covariates")
    
    
  } else {
    
    # -----------------------------------
    # Traditional G-Computation (no severity covariates)
    # -----------------------------------
    if(site_interaction == "TRUE"){
      model_1_formula <- stats::as.formula(paste(laz_var_name, "~",
                                                 abx_var_name, "+", infection_var_name, "+",
                                                 paste(covariate_list, collapse = "+"), "+",
                                                 paste0(site_var_name, "*", abx_var_name), "+",
                                                 paste0(infection_var_name, "*", abx_var_name)))
    } else {
      model_1_formula <- stats::as.formula(paste(laz_var_name, "~",
                                                 abx_var_name, "+", infection_var_name, "+",
                                                 paste(covariate_list, collapse = "+"), "+",
                                                 paste0(infection_var_name, "*", abx_var_name)))
    }
    
    # Fit model with LAZ ~ infection + abx + covariates + infection*abx
    model_1 <- stats::glm(model_1_formula,
                          data = data,
                          family = "gaussian")
    
    caption <- paste0("Summary of model for ", labelled::var_label(model_1$y), " ~ ", labelled::var_label(data[abx_var_name]), 
                      " + \n\n ", labelled::var_label(data[infection_var_name]), " + Baseline Covariates")
    
    # Summary of logistic regression of abx ~ shigella + severity + BL covariates
    abx_logistic <- stats::glm(stats::as.formula(paste(abx_var_name, "~", 
                                         infection_var_name, "+",
                                         paste(covariate_list, collapse = "+"))),
                        data = data,
                        family = "binomial")
    
    caption2 <- paste0("Summary of model for ", labelled::var_label(abx_logistic$y), " ~ ", labelled::var_label(data[infection_var_name]), 
                       " + Baseline Covariates")
  }
    
  tbl_1 <- model_1 %>%
    tbl_regression(estimate_fun = label_style_sigfig(digits = 3),
                   pvalue_fun = label_style_pvalue(digits = 3)) %>%
    modify_caption(caption)
  
  tbl_2 <- abx_logistic %>%
    tbl_regression(exponentiate = TRUE,
                   estimate_fun = label_style_sigfig(digits = 3),
                   pvalue_fun = label_style_pvalue(digits = 3)) %>%
    modify_caption(caption2)
  
  # Summary of logistic regression of shigella ~ BL covariates
  shig_bl_logistic <- stats::glm(stats::as.formula(paste(infection_var_name, "~",
                                           paste(covariate_list, collapse = "+"))),
                          data = data,
                          family = "binomial")
  
  caption3 <- paste0("Summary of model for ", labelled::var_label(shig_bl_logistic$y), " ~ Baseline Covariates")
  
  tbl_3 <- shig_bl_logistic %>% 
    tbl_regression(exponentiate = TRUE,
                   estimate_fun = label_style_sigfig(digits = 3),
                   pvalue_fun = label_style_pvalue(digits = 3)) %>%
    modify_caption(caption3)
  
  # ADDITIONAL DESCRIPTIVE TABLES
  if(!is.null(severity_list)){
    # Summary of appropriate regressions for each severity measure (severity ~ shigella + baseline)
    severity_tbls <- vector("list",length = length(severity_list))
    
    for(i in 1:length(severity_list)){
      if(is.factor(data[,severity_list[i]])){
        severity_fit <- stats::glm(stats::as.formula(paste(severity_list[i], "~", 
                                             infection_var_name, "+",
                                             paste(covariate_list, collapse = "+"))),
                            data = data,
                            family = "binomial")
   
        caption4 <- paste0("Summary of model for ", labelled::var_label(data[severity_list[i]]), " ~ \n\n ", labelled::var_label(data[infection_var_name]), " + Baseline Covariates")
        
        
        severity_tbls[[i]] <- severity_fit %>%
          tbl_regression(exponentiate = TRUE,
                         estimate_fun = label_style_sigfig(digits = 3),
                         pvalue_fun = label_style_pvalue(digits = 3)) %>%
          modify_caption(caption4)
        
      } else{
        severity_fit <- stats::glm(stats::as.formula(paste(severity_list[i], "~", 
                                             infection_var_name, "+",
                                             paste(covariate_list, collapse = "+"))),
                            data = data,
                            family = "gaussian")
        
        
        caption4 <- paste0("Summary of model for ", labelled::var_label(data[severity_list[i]]), " ~ \n\n", labelled::var_label(data[infection_var_name]), " + Baseline Covariates")
        
        
        severity_tbls[[i]] <- severity_fit %>%
          tbl_regression(estimate_fun = label_style_sigfig(digits = 3),
                         pvalue_fun = label_style_pvalue(digits = 3)) %>%
          modify_caption(caption4)
        
      }
    }
    
    # Descriptive statistics of severity measures and baseline covariates by shigella status
    # Relabel the grouping variable
    name_shig <- labelled::var_label(data[[infection_var_name]])
    data[[infection_var_name]] <- factor(data[[infection_var_name]],
                                        levels = c(0, 1),
                                        labels = c("No Shigella", "Shigella"))
    
    # Helper function for table1
    my.render.cont <- function(x) {
      stats <- stats::quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
      sprintf("%s (%s, %s)", round(stats[2], 2), round(stats[1], 2), round(stats[3], 2))
    }
    
   
    caption5 <- paste("Descriptive statistics of severity measures and baseline covariates by ", name_shig)
   
    # if missing any infection_var_name, drop
    data <- data %>%
      dplyr::filter(!is.na(.[[infection_var_name]]))
    
    descriptive_tbl <- table1(stats::as.formula(paste("~", paste(severity_list, collapse = "+"), "+",
                                               paste(covariate_list, collapse = "+"), 
                                               "|" , infection_var_name)),
                              data = data,
                              render.continuous = my.render.cont,
                              footnote = "Median (IQR) or n (%)",
                              overall = FALSE,
                              caption = caption5)
  } else {
    
    # No severity tables
    severity_tbls <- NULL
    
    # Descriptive statistics of severity measures and baseline covariates by shigella status
    # Relabel the grouping variable
    name_shig <- labelled::var_label(data[[infection_var_name]])
    data[[infection_var_name]] <- factor(data[[infection_var_name]],
                                        levels = c(0, 1),
                                        labels = c("No Shigella", "Shigella"))
    
    # Helper function for table1
    my.render.cont <- function(x) {
      stats <- stats::quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
      sprintf("%s (%s, %s)", round(stats[2], 2), round(stats[1], 2), round(stats[3], 2))
    }
    
    caption5 <- paste("Descriptive statistics of baseline covariates by ", name_shig)
    
    # if missing any infection_var_name, drop
    data <- data %>%
      dplyr::filter(!is.na(.[[infection_var_name]]))
    
    descriptive_tbl <- table1(stats::as.formula(paste("~", paste(covariate_list, collapse = "+"), 
                                               "|" , infection_var_name)),
                              data = data,
                              render.continuous = my.render.cont,
                              footnote = "Median (IQR) or n (%)",
                              overall = FALSE,
                              caption = caption5)
  }
  
  return(list(tbl_1 = tbl_1,
              tbl_2 = tbl_2,
              tbl_3 = tbl_3,
              severity_tbls = severity_tbls,
              descriptive_tbl = descriptive_tbl))
}

#' Function to make confidence interval plots of results
#' 
#' @param results `aggcomp_res` object or list of `aggcomp_res` objects
#' @param abx_label label or list of labels to use for antibiotics variable (same length as results list)
#' @param inf_label label or list of labels to use for infection variable (same length as results list)
#' @param outcome_label label or list of labels to use for outcome variable (same length as results list)
#' 
#' @import ggplot2
#' 
#' @return list of confidence interval ggplots
plot_cis <- function(results, 
                     abx_labels = NULL,
                     inf_labels = NULL, 
                     outcome_labels = NULL){
  
  # if single object, put into list for lapply
  if(!is(results, "list")){
    results <- list(results)
  }
  
  plot_results <- lapply(1:length(results), function(i){
    
    res <- results[[i]]
    
    if(is.null(abx_labels)){
      abx_label <- res$parameters$abx_var_name
    } else {
      abx_label <- abx_labels[[i]]
    }
    
    if(is.null(inf_labels)) {
      inf_label <- res$parameters$infection_var_name
    } else{
      inf_label <- inf_labels[[i]]
    }
    
    if(is.null(outcome_labels)) { 
      outcome_label <- res$parameters$laz_var_name
    } else {
      outcome_label <- outcome_labels[[i]]
    }
    
    plot_label <- paste0("Antibiotics: ", abx_label, ", Infection: ", inf_label) 
    
    # make dataframe for plotting
    results_vec <- res$results
    plot_data <- data.frame(
      subgroup = c("Antibiotics", "No Antibiotics"),
      pt_est = c(results_vec['effect_inf_abx'], results_vec['effect_inf_no_abx']),
      lower_ci = c(results_vec['lower_ci_abx'], results_vec['lower_ci_no_abx']),
      upper_ci = c(results_vec['upper_ci_abx'], results_vec['upper_ci_no_abx']),
      label = plot_label
    )
    
    x_loc <- max(plot_data$upper_ci) + 0.01
    
    ggplot2::ggplot(data = plot_data, aes(x = pt_est, y = subgroup, color = subgroup)) +
      geom_point(size = 5) +  # Increase point size
      geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.3) +  # Increase error bar width
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      geom_text(aes(x = x_loc , label = paste0("Point estimate: ", round(pt_est, 4), "\n95% CI: [", round(lower_ci, 4), ", ", round(upper_ci, 4), "]")), 
                hjust = 0, vjust = 0.5, size = 5) +  # Increase text size
      labs(title = paste0("Estimated effect of ", inf_label, " infection on ", outcome_label),
           x = "Estimated effect", y = "Treatment subgroup") +
      theme_minimal(base_size = 14) +  # Increase base font size for all text
      facet_wrap(~ label, scales = "free_y", ncol = 1) +
      labs(color = "Treatment subgroup") +
      theme(
        strip.text = element_text(size = 16),  # Increase facet label size
        axis.title = element_text(size = 18),  # Increase axis titles
        axis.text = element_text(size = 14),  # Increase axis text size
        plot.title = element_text(size = 20)  # Increase plot title size
      ) +
      coord_cartesian(clip = "off")  # Prevent clipping of text labels
  })
  
  return(plot_results)
}


#' Function to make confidence interval plots of results with multiple subplots
#' 
#' @param results `aggcomp_res` object or list of `aggcomp_res` objects
#' @param abx_labels label or list of labels to use for antibiotics variable (same length as results list)
#' @param inf_labels label or list of labels to use for infection variable (same length as results list)
#' @param outcome_labels label or list of labels to use for outcome variable (same length as results list)
#' 
#' @import ggplot2
#' 
#' @return list of confidence interval ggplots
plot_cis_subplots <- function(results, 
                     abx_labels = NULL,
                     inf_labels = NULL, 
                     outcome_labels = NULL) {
  
  # if single object, put into list for lapply
  if (!is(results, "list")) {
    results <- list(results)
  }
  
  # Determine the global x-axis limits
  x_lims <- lapply(results, function(res) {
    results_vec <- res$results
    range(c(results_vec['lower_ci_abx'], results_vec['upper_ci_abx'],
            results_vec['lower_ci_no_abx'], results_vec['upper_ci_no_abx']))
  })
  global_x_min <- min(sapply(x_lims, `[`, 1))
  global_x_max <- max(sapply(x_lims, `[`, 2))
  
  # Add padding for annotations
  annotation_padding <- (global_x_max - global_x_min) * 0.5
  global_x_max <- global_x_max + annotation_padding
  
  # Generate individual plots
  plot_results <- lapply(1:length(results), function(i) {
    res <- results[[i]]
    
    if (is.null(abx_labels)) {
      abx_label <- res$parameters$abx_var_name
    } else {
      abx_label <- abx_labels[[i]]
    }
    
    if (is.null(inf_labels)) {
      inf_label <- res$parameters$infection_var_name
    } else {
      inf_label <- inf_labels[[i]]
    }
    
    if (is.null(outcome_labels)) { 
      outcome_label <- res$parameters$laz_var_name
    } else {
      outcome_label <- outcome_labels[[i]]
    }
    
    plot_label <- paste0("Antibiotics: ", abx_label, ", Infection: ", inf_label) 
    
    # Create dataframe for plotting
    results_vec <- res$results
    plot_data <- data.frame(
      subgroup = c("Antibiotics", "No Antibiotics"),
      pt_est = c(results_vec['effect_inf_abx'], results_vec['effect_inf_no_abx']),
      lower_ci = c(results_vec['lower_ci_abx'], results_vec['lower_ci_no_abx']),
      upper_ci = c(results_vec['upper_ci_abx'], results_vec['upper_ci_no_abx']),
      label = plot_label
    )
    
    x_loc <- global_x_max - annotation_padding / 2
    
    ggplot2::ggplot(data = plot_data, aes(x = pt_est, y = subgroup, color = subgroup)) +
      geom_point(size = 5) +
      geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.3) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      geom_text(aes(x = x_loc, 
                    label = paste0("Point estimate: ", round(pt_est, 4), 
                                   "\n95% CI: [", round(lower_ci, 4), ", ", round(upper_ci, 4), "]")), 
                hjust = 0, vjust = 0.5, size = 5) +
      labs(title = plot_label,
           x = "Estimated effect", y = "Treatment subgroup") +
      theme_minimal(base_size = 14) +
      labs(color = "Treatment subgroup") +
      theme(
        plot.title = element_text(size = 14, hjust = 0.5),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)
      ) +
      coord_cartesian(xlim = c(global_x_min, global_x_max))
  })
  
  # Combine plots with patchwork, keeping individual plot titles and shared x-axis
  combined_plot <- patchwork::wrap_plots(plot_results, ncol = 1) +
    plot_layout(guides = "collect") & 
    theme(legend.position = "bottom")
  
  return(combined_plot)
}

#' Helper function to impute select covariates with mean or mode
#' 
#' @param data dataset to impute
#' @param imp_covariates list of covariates to be imputed
#' @param site_var_name name of site variable if imputing by site
#' @param imp_by_site boolean indicating if imputing by site, default true
#' @param mode_for_binary use mode imputation for binary 0,1 variables
#' 
#' @export
#' 
#' @return dataframe with imputed data and flag columns indicating imputed values
impute_covariates <- function(data,
                              imp_covariates,
                              site_var_name = "site",
                              imp_by_site = TRUE,
                              mode_for_binary = FALSE){
  
  mode_fn <- function(x) {
    x <- na.omit(x)
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  for(cov_name in imp_covariates){
    
    imputed_flag_name <- paste0(cov_name, "_imputed")
    data[[imputed_flag_name]] <- ifelse(is.na(data[[cov_name]]), 1, 0)
    
    if(imp_by_site == TRUE){
      
      if(is.factor(data[[cov_name]]) | (all(unique(data[[cov_name]]) %in% c(0,1, NA)) & mode_for_binary)){
        # Binary var or factor 
        
        # Get mode for each site
        site_modes <- data %>%
          group_by(!!sym(site_var_name)) %>%
          summarise(mode = mode_fn(!!sym(cov_name)))
        
        data <- data %>%
          left_join(site_modes, by = site_var_name) %>%
          mutate(!!cov_name := if_else(is.na(!!sym(cov_name)), mode, !!sym(cov_name))) %>%
          select(-mode)
        
      } else {
        # Continuous var
        site_means <- data %>%
          group_by(!!sym(site_var_name)) %>%
          summarise(site_mean = mean(!!sym(cov_name), na.rm = TRUE))
        
        # replace with site specific mean
        data <- data %>%
          left_join(site_means, by = site_var_name) %>%
          mutate(!!cov_name := if_else(is.na(!!sym(cov_name)), site_mean, !!sym(cov_name))) %>%
          select(-site_mean)
        
      }
      
    } else{

     exit("implement non- site-specific version later")
      
    }
    
  }
  
  return(data)
  
}
  
