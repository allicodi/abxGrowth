% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{one_hot_encode}
\alias{one_hot_encode}
\title{Helper function for one-hot encoding dataset}
\usage{
one_hot_encode(
  data,
  laz_var_name,
  abx_var_name,
  infection_var_name,
  site_var_name,
  site_interaction,
  age_var_name,
  covariate_list,
  severity_list,
  covariate_list_control = NULL,
  id_var_name = "pid"
)
}
\arguments{
\item{data}{dataframe containing dataset used for gcomp}

\item{laz_var_name}{name of growth outcome variable}

\item{abx_var_name}{name of binary antibiotic variable}

\item{infection_var_name}{name of binary infection variable}

\item{site_var_name}{name of site variable in dataset}

\item{site_interaction}{TRUE or FALSE indicating interaction between site and antibiotics}

\item{age_var_name}{name of age covariate in dataset (if applicable to add spline, else NULL)}

\item{covariate_list}{character vector containing names of baseline covariates}

\item{severity_list}{character vector containing names of severity-related covariates (post-infection). If NULL, perform traditional gcomputation. Else, perform longitudinal gcomputation}

\item{covariate_list_control}{character vector containing names of baseline covariates for controls (if applicable, else NULL)}

\item{id_var_name}{name of ID variable in dataset}
}
\value{
dataset with factors one-hot encoded
}
\description{
Helper function for one-hot encoding dataset
}
