% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{plot_cis}
\alias{plot_cis}
\title{Function to make confidence interval plots of results}
\usage{
plot_cis(results, abx_labels = NULL, inf_labels = NULL, outcome_labels = NULL)
}
\arguments{
\item{results}{\code{aggcomp_res} object or \code{agaipw_res} object, or list of \code{aggcomp_res} and/or \code{agaipw_res} objects}

\item{abx_label}{label or list of labels to use for antibiotics variable (same length as results list)}

\item{inf_label}{label or list of labels to use for infection variable (same length as results list)}

\item{outcome_label}{label or list of labels to use for outcome variable (same length as results list)}
}
\value{
list of confidence interval ggplots
}
\description{
Function to make confidence interval plots of results
}
