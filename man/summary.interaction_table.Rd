% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_refinement.R
\name{summary.interaction_table}
\alias{summary.interaction_table}
\title{Create summary of \code{interaction_table}}
\usage{
\method{summary}{interaction_table}(table)
}
\arguments{
\item{table}{An interaction table returned from \link{create_interaction_table}}
}
\value{
A list with the following fields (the first three are taken directly from
the \code{interaction.table}):
\itemize{
\item \code{measure_name}
\item \code{signed}
\item \code{measure_type}
\item \code{number_significant}: The number of significant interactions
\item \code{proportion_negative}: The proprotion of significant interactions which
are negative
}
}