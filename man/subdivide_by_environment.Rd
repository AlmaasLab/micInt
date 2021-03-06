% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phyloseq_subset.R
\name{subdivide_by_environment}
\alias{subdivide_by_environment}
\title{Nest phyloseq objects according to sample variables}
\usage{
subdivide_by_environment(
  physeq,
  variables,
  keep_variables = TRUE,
  keep_empty = TRUE,
  extract_list = FALSE,
  ...
)
}
\arguments{
\item{physeq}{A \linkS4class{phyloseq} object}

\item{variables}{A character vector containing the names of the variables to nest by, these variables
should be discrete (factor, logical, integer, character)}

\item{keep_variables}{Logical, should the names of the sample variables being used, be kept
in the phyloseq objects?}

\item{keep_empty}{Logical, should sample data combinations applying to no samples be kept? If
\code{TRUE}, the coresponding \code{phyloseq} entity will be coded as \code{NULL}.}

\item{extract_list}{Logical, if \code{TRUE}, only the the \code{phyloseq} list column is returned}

\item{...}{Ignored}
}
\value{
A tibble where the variables are in the in the first columns, while the last column is named \code{phyloseq}
and contains the corresponding objects. If \code{extract_list=TRUE}, only this last column is returned.
\code{phyloseq} and contains the phyloseq object correspodning to the specific combination of sample variables.
}
\description{
This function tries to mimic \code{tidyr::nest}, subsetting the \code{phyloseq} object based
on all combinations of sample data.
}
\details{
If only one variable is provided, the element of the \code{phyloseq} column is named after this column
}
\examples{
library(phyloseq)
data("soilrep")
subdivide_by_environment(soilrep,variables=c('Treatment','warmed','clipped'),
keep_variables=FALSE,keep_empty=FALSE)
}
