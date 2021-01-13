#' @title An example dataset
#' @description A time series experiment of water samples from bioreactors with biofilm carriers
#'
#'
#' @format An experiment-level \code{\link{phyloseq}} object with 60 samples, 49 OTUs and
#' 3 sample variable:
#' \describe{
#' \item{Week}{Week of sampling, ranging from 1 to 12}
#' \item{Reactor}{The reactor replicate}
#' \item{Treatment}{Amount of biofilm carriers in the reactor.
#' 1: No carriers
#' 2: Intermediate number of carriers
#' 3: Many carriers
#' }
#'
#' }
"seawater"
