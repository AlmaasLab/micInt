#' @title OTU statistics
#'
#' @description
#' This function returns a summary of OTU statistics for the OTUs of
#' an OTU table
#'
#' @param refined_table The refined OTU table returned from \link{refine_data} or a \code{phyloseq} object.
#' @inheritParams  output_ccrepe_data
#'
#' @return
#' A \code{data.frame} with one row for each OTU with the following
#' elements:
#'
#' \itemize{
#' \item \code{ID} The OTU name
#' \item \code{meanAbundance}
#' \item \code{medianAbundance}
#' \item \code{maxAbundance}
#' \item \code{numberNonZero} The number of samples where the OTU has a
#' non-zero abundance
#' \item \code{taxonomy} collapsed into a single string, may be missing
#' }
#' @export
OTU_stats <- function(refined_table, taxonomy = NULL) {
  if (inherits(refined_table, "phyloseq")) {
    refined_data <- refined_table %>% phyloseq::otu_table() %>% as.matrix()
    if (phyloseq::taxa_are_rows(refined_table)) {
      refined_data <- t(refined_data)
    }
    if (is.null(taxonomy)) {
      taxonomy <- collapse_taxonomy(refined_table)
    }
  }
  else {
    refined_data <- refined_table %>% as.matrix()
  }
  ID <- colnames(refined_data)
  meanAbundance <- base::colMeans(refined_data)
  medianAbundance <- matrixStats::colMedians(refined_data)
  maxAbundance <- matrixStats::colMaxs(refined_data)
  numberNonZero <- base::colSums(refined_data != 0)
  proportionNonZero <- base::colMeans(refined_data != 0)
  res <- data.frame(
    ID, meanAbundance, medianAbundance, maxAbundance,
    numberNonZero, proportionNonZero
  )
  if (!is.null(taxonomy)) {
    # Sometimes, we may risk that the taxonomy table is not the correct order
    taxonomy <- taxonomy[ID]
    res$taxonomy <- taxonomy
  }
  return(res)
}

#' @title countInteractions
#'
#' @description  Counts the number of interactions for each OTU
#'
#' @param IDs The name of each OTU
#'
#' @param interactions_table An \code{interaction_table}
#'
countInteractions <- function(IDs, interactions_table) {
  dataset_OTUs <- c(interactions_table$`OTU_1`, interactions_table$`OTU_2`)
  counts <- vapply(IDs, function(ID) {
    sum(dataset_OTUs == ID)
  }, FUN.VALUE = 1)
  names(counts) <- IDs
  return(counts)
}
# Computes the abunance product for at set of interactions
abundanceProduct <- function(interactions_table, OTU_stat, type = "mean") {
  abundances <- switch(type,
    "mean" = OTU_stat$meanAbundance,
    "median" = OTU_stat$medianAbundance,
    "max" = OTU_stat$maxAbundance
  )
  names(abundances) <- OTU_stat$ID
  abundance_1 <- abundances[interactions_table$OTU_1]
  abundance_2 <- abundances[interactions_table$OTU_2]
  abundance_1 * abundance_2
}

# @title overall_stats
# @description Compiles a summary of significant interactions for the given similarity measures
# @param similarity_measures_significance List of \code{interaction_table} objects returned from
# \link{create_interaction_table}
# @param outputargs The corresponding outputargs list sendt to \link{create_interaction_table}
# @param A \code{data.frame} consisting of the following fields: \itemize{
# \item \code{name} The name of the similarity measure
#
# \item \code{type} The type of similarity measure, equal to the one in the \code{type}
# slot in it \code{sim.measure} object
# \item \code{num_significant} The number of significant interaction
#
# }
# #overall_stats <- function(similarity_measures_significance) {#  overall <- data.frame()#}
