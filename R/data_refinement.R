#' @title Remove metadata from data set
#' @description
#' Removes desired metadata columns from dataset
#' @param OTU_table The raw OTU table
#' @param metadataCols The names (character vector) or position (integer) of the
#' metadata columns to remove from the table
remove_metadata <- function(OTU_table, metadataCols = c("OTU Id", "taxonomy")) {
  if (is.character(metadataCols)) {
    metadata <- which(colnames(OTU_table) %in% metadataCols)
  }
  else {
    metadata <- metadataCols
  }
  # Removes the metadata for the data set if there is any
  if (length(metadata) != 0) {
    OTU_table <- OTU_table[, -metadata]
  }
  refined_table <- as.data.frame(t(OTU_table))
  colnames(refined_table) <- row.names(OTU_table)
  return(refined_table)
}

#' @title cut_abundances
#'
#' @description
#' Removes metadata columns from dataset
#'
#' @param abundance_cutoff
#' Numeric, the threshold cutoff value. If it is \code{NULL}, the filtering process is skipped.
#'
#' @param type
#' The type of measure to base the cutoff on. Can be any
#' of \code{'mean'}, \code{median}, \code{max} which cuts away
#' OTUs based on mean, median and maximum abundance, repectivly
#'
#' @param renormalize
#' Logical, should the abundances be renormalized after the procedyre?
#'
#'
#' @import matrixStats
cut_abundances <- function(refined_table, abundance_cutoff = 0, type = "mean", renormalize = TRUE) {
  if (is.null(abundance_cutoff)) {
    # If we do not set an abundance cutoff, we skip the process entierly
    return(refined_table)
  }
  m_refined_table <- as.matrix(refined_table)
  abundances <- switch(type,
    mean = colMeans(refined_table),
    median = colMedians(m_refined_table),
    max = colMaxs(m_refined_table),
    numberNonZero = colSums(m_refined_table != 0),
    proprotionNonZero = colMeans(m_refined_table != 0)
  )
  refined_table <- refined_table[, abundances > abundance_cutoff]
  if (renormalize) {
    # Renormalizes the table in order for the rows to sum to 1
    refined_table <- renormalize(refined_table)
  }
  return(refined_table)
}
#' @title Refine raw OTU table
#' @description
#' Removes metadata from OTU table and cuts off the least abundant
#' species, defined by the cutoff parameter
#'
#' @param OTU_table The raw OTU table, either as a \code{data.frame} or a \code{matrix} (not phyloseq)
#' @inheritParams cut_abundances
#' @param metadataCols The names (character vector) or position (integer) of the
#' metadata columns to remove from the table
#'
#' @return A data frame with the metadata columns removed, and the OTUs
#' below the cutoff are filtered away. Additionally, in the refined table
#' the OTUs constitute the rows, while the rows are the samples (transposed
#' compared to the original OTU table).
#'
#' @details
#' In order for an OTU-table to be valid, the following criteria must hold:
#'
#' \itemize{
#' \item
#' The data points (sample) are in columns, the abundances for each
#' OTU is in rows.
#' \item
#' The rows may only hold OTU abundances
#' \item
#' There may be as many metadata colums as preferable. However, the all
#' need to be declared in the \code{metadataCols} argument and the column
#' \code{taxonomy} has be there in order for the output file to contain the
#' taxonomy.
#' \item The row names of the table are the OTU names and the column names are the
#' sample names
#' }
#'
#' @export
refine_data <- function(OTU_table, abundance_cutoff = 0, cutoff_type = "mean", renormalize = TRUE, metadataCols = c("OTU Id", "taxonomy")) {
  refined_table <- remove_metadata(OTU_table, metadataCols = metadataCols)
  # Cuts away the least abundant species
  cut_abundances(refined_table, abundance_cutoff, type = cutoff_type, renormalize = renormalize)
}
renormalize <- function(table) {
  matrix <- as.matrix(table)
  res <- diag(1 / rowSums(matrix)) %*% matrix
  table[, ] <- res
  return(table)
}

#' @title output_ccrepe_data
#' @description
#' Takes input from ccrepe and transforms it into a convenient table
#'
#' @param data The results from ccrepe
#'
#' @param taxonomy Named character, with the names being those of the OTUs and the values
#' their taxonomy collapased into a single string (for each OTU)
#'
#' @param threshold.type The type of threshold to be used \code{'q'} denotes q-value,
#'  while \code{'p'} denotes p-value
#'
#' @param threshold.value The critical significance value for including an interaction
#'
#' @param output.file Should the function write the date to a csv-file?
#'
#' @param filename What is the name of the output file (including path). Ignored if \code{output.file} is \code{TRUE}
#'
#' @param return.value Should the function return the interaction table as a return value?
#'
#' @param removeDuplicates If \code{TRUE}, the function ensures that no OTU pair is listed twice in the results.
#' Else, two interacting OTUs may show up twice in reverse order
#' (a pair where OTU A is listed as \code{OTU_1} and OTU B is listed was \code{OTU_2} and vice versa for the other pair)
#'
#' @param csv_option If assigned the value \code{'2'}, the function will use the \code{write.csv2} function which print
#' the numbers if comma as decimal delimer and semicolon as the delimer between numbers.
#' Else, the decimal delimer is point and comma the delimer between numbers.
#'
#' @param score_attributes An object of class \link{sim.measure.attributes} belonging to the similarity measure being used
#'
#'
#'
#' @importFrom utils modifyList write.csv write.csv2
#'
#' @export
output_ccrepe_data <- function(data, taxonomy = NULL, threshold.type = "q", threshold.value = 0.05, output.file = FALSE, filename = NULL,
                               return.value = TRUE, csv_option = "2", removeDuplicates = TRUE, sim.measure.attributes = NULL) {
  significant_interactions <- create_interaction_table(
    data = data, taxonomy = taxonomy, threshold.type = threshold.type,
    threshold.value = threshold.value,
    removeDuplicates = removeDuplicates, score_attributes = sim.measure.attributes
  )
  if (output.file) {
    write.interaction_table(significant_interactions,
      filename = filename,
      csv_option = csv_option
    )
  }
  if (return.value) {
    return(significant_interactions)
  }
  else {
    return(NULL)
  }
}
#' @title Create table of significant interactions
#'
#' @description Performs the actual process of making the table.
#' @return An \code{interaction_table}, being a \code{data.frame} with the following columns
#' (given that they are available in the input):
#' \itemize{
#' \item \code{OTU_1}, \code{OTU_2} The IDs of the two interacting OTUs
#' \item \code{z.stat} The \eqn{z}{z}-value of the interaction
#' \item \code{p.value} The p-value of the interaction (based only on the interaction itself)
#' \item \code{q.value} The q-value of the interaction (corrected for multiple testing)
#' \item \code{taxonomy_1}, \code{taxonomy_2} The taxonomy of the two interacting OTUs, collapsed into a single string.
#' }
#' In addition, the information from the given similarity scores are also provided
#'
#' @seealso \link{collapse_taxonomy} \link{sim.measure.attributes}
#'
#'
#'
#' @export
#'
#' @inheritParams output_ccrepe_data
create_interaction_table <- function(data, taxonomy = NULL, threshold.type = "q", threshold.value = 0.05,
                                     removeDuplicates = TRUE, score_attributes = NULL) {
  p.values <- as.data.frame(data$p.values)
  z.stat <- as.data.frame(data$z.stat)
  sim.score <- as.data.frame(data$sim.score)
  q.values <- as.data.frame(data$q.values)
  if (threshold.type == "q") {
    threshold_matrix <- q.values
  }
  else if (threshold.type == "p") {
    threshold_matrix <- p.values
  }
  else {
    stop("no valid theshold method given, must be 'p' (local p-value)
         or 'q' (familywise false discovery rate)")
  }
  significant_pairs <- which(threshold_matrix < threshold.value, arr.ind = TRUE)
  significant_interactions <- as.data.frame(matrix(colnames(threshold_matrix)[significant_pairs], ncol = 2))
  colnames(significant_interactions) <- c("OTU_1", "OTU_2")
  if (nrow(significant_interactions) == 0) {
    significant_interactions$sim.score <- numeric()
    significant_interactions$p.value <- numeric()
    significant_interactions$q.value <- numeric()
    significant_interactions$z.stat <- numeric()
    if (!is.null(taxonomy)) {
      # Add taxonomy if available
      significant_interactions$taxonomy_1 <- character()
      significant_interactions$taxonomy_2 <- character()
    }
  }
  else {
    # Removing duplicates
    if (removeDuplicates) {
      significant_interactions <- significant_interactions[significant_interactions$OTU_1 > significant_interactions$OTU_2, ]
    }
    if (!is.null(sim.score) && nrow(sim.score) != 0) {
      significant_interactions$sim.score <- apply(significant_interactions, 1, function(x) {
        sim.score[x[1], x[2]]
      })
    }
    if (!is.null(p.values) && nrow(p.values) != 0) {
      significant_interactions$p.value <- apply(significant_interactions, 1, function(x) {
        p.values[x[1], x[2]]
      })
    }
    if (!is.null(q.values) && nrow(q.values) != 0) {
      significant_interactions$q.value <- apply(significant_interactions, 1, function(x) {
        q.values[x[1], x[2]]
      })
    }
    if (!is.null(z.stat) && nrow(z.stat) != 0) {
      significant_interactions$z.stat <- apply(significant_interactions, 1, function(x) {
        z.stat[x[1], x[2]]
      })
    }
    if (!is.null(taxonomy)) {
      # Add taxonomy if available
      significant_interactions$taxonomy_1 <- taxonomy[significant_interactions$OTU_1]
      significant_interactions$taxonomy_2 <- taxonomy[significant_interactions$OTU_2]
    }
    # Sorts by dersired significance column
    if (threshold.type == "q") {
      significant_interactions <- significant_interactions[order(significant_interactions$q.value), ]
    } else {
      significant_interactions <- significant_interactions[order(significant_interactions$p.value), ]
    }
  }
  class(significant_interactions) <- c("interaction_table", class(significant_interactions))
  if (is.null(score_attributes)) {
    attr(significant_interactions, "measure_name") <- NA
    attr(significant_interactions, "signed") <- NA
    attr(significant_interactions, "measure_type") <- NA
  }
  else {
    attr(significant_interactions, "measure_name") <- score_attributes@string
    attr(significant_interactions, "signed") <- signed <- score_attributes@signed
    attr(significant_interactions, "measure_type") <- score_attributes@type_measure
  }
  return(significant_interactions)
}

#'
#' @title write.interaction_table
#'
#' @description Takes an \code{interaction_table} and writes it to a file. For information about parameters,
#' see \link{output_ccrepe_data}
#'
#' @export
write.interaction_table <- function(significant_interactions, filename,
                                     csv_option = "1") {
  if (csv_option == "2") {
    write.csv2(significant_interactions, file = filename, row.names = FALSE)
  }
  else {
    write.csv(significant_interactions, file = filename, row.names = FALSE)
  }
}
#' @name summary.interaction_table
#'
#' @title Create summary of \code{interaction_table}
#'
#' @param table An interaction table returned from \link{create_interaction_table}
#'
#' @return A list with the following fields (the first three are taken directly from
#' the \code{interaction.table}):
#' \itemize{
#' \item \code{measure_name}
#' \item \code{signed}
#' \item \code{measure_type}
#' \item \code{number_significant}: The number of significant interactions
#' \item \code{proportion_negative}: The proprotion of significant interactions which
#' are negative
#' }
#'
#'
#' @export
summary.interaction_table <- function(table) {
  proportion_negative <- ifelse(attr(table,'signed'),
                                sum(table$sim.score < 0) / nrow(table), NA)
  number_significant <- nrow(table)
  c(
    attributes(table)[c("measure_name", "signed", "measure_type")],
    list(number_significant = number_significant, proportion_negative = proportion_negative)
  )
}

#' @export
as.edgelist <- function(x, ...) {
  UseMethod("as.edgelist")
}

#' @name as.edgelist.interaction_table
#' @title Create edge list from \code{interaction_table}
#'
#' @param table An \code{interaction_table}
#'
#' @description Converts an interaction table to a two column character matrix
#' where each row represent an edge between the OTUs. The entities in the each row are the name of the
#' two interacting OTUs.
#'
#' @return The result of this function can be fed directly into \link[igraph]{graph_from_edgelist}
#' for making interacting network
#'
#' @export
as.edgelist.interaction_table <- function(table) {
  as.matrix(table[, c("OTU_1", "OTU_2")])
}

#' @title Collapse taxonomy into single strings
#'
#' @description Given a phyloseq object with a taxonomy table, this function extracts the taxonomy of each OTU and
#' collapse the result into a single string
#'
#' @param phyloseq A \code{\link{phyloseq}} object containing the \code{tax_table} slot. If this slot does not exist,
#' the function will return \code{NULL}.
#'
#' @return A named character vector, where the names are the OTU names and the values are the collapsed taxonomies.

collapse_taxonomy <- function(phyloseq) {
  taxonomy_table <- tax_table(phyloseq, errorIfNULL = FALSE)
  if (is.null(taxonomy_table)) {
    return(taxonomy_table)
  }
  taxonomy_frame <- data.frame(taxonomy_table)
  # Makes a contatented string of the taxonomies of each OTU
  taxonomy <- apply(taxonomy_frame, MARGIN = 1, FUN = function(x) paste(x, collapse = ","))
  names(taxonomy) <- rownames(taxonomy_frame)
  return(taxonomy)
}
