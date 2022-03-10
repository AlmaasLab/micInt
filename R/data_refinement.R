# @title Remove metadata from data set
# @description
# Removes desired metadata columns from dataset
# @param OTU_table The raw OTU table
# @param metadataCols The names (character vector) or position (integer) of the
# metadata columns to remove from the table
#
remove_metadata <-
  function(OTU_table,
           metadataCols = c("OTU Id", "taxonomy")) {
    if (is.character(metadataCols)) {
      metadata <- which(colnames(OTU_table) %in% metadataCols)
    }
    else {
      metadata <- metadataCols
    }
    # Removes the metadata for the data set if there is any
    if (length(metadata) != 0) {
      OTU_table <- OTU_table[,-metadata]
    }
    refined_table <- as.data.frame(t(OTU_table))
    colnames(refined_table) <- row.names(OTU_table)
    return(refined_table)
  }

# @title cut_abundances
#
# @description
# Removes metadata columns from dataset
#
# @param refined_table A \code{data.frame} where OTUs are in columns and samples in rows with metadata removed
#
# @param abundance_cutoff
# Numeric, the threshold cutoff value. If it is \code{NULL}, the filtering process is skipped.
#
# @param cutoff_type
# The type of measure to base the cutoff on. Can be any
# of \code{'mean'}, \code{median}, \code{max} which cuts away
# OTUs based on mean, median and maximum abundance, repectivly
#
# @param renormalize
# Logical, should the abundances be renormalized after the procedyre?
#
#
cut_abundances <-
  function(refined_table,
           abundance_cutoff = 0,
           cutoff_type = "mean",
           raw_value_cutoff = TRUE,
           renormalize = FALSE) {
    if (is.null(abundance_cutoff)) {
      # If we do not set an abundance cutoff, we skip the process entierly
      return(refined_table)
    }
    m_refined_table <- as.matrix(refined_table)
    if(!raw_value_cutoff){
      m_refined_table <- renormalize(m_refined_table)
    }
    abundances <- switch(
      cutoff_type,
      mean = base::colMeans(refined_table),
      median = matrixStats::colMedians(m_refined_table),
      max = matrixStats::colMaxs(m_refined_table),
      numberNonZero = base::colSums(m_refined_table != 0),
      proprotionNonZero = base::colMeans(m_refined_table != 0)
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
#' @param OTU_table The raw OTU table, either as a \code{data.frame}, a \code{matrix} or a \code{phyloseq} object
#'
#' @param abundance_cutoff
#' Numeric, the threshold cutoff value. If it is \code{NULL}, the filtering process is skipped.
#'
#' @param cutoff_type
#' The type of measure to base the cutoff on. Can be any
#' of \code{'mean'}, \code{median}, \code{max} which cuts away
#' OTUs based on mean, median and maximum abundance, repectivly
#'
#' @param raw_value_cutoff
#'
#' Logical, should filtering be based on the raw abundances? If not,
#' the sample-wise relative abundances are used for filtering. Note that this
#' parameter does not determine whether the \emph{results} of the function are
#' relative abundances.
#'
#' @param renormalize
#' Logical, should the abundances be renormalized (sample-wise) after the procedure?
#'
#' @param metadataCols The names (character vector) or position (integer) of the
#' metadata columns to remove from the table
#'
#' @return A data frame (always a base data frame even though a tibble is supplied)
#' with the metadata columns removed, and the OTUs
#' below the cutoff are filtered away. Additionally, in the refined table
#' the OTUs constitute the rows, while the rows are the samples (transposed
#' compared to the original OTU table).
#'
#' @details
#' \subsection{Critera for OTU tables}{
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
#' Of course, this does not apply to the case when a \code{phyloseq} object is provided as everything
#'  will be handled automatically
#' }
#' \subsection{Tibbles}{
#' Tibbles are troublesome in this context as the do not support rownames. In order to obtain the OTU IDs, the function
#' looks through the following in order and proceed with the next until successfull:
#' \enumerate{
#' \item If the original OTU table has rownames, they are used as the OTU IDs
#' \item If one (the first) of \code{c('OTU Id','#OTU ID','OTU ID','OTU_ID')} reside in the colnames of the OTU table,
#' the corresponding column is used as the OTU IDs
#' \item If \code{metadataCols} is not empty, the first of the specified columns regarded as the OTU ID (gives warning)
#' \item The numeric row indicies are treated as the OTU IDs (gives warning)
#' }
#' }
#'
#' @examples
#' library(micInt)
#' data("seawater")
#' refine_data(seawater,abundance_cutoff = 1e-3,cutoff_type = "max")

#'

#'
#' @export
refine_data <-
  function(OTU_table,
           abundance_cutoff = 0,
           cutoff_type = "mean",
           raw_value_cutoff = TRUE,
           renormalize = FALSE,
           metadataCols = c("OTU Id", "taxonomy")) {
    if (inherits(OTU_table, 'phyloseq') || inherits(OTU_table,'otu_table'))
      {
      refined_table <- phyloseq::otu_table(OTU_table) %>% data.frame(check.names = FALSE)
      if (phyloseq::taxa_are_rows(OTU_table)) {
        refined_table %<>% t()
      }
    }
    else{
      refined_table <-
        remove_metadata(OTU_table, metadataCols = metadataCols)
      if(inherits(OTU_table,'tbl_df')){
        # If we are provided a tibble as the OTU table, we should proceed with care,
        # the refined_table
        # MUST be a base data frame (tidy data just makes a mess)
        # First of all, we need to find the column in the OTU table containing the OTU Id
        if(tibble::has_rownames(OTU_table)){
        # Our first guess is that they already reside inside the OTU table (but are removed upon subsetting)
          row_names <- rownames(OTU_table)
        }
        else{
          if(length(metadataCols) == 0){
            warning('No rownames of the OTU table are found, nor does the data have any metadata,
                    using the numeric indicies of the rows as OTU IDs')
            row_names <- seq_len(nrow(OTU_table))
          }
          else{
          candidate_names <- c('OTU Id','#OTU ID','OTU ID','OTU_ID')
          matched_candidates <- candidate_names[candidate_names %in% names(OTU_table)]
          if(length(matched_candidates) == 0){
            # Of cource, the column could have another name, so we assume it to be the first of the metadata columns
            warning('OTU ID column could not be automatically determined,
                    picking the first metadata column of the OTU table')

          ID_col <- metadataCols[1]
          }
          else{
            ID_col <- matched_candidates[1]
          }
          row_names <- OTU_table[[ID_col]]
          }
        }
        refined_table <- data.frame(refined_table,row.names = names(OTU_table),check.names = FALSE,
                                    stringsAsFactors = FALSE)
        names(refined_table) <- row_names
      }
    }
  # Cuts away the least abundant species
  cut_abundances(refined_table,
                abundance_cutoff,
                raw_value_cutoff = raw_value_cutoff,
                cutoff_type = cutoff_type,
                renormalize = renormalize)
  }
renormalize <- function(table) {
  matrix <- as.matrix(table)
  res <- diag(1 / rowSums(matrix)) %*% matrix
  table[,] <- res
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
#' @param score_attributes An object of class \linkS4class{sim.measure.attributes} belonging to the similarity measure being used
#'
#' @examples
#' library(micInt)
#' sim.scores <- similarity_measures(subset= c("spearman","pearson"))
#' res <- runAnalysis(OTU_table = seawater,
#' sim.scores = sim.scores,returnVariable = 'ccrepe_res',
#' iterations = 100,parallel = FALSE)
#' output_ccrepe_data(res$ccrepe_res$spearman$res)
#'
#'
#' @importFrom utils modifyList write.csv write.csv2
#'
#' @export
output_ccrepe_data <-
  function(data,
           taxonomy = NULL,
           threshold.type = "q",
           threshold.value = 0.05,
           output.file = FALSE,
           filename = NULL,
           return.value = TRUE,
           csv_option = "2",
           removeDuplicates = TRUE,
           score_attributes = NULL) {
    significant_interactions <- create_interaction_table(
      data = data,
      taxonomy = taxonomy,
      threshold.type = threshold.type,
      threshold.value = threshold.value,
      removeDuplicates = removeDuplicates,
      score_attributes = score_attributes
    )
    if (output.file) {
      write.interaction_table(significant_interactions,
                              filename = filename,
                              csv_option = csv_option)
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
#' \item \code{sim.score} The similarity measure score for the interaction
#' \item \code{z.stat} The \eqn{z}{z}-value of the interaction
#' \item \code{p.value} The p-value of the interaction (based only on the interaction itself)
#' \item \code{q.value} The q-value of the interaction (corrected for multiple testing)
#' \item \code{taxonomy_1}, \code{taxonomy_2} The taxonomy of the two interacting OTUs, collapsed into a single string.
#' }
#' In addition, the information from the given similarity scores are also provided
#'
#' @seealso \link{collapse_taxonomy} \linkS4class{sim.measure.attributes}
#'
#' @examples
#' sim.scores <- similarity_measures(subset= c("spearman","pearson"))
#' res <- runAnalysis(OTU_table = seawater,sim.scores = sim.scores,
#' returnVariables = 'ccrepe_res',iterations = 100,parallel = FALSE)
#' collapsed_taxonomy <- collapse_taxonomy(seawater)
#' create_interaction_table(res$ccrepe_res$spearman$res,taxonomy = collapsed_taxonomy,
#' threshold.value = 0.05, threshold.type = "q",
#'  score_attributes = sim.measure.attributes(sim.scores[["spearman]]))
#'
#' @export
#'
#' @inheritParams output_ccrepe_data
create_interaction_table <-
  function(data,
           taxonomy = NULL,
           threshold.type = "q",
           threshold.value = 0.05,
           removeDuplicates = TRUE,
           score_attributes = NULL) {
    p.values <- data$p.values
    z.stat <- data$z.stat
    sim.score <- data$sim.score
    q.values <- data$q.values
    if (threshold.type == "q") {
      threshold_matrix <- q.values
    }
    else if (threshold.type == "p") {
      threshold_matrix <- p.values
    }
    else {
      stop(
        "no valid theshold method given, must be 'p' (local p-value)
        or 'q' (familywise false discovery rate)"
      )
    }
    significant_pairs <-
      which(threshold_matrix < threshold.value, arr.ind = TRUE)
    # Removing duplicates
    if (removeDuplicates) {
      significant_pairs <-
        significant_pairs[significant_pairs[,1] > significant_pairs[,2],]
    }
    significant_matrix <- matrix(colnames(threshold_matrix)[significant_pairs],
           ncol = 2)
    significant_interactions <-
      as.data.frame(significant_matrix, stringsAsFactors=FALSE)
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
      subscript_matrix <- significant_pairs
      if (!is.null(sim.score) && nrow(sim.score) != 0) {
        significant_interactions$sim.score <-
          sim.score[subscript_matrix]
      }
      if (!is.null(p.values) && nrow(p.values) != 0) {
        significant_interactions$p.value <-
          p.values[subscript_matrix]
      }
      if (!is.null(q.values) && nrow(q.values) != 0) {
        significant_interactions$q.value <-
          q.values[subscript_matrix]
      }
      if (!is.null(z.stat) && nrow(z.stat) != 0) {
        significant_interactions$z.stat <-
          z.stat[subscript_matrix]
      }
      if (!is.null(taxonomy)) {
        # Add taxonomy if available
        significant_interactions$taxonomy_1 <-
          taxonomy[significant_interactions$OTU_1]
        significant_interactions$taxonomy_2 <-
          taxonomy[significant_interactions$OTU_2]
      }
      # Sorts by dersired significance column
      if (threshold.type == "q") {
        significant_interactions <-
          significant_interactions[order(significant_interactions$q.value),]
      } else {
        significant_interactions <-
          significant_interactions[order(significant_interactions$p.value),]
      }
    }
    class(significant_interactions) <-
      c("interaction_table", class(significant_interactions))
    if (is.null(score_attributes)) {
      attr(significant_interactions, "measure_name") <- NA
      attr(significant_interactions, "signed") <- NA
      attr(significant_interactions, "measure_type") <- NA
    }
    else {
      attr(significant_interactions, "measure_name") <-
        score_attributes@string
      attr(significant_interactions, "signed") <-
         score_attributes@signed
      attr(significant_interactions, "measure_type") <-
        score_attributes@type_measure
    }
    return(significant_interactions)
    }

#'
#' @title write.interaction_table
#'
#' @description Takes an \code{interaction_table} and writes it to a file. For information about parameters,
#' see \link{output_ccrepe_data}
#'
#' @param significant_interactions An \code{interaction_table} returned from \code{\link{create_interaction_table}}
#'
#' @inheritParams output_ccrepe_data
#' @examples
#' library(micInt)
#' data("seawater")
#' sim.scores <- similarity_measures(subset= c("spearman","pearson"))
#' res <- runAnalysis(OTU_table = seawater,iterations = 100,
#' sim.scores = sim.scores,parallel = FALSE)
#' int_table <- res$similarity_measures_significance$spearman
#' write.interaction_table(int_table,"ccrepe_results.csv",csv_option = "1")
#' # Ensures the file is deleted after running example
#' unlink("ccrepe_results.csv")

#' @export
write.interaction_table <-
  function(significant_interactions,
           filename,
           csv_option = "2") {
    if (csv_option == "2") {
      write.csv2(significant_interactions,
                 file = filename,
                 row.names = FALSE)
    }
    else {
      write.csv(significant_interactions,
                file = filename,
                row.names = FALSE)
    }
  }
#' @name summary.interaction_table
#'
#' @title Create summary of \code{interaction_table}
#'
#' @param object An interaction table returned from \link{create_interaction_table}
#'
#' @param ... further arguments passed to or from other methods
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
#' library(micInt)
#' sim.scores <- similarity_measures(subset= c("spearman","pearson"))
#' res <- runAnalysis(OTU_table = seawater,sim.scores = sim.scores,
#' iterations = 100, parallel = FALSE)
#' int_table <- res$similarity_measures_significance$spearman
#' summary(int_table)
#'
#'
#' @export
summary.interaction_table <- function(object, ...) {
  signed_attribute <- attr(object, 'signed')
  proportion_negative <- ifelse(is.null(signed_attribute) || signed_attribute,
                                sum(object$sim.score < 0) / nrow(object),
                                NA)
  number_significant <- nrow(object)
  res <- c(
    attributes(object)[c("measure_name", "signed", "measure_type")],
    list(
      number_significant = number_significant,
      proportion_negative = proportion_negative
    )
  )
  lapply(res, function(attr) if(is.null(attr)){
    NA
  }
  else{
    attr
  })
}

#' @name edgelist
#' @title Convert an object to a list of edges
#'
#' @description This is a S3 generic to convert a suitable object represetning (on some way) a graph to an
#' edgelist which can be inported to \code{igraph} by the function \link{graph_from_edgelist}.
#'
#' @param x The R object to convert to an edgelist
#'
#' @param ... further arguments passed to or from other methods.
#'
#' @examples
#' library(micInt)
#' library(igraph)
#' sim.scores <- similarity_measures(subset= c("spearman","pearson"))
#' res <- runAnalysis(OTU_table = seawater,iterations = 100,sim.scores = sim.scores,parallel = FALSE)
#' int_table <- res$similarity_measures_significance$spearman
#' edgelist <- as.edgelist(int_table)
#' graph_from_edgelist(edgelist)
#' @export
as.edgelist <- function(x, ...) {
  UseMethod("as.edgelist")
}

#' @name as.edgelist.interaction_table
#' @title Creates an \code{igraph} object from \code{interaction_table}
#'
#' @param x An \code{interaction_table}
#'
#' @param ... Ignored
#'
#' @description Converts an interaction table to a two column character matrix
#' where each row represent an edge between the OTUs. The entities in the each row are the name of the
#' two interacting OTUs.
#'
#' @return The result of this function can be fed directly into \link[igraph]{graph_from_edgelist}
#' for making interacting network
#'
#' @examples
#' library(micInt)
#' library(igraph)
#' sim.scores <- similarity_measures(subset= c("spearman","pearson"))
#' res <- runAnalysis(OTU_table = seawater,sim.scores = sim.scores, iterations = 100,parallel = FALSE)
#' int_table <- res$similarity_measures_significance$spearman
#' edgelist <- as.edgelist(int_table)
#' graph_from_data_frame(edgelist)
#' @export
as.edgelist.interaction_table <- function(x, ...) {
  as.matrix(x[, c("OTU_1", "OTU_2")])
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
#' @examples
#' library(micInt)
#' collapsed <- collapse_taxonomy(seawater)
#' collapsed
#' @export

collapse_taxonomy <- function(phyloseq) {
  taxonomy_table <- tax_table(phyloseq, errorIfNULL = FALSE)
  if (is.null(taxonomy_table)) {
    return(taxonomy_table)
  }
  taxonomy_frame <- data.frame(taxonomy_table)
  # Makes a contatented string of the taxonomies of each OTU
  taxonomy <-
    apply(
      taxonomy_frame,
      MARGIN = 1,
      FUN = function(x)
        paste(x, collapse = ",")
    )
  names(taxonomy) <- rownames(taxonomy_frame)
  return(taxonomy)
}


#' @title Scale a phyloseq object
#'
#' @description This function scales the abundances in each sample by a constant, making it suitable to transform from
#' relative data to absolute data.
#'
#' @param object A \code{phyloseq} object, where the abundances in its OTU table are to be scaled
#' @param column Either a character vector of length 1 with the name of the column in the sample data to scale by or a
#' numeric vector which length is equal to the number of samples.
#'
#' @return A \code{phyloseq} object where the abundances in the OTU table are scaled
#' @examples
#' library(micInt)
#' library(phyloseq)
#' data("soilrep")
#' s <- 1:nsamples(soilrep)
#' scale_by_column(soilrep, s)
#' @export
scale_by_column <- function(object, column) {
  if (is.character(column)) {
    if (length(column) != 1) {
      stop("The input must be a character vector of length 1")
    }
    if (column %in% phyloseq::sample_variables(object) %>% magrittr::not()) {
      stop("The selected column does not exist in the sample data")
    }
    scale_vector <- sample_data(object)[[column]]
  }
  else{
    if (length(column) != nsamples(object)) {
      stop("The length of the input must equal the number of samples")
    }
    scale_vector <- column
  }
  if (phyloseq::taxa_are_rows(object)) {
    OTU_table <-
      scale(
        phyloseq::otu_table(object),
        center = FALSE,
        scale = 1 / scale_vector
      ) %>%
      phyloseq::otu_table(taxa_are_rows = TRUE)
  }
  else{
    OTU_table <-
      t(scale(
        t(phyloseq::otu_table(object)),
        center = FALSE,
        scale = 1 / scale_vector
      )) %>%
      phyloseq::otu_table(taxa_are_rows = FALSE)
  }
  otu_table(object) <- OTU_table
  return(object)
}
