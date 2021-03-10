#' @title runAnalysis
#'
#' @description
#' Runs an automized processing of the OTU table or phyloseq object, passes the jobs to \code{ccrepe} and saves the results
#'
#' @param OTU_table The raw OTU table (if a \code{data.frame} or a \code{phyloseq} \code{otu_table} is supplied) to be treated
#' \emph{or} an experiment level \code{phyloseq} object containing the data (the latter is recommended).
#' Note that in the case of a \code{phyloseq} \code{otu_table}, no taxonomy can be handled.
#'
#' @param abundance_cutoff The mean abundance cutoff for the OTUs. If it is \code{NULL}, the there will be not filtering.
#'
#' @param parallel Should the analysis be run in parallel?
#'
#' @param sim.scores The similarity measures of class \link{sim.measure}
#' to use. If it is \code{NULL}, all measures available in the package
#' will be used (recommanded for most purposes).
#'
#' @param subset Character, the subset of similarity measures to use, denoted by
#' the its name in the list (not necessarly its string) returned from \link{similarity_measures} or similarity
#' measure modiftying function such as \link{noisify}
#' If \code{NULL}, all available measures will be used
#'
#' @param file Should the tables of significant interactions be written to a file? If so, they are printed to \code{csv}-files containing the name of the similarity measure
#' @param returnVariables Which variables should the function return (character vector)?
#' Available options are: \itemize{
#' \item \code{similarity_measures_significance}:
#' The \code{interaction_table} of significant interactions
#' \item \code{refined_table}: The processed  OTU table
#' \item \code{min_dataset}: The smallest non-zero entity in the refined table
#' \item \code{taxonomy}: A named numberic containing the taxonomy of each OTU (collapsed into a single string)
#' \item \code{outputargs}: A list (with a element for each similarity measure) comtaining the arguments to be passed
#' to \link{output_ccrepe_data} for each similarity measure
#' \item \code{common_outputargs}: Like \code{outputargs}, but these arguments stay the same for all similarity
#' measure in order to avoid duplicates.
#' }
#'
#' In addition, all paramerters for this function are available. Other internal variables
#' found upon inspection of the source code may also be returned,
#' but they are for advanced users only.
#' If \code{NULL}, the listed parameters in this section in addtion to an echo of the parameters are retuned.
#' Note: If \code{OTU_table} is a \code{phyloseq} object, the returned variables is a data frame corresponding to the
#' \code{phyloseq} object. This is due to the fact that it is intervally converted into a data frame.
#'
#' @param magnitude_factor When making noisified functions, the magnitude of the noise
#' will be this number multiplied with \code{min_dataset}
#'
#' @param q_crit Numeric, the q-value cutoff when construction interaction tables
#'
#' @param prefix The prefix of the file names being written. Ignored if \code{file=FALSE}.
#'
#' @param postfix The postfix of the file names being written. Ignored if \code{file=FALSE}.
#'
#' @param metadataCols The names (character vector) or position (integer) of the
#' metadata columns to remove from the table before analyzing it. Ignored if a \code{phyloseq} object is supplied
#'
#' @param renormalize Should the data be renormalized during filtering process and permutation?
#' Should be \code{TRUE} when used on relative abundances, but must be \code{FALSE} if absolute abundances are used.
#'
#' @param iterations Integer of length one, the number of iterations to run
#'
#' @return A list of the variables requested from the parameter \code{returnVariables}.
#'
#' @details
#' If the function is told to output a file and no prefix is given, the csv-files will all share a common prefix of the form:
#' \code{q_crit=(critical q-value)_cutoff=(the mean abundance cutoff)_magfac=(the magnitude factor)},
#' where all numbers are in scientific notation. Then the sim.score name follows, then the postfix and finally the csv
#' extention.
#' The postfix is by default empty.
#'
#' In order for an OTU-table to be valid when the argument \code{OTU_table} is a \code{data.frame},
#' the following criteria must hold:
#'
#' \itemize{
#' \item
#' The data points (sample) are in columns, the abundances for each
#' OTU is in rows.
#' \item
#' The rows may only hold OTU abundances
#' \item
#' There may be as many metadata colums as you like. However, they all
#' need to be declared in the \code{metadataCols} argument and the column
#' \code{taxonomy} has be there in order for the output file to contain the
#' taxonomy.
#' \item The row names of the table are the OTU names and the column names are the
#' sample names
#' }
#' For \code{phyloseq} objects (both experiment level and \code{otu_table}), you do not need to care about this, it is automatically handeled
#'
#' @seealso \link{output_ccrepe_data}
#'
#' @examples
#' library(micInt)
#' data(seawater)
#' sim.scores <- similarity_measures(subset= c("spearman","pearson"))
#' runAnalysis(OTU_table = seawater, sim.scores = sim.scores, parallel = FALSE,
#' iterations = 100)
#'
#' @import phyloseq
#' @export
runAnalysis <- function(OTU_table, abundance_cutoff = 1e-04, q_crit = 0.05, parallel = TRUE,
                        returnVariables = NULL, subset = NULL, sim.scores = NULL, file = FALSE, magnitude_factor = 10, prefix = NULL,
                        metadataCols = c("OTU Id", "taxonomy"),
                        postfix = "", renormalize=TRUE,iterations = 1000) {
  if (is.null(prefix)) {
    prefix <- create_prefix(q_crit = q_crit, cutoff = abundance_cutoff, magfac = magnitude_factor)
  }
  if (inherits(OTU_table, "phyloseq")) {
    phyloseq_object <- OTU_table
    taxonomy <- collapse_taxonomy(phyloseq_object)
    OTU_table <- phyloseq::otu_table(OTU_table) %>% data.frame()
    if (!phyloseq::taxa_are_rows(phyloseq_object)) {
      OTU_table <- t(OTU_table) %>% data.frame()
    }
    metadataCols <- NULL
  }
  else if (inherits(OTU_table,"otu_table")){
    phyloseq_otu_table <- OTU_table
    taxonomy <- NULL
    metadataCols <- NULL
    OTU_table <- OTU_table %>% data.frame()
    if (!phyloseq::taxa_are_rows(phyloseq_otu_table)) {
      OTU_table <- t(OTU_table) %>% data.frame()
    }
  }
  # In case we are provided with a data frame instead
  else {
    # If we are provided a data frame, we must see if the taxonomy column exists
    if ("taxonomy" %in% colnames(OTU_table)) {
      taxonomy <- OTU_table$taxonomy
      names(taxonomy) <- rownames(OTU_table)
    }
    else {
      taxonomy <- NULL
    }
  }
  refined_table <- refine_data(OTU_table, abundance_cutoff = abundance_cutoff,
                               metadataCols = metadataCols,renormalize=renormalize)
  # The smallest non-zero value in the data set
  min_dataset <- min(apply(refined_table, MARGIN = 2, function(x) min(x[x > 0])))
  magnitude <- magnitude_factor * min_dataset
  if (is.null(sim.scores)) {
    sim.scores <- noisify(magnitude = magnitude)
  }
  ccrepe_job <- create_ccrepe_jobs(sim.scores = sim.scores,
      prefix = prefix, postfix = paste0(postfix, ".csv")
    )
  ccrepe_commonargs <- list(x = refined_table, min.subj = 10,
                            verbose = FALSE,iterations = iterations,
                            renormalize=renormalize)
  if (!is.null(subset)) {
    sim.scores <- sim.scores[subset]
    ccrepe_job <- ccrepe_job[subset]
  }
  stringlist <- lapply(ccrepe_job, function(x) list(string = x$string))
  ccrepe_res <- ccrepe_analysis(ccrepe_job=lapply(ccrepe_job,function(x)
    x$ccrepe_args),
                                                  commonargs=ccrepe_commonargs, parallel = parallel)
  for (i in 1:length(ccrepe_job)) {
    ccrepe_res[[i]] <- c(list(res=ccrepe_res[[i]]), ccrepe_job[[i]]$output_args)
  }
  outputargs <- add_outputargs(ccrepe_res)
  common_outputargs=list(threshold.value=q_crit,return.value=TRUE,taxonomy=taxonomy,output.file=file,
                         threshold.type='q',removeDuplicates = TRUE,csv_option='1')
  similarity_measures_significance <- lapply(
    outputargs,
    function(x) (do.call(
        what = micInt::output_ccrepe_data,
        args = c(x,common_outputargs)
      ))
  )
  if (is.null(returnVariables)) {
    parameters <- formalArgs('runAnalysis')
    API_variables <- c('similarity_measures_significance','outputargs','common_outputargs',
                       'refined_table','min_dataset','taxonomy')
    mget(c(parameters,API_variables)) %>% return()
  }
  else {
    return(mget(returnVariables))
  }
}


# Diagnostic plots --------------------------------------------------------
#' @title autoplot.interaction_table
#'
#' @param object An \code{interaction_table}
#'
#' @param OTU_stat OTU statistics obtained from the function \link{OTU_stats}
#'
#' @param type One of the following: \itemize{
#' \item \code{'num_int'} Plots the number of interactions for each OTU against
#' its abundance
#' \item \code{'ab_prod'} Plots the product of the abundances of each significant pairs of
#' OTUs against the q- or p-value of the interaction
#' \item \code{'sim_cor'} Plots the absolute value of \code{sim.score} of each significant interactions against
#' the q- or p-value of the interaction
#' }
#'
#' @param cutoff_type
#' One of the following:
#' \itemize{
#' \item \code{'q'} Plots according to the q-values
#' \item \code{'p'} Plots according to the p-values
#' }
#' This parameter is ignored if \code{type='num_int'}
#'
#' @param abundance_type The type of abundance to use in the plots. One of \code{c('mean','median','max')}.
#'
#'
#' @param ... Currently ignored
#'
#'
#' @description Makes a diagnostic plot of an \code{interaction_table}
#'
#' @details Note that if \code{type='sim_cor'}, a column named \code{sign}, indicating the sign of the interaction,
#' is available even though the plotting function does not use this feature.
#'
#' @import ggplot2
#' @importFrom rlang .data
#'
#' @return A \link{ggplot} object showing the desired diagnostic plot
#'
#' @examples
#' library(micInt)
#' data("seawater")
#' sim.scores <- similarity_measures(subset= c("spearman","pearson"))
#' res <- runAnalysis(OTU_table = seawater,sim.scores = sim.scores,returnVariables = 'ccrepe_res',iterations = 100,parallel = FALSE)
#' int_table <- create_interaction_table(res$ccrepe_res$spearman$res)
#' stats <- OTU_stats(seawater)
#' autoplot(int_stable,stats,type = "ab_prod",cutoff_type = "p",abundance_type = "max")
#' @export
autoplot.interaction_table <- function(object, OTU_stat, type = "num_int", cutoff_type = "q",
                                       abundance_type = "mean", ...) {
  plot_after <- switch(abundance_type,
    "mean" = OTU_stat$meanAbundance,
    "max" = OTU_stat$maxAbundance,
    "median" = OTU_stat$medianAbundance
  )
  if (type == "num_int") {
    num_int <- countInteractions(OTU_stat$ID, object)
    xlab <- switch(abundance_type,
      "mean" = "Mean abundance",
      "max" = "Max abundance",
      "median" = "Median abundance"
    )
    return(ggplot() + geom_point(aes(x = plot_after, y = num_int)) + scale_x_continuous(name = xlab, trans = "log10") +
      scale_y_continuous(name = "Number of interactions"))
  }
  else if(type == "ab_prod" || type == "sim_cor"){
    if (cutoff_type == "q") {
      description <- "q-value"
      valueColumn <- "q.value"
    }
    else {
      description <- "p-value"
      valueColumn <- "p.value"
    }
    if(type == 'ab_prod'){
    abundance_product <- abundanceProduct(object, OTU_stat, type = abundance_type)
    y <- abundance_product
    x <- object[[valueColumn]]
    return(ggplot() + geom_point(aes(x = x, y = y)) + xlab(description) + ylab("Abundance product") +
      scale_x_continuous(trans = "log10") + scale_y_continuous(trans = "log10"))
    }
    else{
      plot_frame <- object %>% dplyr::transmute(.data$sim.score, !! rlang::sym(valueColumn), sign = sign(.data$sim.score))
      return(ggplot(plot_frame) + geom_point(aes(x=!! rlang::sym(valueColumn), y = abs(.data$sim.score))) + xlab(description) +
        ylab('Absolute value of sim.score')+scale_x_log10())
    }
  }

  else {
    stop('Enter a valid plot type')
  }
}

#' @title ratio_shared_interactions
#'
#' @description
#' Finds the Jaccard index of number of shared interactions between the similarity measures
#'
#' @param similarity_measures_significance A list objects of the class \code{interaction_table}
#'
#' @return A matrix showing the Jaccard indecies of the number of significant interactions.
#' Measures with no significant interactions are ignored.
#' @examples
#' library(micInt)
#' data("seawater")
#' sim.scores <- similarity_measures(subset= c("spearman","pearson"))
#' res <- runAnalysis(OTU_table = seawater,sim.scores = sim.scores,
#' iterations = 100,parallel = FALSE)
#' int_tables <- res$similarity_measures_significance
#' ratio_shared_interactions(int_tables)
#' @export
ratio_shared_interactions <- function(similarity_measures_significance) {
  # Ignores tables with zero columns
  to_include <- unlist(lapply(similarity_measures_significance, function(x) nrow(x) != 0))
  OTU_pairs <- lapply(similarity_measures_significance[to_include], extract_OTUs)
  n_sim_measures <- length(OTU_pairs)
  res <- matrix(data = NA, nrow = n_sim_measures, ncol = n_sim_measures)
  dimnames(res) <- list(names(OTU_pairs), names(OTU_pairs))
  for (i in 1:(n_sim_measures - 1)) {
    for (j in (i + 1):n_sim_measures) {
      res[i, j] <- count_matches(OTU_pairs[[i]], OTU_pairs[[j]]) / min(length(OTU_pairs[[i]]), length(OTU_pairs[[j]]))
      res[j, i] <- res[i, j]
    }
  }
  return(res)
}
# Ensures that each pair of OTUs are ordered by the lowest index first
extract_OTUs <- function(sim_result) {
  OTUs <- sim_result[c("OTU_1", "OTU_2")]
  OTUs_orded <- as.data.frame(apply(OTUs, 1, function(x) c(min(x), max(x))))
  OTUs_list <- as.list(OTUs_orded)
  return(OTUs_list)
}
# Count the number of identical (not necessarly contigious) rows
count_matches <- function(OTUs_1, OTUs_2) {
  length(intersect(OTUs_1, OTUs_2))
}
