#' @title Get premade similarity measure objects
#'
#' @description
#' This function allows retriving the ten ready-made similarity measure
#'provided by this package. Custom similarity measures can also be made \link{sim.measure-class}.
#'
#' @param subset The subset of the similarity measures to be returned.
#' By default, all similarity measures are returned
#'
#' @return A list of the \link{sim.measure} objects defined by the function
#'
#' @import vegan
#' @import infotheo
#' @importFrom stats cor
#' @export
similarity_measures <- function(subset = c("pearson", "spearman", "kendall",
                                           "bray_curtis", "jaccard",  "gen_jaccard",
                                           "mutual_information", "nc.score", "euclidean",
                                           "cosine")) {
  measures <- list()
  #************************************************************
  # Pearson linear correlation
  #***********************************************************
  pearson_cor <- function(x, y = NULL) {
    cor(x = x, y = y, method = "pearson")
  }
  measures$pearson <- sim.measure(
    FUN = pearson_cor, string = "pearson", signed = TRUE,
    type = "parametric"
  )
  #************************************************************
  # Spearman correlation
  #***********************************************************
  spearman_cor <- function(x, y = NULL) {
    cor(x = x, y = y, method = "spearman")
  }
  measures$spearman <- sim.measure(
    FUN = spearman_cor, string = "spearman", mean_scaleable = FALSE,
    signed = TRUE, type = "non-parametric"
  )
  #************************************************************
  # Kendall's tau
  #***********************************************************
  kendall_cor <- function(x, y = NULL) {
    cor(x = x, y = y, method = "kendall")
  }
  measures$kendall <- sim.measure(
    FUN = kendall_cor, string = "kendall", mean_scaleable = FALSE,
    signed = TRUE, type = "non-parametric"
  )
  #************************************************************
  # Bray-Curtis
  #*************************************************************
  # We must transform the disimilarities into similarities
  bray_curtis_cor <- function(x, y = NULL) {
    if (is.null(y)) {
      ones <- rep(1, dim(x)[2])
      res <- ones %*% t(ones) - as.matrix(vegdist(t(x), method = "bray", upper = TRUE))
    }
    else {
      # The vegdist-function accepts matrix inputs only, but this can be fixed
      # the following way:
      res <- bray_curtis_cor(cbind(x, y))[1, 2]
    }
    return(res)
  }
  measures$bray_curtis <- sim.measure(
    FUN = bray_curtis_cor, string = "bray_curtis", mean_scaleable = TRUE,
    signed = FALSE, type = "parametric"
  )
  #************************************************************
  # Jaccard index
  #*************************************************************
  jaccard_cor <- function(x, y = NULL) {
    if (is.null(y)) {
      ones <- rep(1, dim(x)[2])
      res <- ones %*% t(ones) - as.matrix(vegdist(t(x), method = "jaccard", binary = TRUE, upper = TRUE))
    }
    else {
      res <- jaccard_cor(cbind(x, y))[1, 2]
    }
    return(res)
  }
  measures$jaccard <- sim.measure(
    FUN = jaccard_cor, string = "jaccard_index", categorical = TRUE,
    signed = FALSE, type = "presence-absence"
  )
  #************************************************************
  # Generalized Jaccard
  #*************************************************************
  gen_jaccard_cor <- function(x, y = NULL) {
    if (is.null(y)) {
      ones <- rep(1, dim(x)[2])
      res <- ones %*% t(ones) - as.matrix(designdist(t(x), method = "1-J/(A+B-J)", terms = "minimum"))
    }
    else {
      res <- gen_jaccard_cor(cbind(x, y))[1, 2]
    }
    return(res)
  }
  measures$gen_jaccard <- sim.measure(
    FUN = gen_jaccard_cor, string = "generalized_jaccard_index",
    mean_scaleable = TRUE, signed = FALSE, type = "parametric"
  )
  #************************************************************
  # Mutual information
  #*************************************************************
  mutual_information <- function(x, y = NULL) {
    if (is.null(y)) {
      disc <- discretize(x)
      res <- mutinformation(disc)
    }
    else {
      res <- mutual_information(cbind(x, y))[1, 2]
    }
    return(res)
  }
  measures$mutual_information <- sim.measure(
    FUN = mutual_information, string = "mutual_information", mean_scaleable = FALSE,
    signed = FALSE, type = "non-parametric"
  )
  #*******************************************************************************
  # nc.score
  #******************************************************************************
  measures$nc.score <- sim.measure(FUN = function(x, y = NULL) {
    nc.score(x, y)
  }, string = "nc_score", mean_scaleable = FALSE, signed = TRUE, type = "non-parametric")
  #*******************************************************************************
  # Euclidean distance
  #******************************************************************************
  euclidean_similarity <- function(x, y = NULL) {
    if (is.null(y)) {
      ones <- rep(1, dim(x)[2])
      res <- ones %*% t(ones) - as.matrix(designdist(t(x), method = "((A+B-2*J)/P)/(1+(A+B-2*J)/P)", terms = "quadratic"))
    }
    else {
      res <- euclidean_similarity(cbind(x, y))[1, 2]
    }
    return(res)
  }
  measures$euclidean <- sim.measure(
    FUN = euclidean_similarity, string = "squared_euclidean", mean_scaleable = TRUE,
    signed = FALSE, type = "parametric"
  )
  #*******************************************************************************
  # Cosine distance
  #******************************************************************************
  cosine_similarity <- function(x, y = NULL) {
    if (is.null(y)) {
      # Note: As this is already a similarity, we do not need to subtract from a
      # matrix of ones
      res <- as.matrix(designdist(t(x), method = "J/sqrt(A*B)", terms = "quadratic"))
    }
    else {
      res <- cosine_similarity(cbind(x, y))[1, 2]
    }
    return(res)
  }
  measures$cosine <- sim.measure(FUN = cosine_similarity, string = "cosine", signed = TRUE, type = "parametric")
  if(is.null(subset)){
    subset <- c("pearson", "spearman", "kendall", "bray_curtis", "jaccard",
                "gen_jaccard", "mutual_information", "nc.score", "euclidean",  "cosine")
  }
  measures <- measures[subset]
  return(measures)
}

#' @title create_ccrepe_jobs
#' @description  Creates ccrepe jobs out of sim.score functions
#' @inheritParams runAnalysis
create_ccrepe_jobs <- function(sim.scores = similarity_measures(), prefix =
                                 "significant_interactions", postfix = ".csv") {
  jobs <- lapply(sim.scores, function(sim.score) list(
      ccrepe_args =list(
        sim.score = sim.score@FUN
      ),
      output_args = list(
        filename = paste0(prefix, "_", sim.score@string, postfix),
        score_attributes = sim.measure.attributes(sim.score)
      )
    ))
}
#*************************************************************
# Creates similarity measure functions which adds random noise
# at certain magnitude to the data before calculating
# the similarity score
#*************************************************************
#' @title Noisify
#'
#' @description
#' Creates similarity measure functions which adds random noise
#' at certain magnitude to the data before calculating
#' the similarity score
#'
#' @param sim.scores The list of objects of class \link{sim.measure} to nosify
#'
#' @param magnitude The magnitude of the noise to add. This is: The standard deviation for normal noise and the
#' radius for the interval for uniform noise
#'
#' @param noise The type of noise to add.
#' \code{'none'} just returns the similarity measure without doing anything.
#' \code{'uniform'} adds uniformly distributed noise
#' \code{'normal'} adds normaly distributed  noise
#' Any combination of the options can be passes as a character vector.
#' In this case, sim.score functions are returned for the selected types of noise.
#'
#' @return A list of \link{sim.measure} objects corresponding to the input where noise have been
#' added
#'
#' @importFrom stats rnorm runif
#'
#' @export
noisify <- function(sim.scores = mean_scale(), magnitude = 1e-5, noise = c("none", "uniform", "normal")) {
  res <- list()
  if ("none" %in% noise) {
    res <- c(res, sim.scores)
  }
  # Keeps only the similarity measures which can be noisified (presence-absence
  # detecting functions excluded)
  is.noisifiable <- lapply(sim.scores, function(x) !x@categorical)
  sim.scores <- sim.scores[unlist(is.noisifiable)]
  if ("uniform" %in% noise) {
    noiseFUN <- function(n) {
      runif(n, min = -magnitude, max = magnitude)
    }
    uniform_functions <- lapply(
      sim.scores,
      function(sim.score) sim.measure(
          FUN = noisificationTemplate(sim.score@FUN, noiseFUN),
          string = paste0(sim.score@string, "_uniform"),
          mean_scaleable = FALSE, signed = sim.score@signed,
          type = sim.score@type
        )
    )

    names(uniform_functions) <- paste0(names(sim.scores), "_uniform")
    res <- c(res, uniform_functions)
  }
  if ("normal" %in% noise) {
    noiseFUN <- function(n) {
      rnorm(n, mean = 0, sd = magnitude)
    }
    normal_functions <- lapply(
      sim.scores,
      function(sim.score) sim.measure(
          FUN = noisificationTemplate(sim.score@FUN, noiseFUN),
          string = paste0(sim.score@string, "_normal"),
          mean_scaleable = FALSE, signed = sim.score@signed,
          type = sim.score@type
        )
    )
    names(normal_functions) <- paste0(names(sim.scores), "_normal")
    res <- c(res, normal_functions)
  }
  return(res)
}

#*************************************************************
# Helper function for nosify, creates the acutual noised
# functions
#*************************************************************
noisificationTemplate <- function(FUN, noiseFUN) {
  addNoise <- function(x) {
    abs(x + noiseFUN(n = length(x)))
  }
  noisifiedFunction <- function(x, y = NULL)
    if (is.null(y)) {
      noisedFrame <- apply(x,
        MARGIN = 2,
        FUN = addNoise
      )
      return(FUN(noisedFrame))
    }
    else {
      noisedX <- addNoise(x)
      noisedY <- addNoise(y)
      return(FUN(noisedX, noisedY))
    }
}

#' @title mean_scale
#'
#' @param sim.scores A list of \code{sim.measure} objects
#'
#' @param append Logical, should the \code{sim.scores} passed to the function also be returned?
#'
#' @description Creates \code{sim.measure} objects which devide each compontent of by its mean before performing
#' the calculations
#' @export
mean_scale <- function(sim.scores = similarity_measures(), append = TRUE) {
  scaled_scores <- list()
  is.mean_scalable <- lapply(sim.scores, function(x) x@mean_scaleable)
  measures_to_scale <- sim.scores[unlist(is.mean_scalable)]
  scaled_scores <- lapply(
    measures_to_scale,
    function(sim.score) sim.measure(
        FUN = scalationTemplate(sim.score@FUN),
        string = paste0(sim.score@string, "_scaled"),
        mean_scaleable = FALSE,
        signed = sim.score@signed,
        type = sim.score@type
      )
  )
  names(scaled_scores) <- paste0(names(measures_to_scale), "_scaled")
  if (append) {
    return(c(sim.scores, scaled_scores))
  }
  else {
    return(scaled_scores)
  }
}

scalationTemplate <- function(FUN) {
  function(x, y = NULL) {
    if (is.null(y)) {
      x <- x / colMeans(x)
      res <- FUN(x)
    }
    else {
      x <- x / mean(x)
      y <- y / mean(y)
      res <- FUN(x, y)
    }
    return(res)
  }
}
