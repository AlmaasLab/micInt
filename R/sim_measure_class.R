#' @name sim.measure-class
#' @title Similarity measure
#' @description
#' A similarity measure structure with the function ifself, and metadata fields
#'
#' @slot FUN The similarity measure function to use, of the form \code{f(x,y)}.
#' The function must satisfy the following: \itemize{
#' \item If \code{x} and \code{y} are both given, they are treated as vectors and the
#' similarity score between them is returned
#' \item If only \code{x} is given, it is treated as a matrix (or data frame) where the
#' features are in columns and sample in the rows. The function must then return
#' the matrix of pairwise similarities
#' \item The function must return a similarity score in the interval [0,1], where 1 is
#' perfect similarity and 0 is perfect dissimilarity. Alternativly, signed similarity scores
#' in the interval [-1,1] can also be used.
#' }
#' @slot string The similarity score's human readable name.
#'
#' @slot categorical Logical value indicating whether the function only
#' considers presence-absence of an OTU. If \code{TRUE}, the function will not be
#' noisified as there is no point in adding noise to the data.
#'
#'
#' @slot mean_scaleable Logical value telling whether it makes sense to scale the input by its
#' mean before applying this function
#'
#' @slot signed If \code{TRUE}, the similarity measures gives signed results in the interval [-1,1],
#' else they are in the interval [0,1].
#'
#' @slot type Character; the way the similarity measure processes the information. On of: \itemize{
#' \item \code{"presence-absence"} Measures the similarity of OTUs solely by their presence-absence
#' pattern
#' #' \item \code{"non-parametric"} Measures the similarity based on the rank pattern, but besides this,
#' no other information on the abudandances are used
#' \item \code{"parameteric"} The actual values of the abundances matter in the calulations
#' }
#'

#'
#'
#' @seealso noisify
#'
#' @export
setClass(Class = "sim.measure", slots = c(
  FUN = "function", string = "character",
  categorical = "logical", mean_scaleable = "logical", signed = "logical", type =
    "character"
))

#' @rdname sim.measure-class
#' @param FUN The similarity measure function to use, of the form \code{f(x,y)}.
#' The function must satisfy the following: \itemize{
#' \item If \code{x} and \code{y} are both given, they are treated as vectors and the
#' similarity score between them is returned
#' \item If only \code{x} is given, it is treated as a matrix (or data frame) where the
#' features are in columns and sample in the rows. The function must then return
#' the matrix of pairwise similarities
#' \item The function must return a similarity score in the interval [0,1], where 1 is
#' perfect similarity and 0 is perfect dissimilarity. Alternativly, signed similarity scores
#' in the interval [-1,1] can also be used.
#' }
#' @param string The similarity score's human readable name.
#'
#' @param categorical Logical value indicating whether the function only
#' considers presence-absence of an OTU. If \code{TRUE}, the function will not be
#' noisified as there is no point in adding noise to the data.
#'
#'
#' @param mean_scaleable Logical value telling whether it makes sense to scale the input by its
#' mean before applying this function
#'
#' @param signed If \code{TRUE}, the similarity measures gives signed results in the interval [-1,1],
#' else they are in the interval [0,1].
#'
#' @param type Character; the way the similarity measure processes the information. On of: \itemize{
#' \item \code{"presence-absence"} Measures the similarity of OTUs solely by their presence-absence
#' pattern
#' \item \code{"non-parametric"} Measures the similarity based on the rank pattern, but besides this,
#' no other information on the abudandances are used
#' \item \code{"parameteric"} The actual values of the abundances matter in the calulations
#' }
#' @examples
#' # This measures how well the arguments agree with the signs
#' # This primary function accept two vectors only
#' my_measure_simple <- function(x,y){
#'  (sum((x > 0) & (y > 0)) + sum((x < 0) & (y < 0))) /
#'   (sum(x > 0) + sum(x < 0))
#' }
#' # We extend this function into acception a matrix argument when y is now supplied
#' my_measure <- function(x,y=NULL){
#' if(is.null(y)){
#' n_vectors <- ncol(x)
#' res <- matrix(NA_real_,n_vectors,n_vectors)
#' for (i in seq_len(n_vectors)){
#'     for(j in seq_len(n_vectors)){
#'     res[i,j] <- my_measure_simple(x[,i],x[,j])
#'     }
#' }
#'     res
#'
#' }
#' else{
#' my_measure_simple(x,y)
#' }
#' }
#'
#' my_sim_score <-sim.measure(my_measure,string= "My measure", categorical = FALSE,
#'  mean_scaleable = FALSE, signed = FALSE, type = "non-paramtric")
#' sim_fun <- sim_measure_function(my_sim_score)
#' sim_fun(c(1,0,-1),c(2,-1,-1))
#' sim.measure.attributes(my_sim_score)
#'
#' @export
sim.measure <- function(FUN, string, categorical = FALSE, mean_scaleable = FALSE, signed = FALSE, type = "parametric") {
  new("sim.measure",
    FUN = FUN, string = string, categorical = categorical,
    mean_scaleable = mean_scaleable, signed = signed, type = type
  )
}

#' @rdname sim.measure.attributes-class
#' @title Attributes of similarity measures
#' @description
#' Accessor and construction functions for dealing with the
#' \code{sim.measure} attributes of similarity measures. The \code{sim.measure.attributes}
#' extracts or constructs an object showing the attributes, whereas \code{sim_measure_function} extracts
#' the similarity function itself.
#' @slot string Character, the similarity score's human readable name
#' @slot type_measure Character, what type to similarity measure is this (equvialent to \code{type} field
#' in \link{sim.measure})
#' @slot signed Logical, does this similarity measure return signed results?
#' @examples
#' #' @examples
#' # This measures how well the arguments agree with the signs
#' # This primary function accept two vectors only
#' my_measure_simple <- function(x,y){
#'  (sum((x > 0) & (y > 0)) + sum((x < 0) & (y < 0))) /
#'   (sum(x > 0) + sum(x < 0))
#' }
#' # We extend this function into acception a matrix argument when y is now supplied
#' my_measure <- function(x,y=NULL){
#' if(is.null(y)){
#' n_vectors <- ncol(x)
#' res <- matrix(NA_real_,n_vectors,n_vectors)
#' for (i in seq_len(n_vectors)){
#'     for(j in seq_len(n_vectors)){
#'     res[i,j] <- my_measure_simple(x[,i],x[,j])
#'     }
#' }
#'     res
#'
#' }
#' else{
#' my_measure_simple(x,y)
#' }
#' }
#'
#' my_sim_score <-sim.measure(my_measure,string= "My measure", categorical = FALSE,
#'  mean_scaleable = FALSE, signed = FALSE, type = "non-paramtric")
#' sim_fun <- sim_measure_function(my_sim_score)
#' sim_fun(c(1,0,-1),c(2,-1,-1))
#' sim.measure.attributes(my_sim_score)
#' @export
setClass(Class = "sim.measure.attributes", slots = c(string = "character", type_measure = "character", signed = "logical"))

#' @rdname sim.measure.attributes-class
#' @export
setGeneric("sim.measure.attributes", def = function(x, ...) standardGeneric("sim.measure.attributes"))
#' @rdname sim.measure.attributes-class
#' @export
setGeneric("sim_measure_function", def = function(x,...) standardGeneric("sim_measure_function"))

#' @rdname sim.measure.attributes-class
#' @param x An \link{sim.measure}
#' @export
setMethod("sim.measure.attributes", "sim.measure",
  definition = function(x,...) {
    new("sim.measure.attributes", string = x@string, type_measure = x@type, signed = x@signed)
  }
)

#' @rdname sim.measure.attributes-class
#' @export
setMethod("sim_measure_function","sim.measure",definition = function(x,...)
            x@FUN
            )

#'
#' @rdname sim.measure.attributes-class
#' @param type_measure Character, the type of similarity measure which should be repesented
#' @param signed Logical, is the similarity measure signed?
#' @param ... Ignored
#' @export
setMethod(
  "sim.measure.attributes", "character",
  function(x, type_measure, signed) {
    new("sim.measure.attributes", string = x, type_measure = type_measure, signed = signed)
  }
)
