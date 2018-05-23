#' @name sim.measure
#' @title sim.measure
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
#' @seealso noisify
#'
#' @export
setClass(Class="sim.measure",slots=c(FUN="function",string="character",
                                     categorical="logical",mean_scaleable="logical",signed="logical",type=
                                       "character"))

#' @name sim.measure
#'@rdname sim.measure
sim.measure=function(FUN,string,categorical=FALSE,mean_scaleable=FALSE,signed=FALSE,type="parametric"){
  new("sim.measure",FUN=FUN,string=string,categorical=categorical,
      mean_scaleable=mean_scaleable,signed=signed,type=type)
}
#' @name sim.measure.attributes
#' @title sim.measure.attributes
#' @description
#' The useful metadata fields for a \code{sim.measure}
#' @slot string Character, the similarity score's human readable name
#' @slot type_measure Character, what type to similarity measure is this (equvialent to \code{type} field
#' in \link{sim.measure})
#' @slot signed Logical, does this similarity measure return signed results?
#' @export
setClass(Class="sim.measure.attributes",slots = c(string="character",type_measure="character",signed="logical"))

setGeneric("sim.measure.attributes",def=function(x,...) standardGeneric("sim.measure.attributes")
)

#' @rdname sim.measure.attributes
#' @export
setMethod("sim.measure.attributes","sim.measure",
          definition = function(x){
            new("sim.measure.attributes",string=x@string,type_measure=x@type,signed=x@signed)
          }
)
#' @rdname sim.measure.attributes
#' @export
setMethod("sim.measure.attributes","character",
          function(x,type_measure,signed){
            new("sim.measure.attributes",string=x,type_measure=type_measure,signed=signed)
          }
)
