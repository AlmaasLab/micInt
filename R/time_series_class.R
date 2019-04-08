#' @name OTU_time_series-class
#'
#' @title Time series class for OTU tables
#'
#' @description This class provides a way to integrate measurements of time
#' to OTU tables, useful for time series analysis
#'
#' @slot table The refined OTU table representing the community
#'
#'
#' @slot time_points Numeric, the time from the start of the
#' experiment or any other arbitrary time reference
#'
#' @param table A refined table returned form \code{\link{refine_data}} or a \code{phyloseq} object
#'
#' @param time_points A numeric vector of the start point. If a \code{phyloseq} object is provided, it can be a
#' character vector of length one, specifying which column in the sample data to use for the time points.
#'
#' @details If you provide a \code{phyloseq} object, you might want to split it into groups before
#' using the constructor function.
#'
#' @export
setClass(Class = "OTU_time_series", slots = c(table = "data.frame", time_points = "numeric"))

#' @rdname OTU_time_series-class
setGeneric('OTU_time_series',def = function(table,time_points){
  standardGeneric('OTU_time_series')
}
)

#' @rdname OTU_time_series-class
#' @export
setMethod('OTU_time_series',signature = c(table='data.frame','numeric'),definition = function(table, time_points){
  # Ensures the observations come in the correct order
  sorted_table <- table[order(time_points), ]
  res <- new("OTU_time_series", table = sorted_table, time_points = sort(time_points))
  stopifnot(validObject(res))
  return(res)
})

#' @rdname OTU_time_series-class
#' @export
setMethod('OTU_time_series',signature = c(table='phyloseq','numeric'),
          definition = function(table, time_points){
  table_refined <- phyloseq::otu_table(table) %>% data.frame()
  if (phyloseq::taxa_are_rows(table)) {
    table_refined <- t(table_refined) %>% data.frame()
  }
  rownames(table_refined) <- phyloseq::sample_names(table)
  colnames(table_refined) <- phyloseq::taxa_names(table)
  OTU_time_series(table_refined,time_points = time_points)
})

#' @rdname OTU_time_series-class
#' @export
setMethod('OTU_time_series',signature = c(table='phyloseq','character'), definition = function(table,time_points){
if(length(time_points) != 1){
  stop('Only one column from the sample data column can be selected')
}
  time_points <- phyloseq::sample_data(table)[[time_points]] %>% as.numeric()
  OTU_time_series(table, time_points)
}
  )

setValidity("OTU_time_series", method = function(object) {
  retval <- NULL
  if (nrow(object@table) != length(object@time_points)) {
    retval <- c(retval, "Number of rows in table and number of time points differ")
  }
  if (is.null(retval)) {
    return(TRUE)
  }
  else {
    return(retval)
  }
})
