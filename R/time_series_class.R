#' @name OTU_time_series
#'
#' @title Time series class for OTU tables
#'
#' @description This class provides a way to integrate measurements of time
#' to OTU tables, useful for time series analysis
#'
#' @slot table The refined OTU table returned from \link{refine_data}
#'
#' @slot time_points Numeric, the time from the start of the
#' experiment or any other arbitrary time reference
#' @export
setClass(Class = "OTU_time_series", slots = c(table = "data.frame", time_points = "numeric"))

#' @rdname OTU_time_series
#' @export
OTU_time_series <- function(table, time_points) {
  # Ensures the observations come in the correct order
  sorted_table <- table[order(time_points), ]
  res <- new("OTU_time_series", table = sorted_table, time_points = sort(time_points))
  stopifnot(validObject(res))
  return(res)
}

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
