#' @title create_prefix
#' @description Standard function to create systematic prefixes from
#' parameters. The output will be of the form
#' \code{q_crit=(critical q-value)_cutoff=(the mean abundance cutoff)_magfac=(the magnitude factor)}
#' in scientific notation.
#'
#' @param q_crit critical q-value. Ignored if not set
#' @param p_crit critical p-value, ignored if \code{q-crit} is set.
#' @param cutoff The abundance cutoff. Ignored if not set
#' @param magfac The magnitude factor. Ignored if not set
#' @param digits The number of significant digits in the resulting string
#' @return A character vector of length one
#' @examples
#' create_prefix(q_crit = 1e-5,cutoff = 1e-4,magfac = 10, digits = 3)
#' create_prefix(p_crit = 1e-5,cutoff = 1e-4,magfac = 10, digits = 5)
#' @export
create_prefix <- function(q_crit = NULL, p_crit = NULL, cutoff = NULL, magfac = NULL,
                          digits = 2) {
  if (!is.null(q_crit)) {
    prefix <- paste0("q_crit=", format(q_crit, scientific = TRUE, digits = digits))
  }
  else if (!is.null(p_crit)) {
    prefix <- paste0("p_crit=", format(p_crit, scientific = TRUE, digits = digits))
  }
  else {
    prefix <- ""
  }
  if (!is.null(cutoff)) {
    if (prefix != "") {
      prefix <- paste0(prefix, "_")
    }
    prefix <- paste0(prefix, "cutoff=", format(cutoff, scientific = TRUE, digits = digits))
  }
  if (!is.null(magfac)) {
    if (prefix != "") {
      prefix <- paste0(prefix, "_")
    }
    prefix <- paste0(prefix, "magfac=", format(magfac, scientific = TRUE, digits = digits))
  }
  return(prefix)
}

# @title
# Add outputargs
#
# @description
# Add  output arguments being used by \code{\link{output_ccrepe_data}} in addition to
# the ones returned from \code{\link{ccrepe_analysis}}.
#
# @param ccrepe_res A list of results from \code{\link{ccrepe}}
#
# @inheritParams output_ccrepe_data
#
add_outputargs <- function(ccrepe_res, taxonomy = NULL, output.file = NULL, return.value = NULL,
                           threshold.type = NULL, threshold.value = NULL,
                           csv_option = NULL, removeDuplicates =  NULL) {
  ccrepe_res <- lapply(ccrepe_res, function(i) {
    n <- names(i)
    # Renames the list elements such that entity name "res" becomes "data"
    if ("res" %in% n) {
      names(n) <- n
      n["res"] <- "data"
      names(i) <- n
      return(i)
    }
  })
  output_defaultargs <- list(
    taxonomy = taxonomy, output.file = output.file,
    return.value = return.value, threshold.value = threshold.value,
    threshold.type = threshold.type,
    csv_option = csv_option, removeDuplicates = removeDuplicates
  )
  outputargs <- lapply(ccrepe_res, function(x) modifyList(x, output_defaultargs))
}
