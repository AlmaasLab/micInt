#' @title Make Lotka-Volterra coefficient system
#' @description This function accepts time series of OTU measurements
#' and the linear equation system for its coefficients, as described in cited article.
#'
#' @param time_series An \code{OTU_time_series} or a list of them
#' If a list is passed, the OTUs must be the same and in the same order in each of the sub-objects.
#'
#' @param kind Charachter, one of \code{c('integral','log_integral')} Choose whether use the integral or log-integral approach
#' described in the article cited.
#'
#' @param removeZeros Logical, for the \code{'integral'} option, remove the equations which
#' produce zero rows. This is ignored for the \code{'log_integral'} because in that
#' case, all two consecutive measurements which contains a zero for the OTU in question,
#' must be discared anyway.
#'
#'
#' @return A list of lists (one for each OTU) with two elements:
#' \itemize{
#' \item \code{A} The coefficient matrix for the system
#' \item \code{b} The right side of the system
#' }
#'
#' @references
#' P. H. Kloppers and J. C. Greeff. ``Lotka-Volterra model parameter
#' estimation using experiential data''. In: \emph{Appl. Math. Comput. 224}
#' (Nov. 2013), pp. 817–825. ISSN: 0096-3003. DOI: \url{https://doi.org/10.1016/j.amc.2013.08.093}
#' @examples
#' library(micInt)
#' library(phyloseq)
#' data("seawater")
#' subsetted_seawater <- subset_samples(seawater, Reactor == 2)
#' integralSystem(OTU_time_series(subsetted_seawater,"Week"), kind = "log_integral")
#' @export
integralSystem <- function(time_series, kind = "integral", removeZeros = TRUE) {
  if (inherits(time_series,'OTU_time_series')) {
    time_series <- list(time_series)
  }
  if (kind == "integral") {
    equations <- lapply(time_series, function(x) extract_integral_system(x, removeZeros = removeZeros))
  }
  else {
    equations <- lapply(time_series, extract_log_integral_system)
  }
  # Now we have lists of equation systems for each time series. We now want to combine them
  stacked_equations <- stack_equations(equations)
  return(stacked_equations)
}


extract_integral_system <- function(time_series, removeZeros) {
  force(time_series)
  timediff <- diff(time_series@time_points)
  # Note: OTUs in columns, data points in rows as in the refined table
  OTU_matrix <- as.matrix(time_series@table)
  OTU_diff <- diff(OTU_matrix)
  n_timesteps <- length(timediff)
  n_OTUs <- ncol(OTU_matrix)
  # The matrix holding the OTU abundances at the start of each timestep
  start_matrix <- OTU_matrix[1:n_timesteps, ]
  # The matrix holding the OTU abundances at the end of each timestep
  end_matrix <- OTU_matrix[1:n_timesteps + 1, ]
  # Adds the constant column to each of the matrices
  ones <- rep(1, n_timesteps)
  start_matrix <- cbind(ones, start_matrix)
  end_matrix <- cbind(ones, end_matrix)
  # Tolerance for determining zero rows
  tol <- 2 * .Machine$double.eps
  # Finds the system for each OTU
  systems <- lapply(seq_len(n_OTUs), FUN = function(i) {
    # The coefficient matrix of the system
    # Notice that the first column of start_matrix and end_matrix are columns with ones,
    # so we have to increment the index by one when adressing those matrices
    A <- (diag(timediff) / 2) %*% (diag(start_matrix[, i + 1]) %*% start_matrix +
      diag(end_matrix[, i + 1]) %*% end_matrix)
    colnames(A) <- c("self", colnames(OTU_matrix))
    # The right side of the system
    b <- OTU_diff[, i]
    if (removeZeros) {
      nonZeroRows <- which(apply(X = A, MARGIN = 1, FUN = function(row) any(abs(row) > tol)))
      A <- A[nonZeroRows, ]
      b <- b[nonZeroRows]
    }
    return(list(A = A, b = b))
  })
  names(systems) <- colnames(OTU_matrix)
  return(systems)
}

extract_log_integral_system <- function(time_series) {
  timediff <- diff(time_series@time_points)
  # Note: OTUs in columns, data points in rows as in refined table
  OTU_matrix <- as.matrix(time_series@table)
  # The matrix describing the changes between consequtive time points
  OTU_diff <- diff(OTU_matrix)
  # The matrix describing the change in log abundances between consequtive time points
  # Note that this produces infinite and NaN values when taking the logarithm of zeros
  log_OTU_diff <- suppressWarnings(
    diff(log(OTU_matrix))
  )
  n_OTUs <- ncol(OTU_matrix)
  n_timesteps <- length(timediff)
  # The coefficient matrix is the same in this case, assuming we take all equations into
  # account. We have to exclude some of them later due to the NaN-value on the right side
  twos <- rep(2, n_timesteps)
  A_full <- (diag(timediff) / 2) %*% cbind(
    twos,
    OTU_matrix[seq_len(n_timesteps), ] +
      OTU_matrix[seq_len(n_timesteps) + 1, ]
  )
  colnames(A_full) <- c("self", colnames(OTU_matrix))
  # Finds the system for each OTU
  systems <- lapply(seq_len(n_OTUs), FUN = function(i) {
    # Finds the applicable equations (where the right side is finite)
    valid_equations <- which(is.finite(log_OTU_diff[, i]))
    # The right side of the system
    b <- log_OTU_diff[valid_equations, i]
    # The coefficient matrix of the system
    A <- A_full[valid_equations, ]
    return(list(A = A, b = b))
  })
  names(systems) <- colnames(OTU_matrix)
  return(systems)
}

#' @title Stack equations from different time series
#'
#' @param equations A list of lists of equations (coefficient matrix
#' \code{A} and right side \code{b}), where the second list layer contains the
#' name of the respective OTUs
#'
#' @description For each OTU, all equations will be stacked together into a single
#' system
#'
stack_equations <- function(equations) {
  # First, we find all available OTU names
  OTU_names <- lapply(equations, names) %>% unlist() %>% unique()
  stacked_equations <- lapply(OTU_names, FUN = function(name) {
    equations_to_stack <- lapply(equations, function(x) x[[name]])
    right_sides_to_stack <- lapply(equations_to_stack, function(x) x$b)
    matrices_to_stack <- lapply(equations_to_stack, function(x) x$A)
    stacked_right_sides <- do.call(c, right_sides_to_stack)
    stacked_matrices <- do.call(rbind, matrices_to_stack)
    list(A = stacked_matrices, b = stacked_right_sides)
  })
  names(stacked_equations) <- OTU_names
  return(stacked_equations)
}
