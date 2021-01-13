#' @title Solve OTU equation systems by ridge regularization
#'
#' @description When determining the Lotka-Volterra coefficients of an OTU system, this function solves the linear
#' equations using ridge regularization (l^2).
#'
#'
#'
#' @param equations Equations returned from \link{integralSystem}
#'
#' @param weights Named numeric, a vector with two entries. The first one is named \code{'self'} and represents the
#' magnitude of penalizing the maximal growth rate of an OTU. The second is named \code{'interaction'} and represents
#' the magnitude of penalizing the interactions between the OTUs. Ideally, the weights should be found
#' by cross-validation
#'
#' @details
#'
#' The equations are fitted individually, this is: For each OTU, all coefficients determining it abundance
#' is fitted independently from the other coeffcients. The entity
#'
#' \deqn{\| A\beta-b \|^2+\lambda_{\mathrm{self}}\beta_{\mathrm{self}}^2+
#' \lambda_{\mathrm{interaction}}\|\beta_{\mathrm{interaction}}\|^2}
#' will be minimized where \eqn{\beta} is the estimated solution, \eqn{\beta_{\mathrm{self}}} is the estimated
#' maximal growth rate of the OTU and \eqn{\beta_{\mathrm{inteaction}}} are the estimated effect of the interactions
#' with the other OTUs (and itself).
#'
#'
#' @return
#'
#' An n_OTU times (n_OTU+1) numeric matrix where each row contains the coefficients determining
#' the growth of the selected OTU. The first column contains the estimated growth rates, while the other
#' rows contain the estimates for the interaction coefficients. This object has
#' the S3 class attribute \code{LV}.
#'
#' @references
#'
#' Richard R. Stein et al. ``Ecological Modeling from Time-Series Inference: Insight into Dy-
#' namics and Stability of Intestinal Microbiota.'' In: \emph{PLoS Comput. Biol. 9.12 (desember
#' 2013)}. issn: 1553-7358. DOI: \url{https://doi.org/10.1371/journal.pcbi.1003388}.
#'
#'
#'
#' @seealso \code{\link{integralSystem}}
#' @examples
#' library(micInt)
#' library(phyloseq)
#' data("seawater")
#' subsetted_seawater <- subset_samples(seawater, Reactor == 2)
#' systems <- integralSystem(OTU_time_series(subsetted_seawater,"Week"), kind = "log_integral"))
#' fit <- ridge_fit(systems,weights = c(self = 10,interaction = 1))
#' @export

ridge_fit <- function(equations, weights) {
  n_coefficients <- ncol(equations[[1]]$A)
  # Ridge regularization matrix
  ridge_matrix <- c(weights["self"], rep(weights["interaction"], n_coefficients - 1)) %>% diag()
  solutions <- lapply(equations, function(equation) {
    solution <- tryCatch(solve(t(equation$A) %*% equation$A + ridge_matrix, t(equation$A) %*% equation$b),
      error = function(e) {
        stop("Equation system singular")
      }
    )
    solution %>% matrix(nrow = 1)
  })
  solution_matrix <- do.call(rbind, solutions)
  rownames(solution_matrix) <- names(solutions)
  colnames(solution_matrix) <- c("self", names(solutions))
  class(solution_matrix) <- c("LV", class(solution_matrix))
  return(solution_matrix)
}
