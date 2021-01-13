#' @title Predict trajectory of a Lotka-Volterra system from estimated
#' coefficients
#'
#' @description This function utilizes the \code{\link{deSolve}} package in order to
#' simulate a Lotka-Volterra system with fitted coefficients. Hence, it can be used
#' for predictions.
#'
#' @param object An object of class \code{LV_fit} returned from a suitable
#' fitting function such as \code{\link{ridge_fit}}, containing the Lotka-Volterra coefficients
#'
#'
#' @param start The vector of the microbial community at the start
#' of the simulation. The OTUs must be in the same order as in the
#' \code{fit} argument.
#'
#' @param times The time points to include in the simulation
#' @param ... Additional parameters to the function \code{\link{ode}} used to solve the system
#' @details By default, the \code{lsoda} is the method used, but this can be changed by passing an \code{method} argument
#' @return An object of type \code{\link{OTU_time_series}}, where the
#' table is a matrix of type \code{\link{deSolve}}
#' @seealso \link{ridge_fit}, \link{deSolve-package}
#' @importFrom deSolve ode
#' @examples
#' library(micInt)
#' library(phyloseq)
#' data("seawater")
#' subsetted_seawater <- subset_samples(seawater, Reactor == 2)
#' systems <- integralSystem(OTU_time_series(subsetted_seawater,"Week"), kind = "log_integral"))
#' fit_1 <- ridge_fit(systems,weights = c(self = 10,interaction = 1))
#' fit_2 <- ridge_fit(systems,weights = c(self = 1,interaction = 1))
#' prediction_1 <- predict(fit_1,start = rep(1,ncol(fit)),times = c(1,2,3,4,5))
#' prediction_2 <- predict(fit_2,start = rep(1,ncol(fit)),times = c(1,2,3,4,5))
#' plot_trajectory(list(fit_1=prediction_1,fit_2 =prediction_2))
#'
#'
#' @export
predict.LV <- function(object, start, times, ...) {
  sol <- ode(y = start, func = system_equation, parms = object, times = times, ...)
  OTU_time_series(table = sol[, -1] %>% as.data.frame(), time_points = sol[, 1])
}

# @title The ODE function for the Lotka-Volterra system
# @param t Numeric, the time point a which the function is evaluated, ignored as the equation is
# autonomous
# @param y Numeric, the current state of the system
# @param parms Matrix, the Lotka-Volterra parameters. The maximum growth rate is in the first column,
# while the interaction-based parameters consitute the rest of the matrix.
system_equation <- function(t, y, parms) {
  self_growth_rate <- parms[, 1]
  interaction_matrix <- parms[, -1]
  # The right hand side as it would be if devided by the concentration of
  # each OTU beforehand
  right_side <- self_growth_rate + rowSums(interaction_matrix %*% y)
  # Finally, multiply with the concentraions
  dxdt <- y * right_side
  list(dxdt)
}
