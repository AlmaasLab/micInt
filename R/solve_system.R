#' @title Predict trajectory of a Lotka-Volterra system from estimated
#' coefficients
#'
#' @description This function utilizes the \code{\link{deSolve}} package in order to
#' simulate a Lotka-Volterra system with fitted coefficients. Hence, it can be used
#' for predictions.
#'
#' @param fit An object of class \code{LV_fit} returned from a suitable
#' fitting function such as \code{\link{ridge_fit}}, containing the Lotka-Volterra coefficients
#'
#'
#' @param start The vector of the microbial community at the start
#' of the simulation. The OTUs must be in the same order as in the
#' \code{fit} argument.
#'
#' @param times The time points to include in the simulation
#'
#' @return Matrix of class \code{deSolve}. For more information, see the description in \link{ode}.
#'
#'
#' @seealso \link{ridge_fit}, \link{deSolve-package}
#'
#'
#'
#'
#' @export
predict.LV=function(fit,start,times){
ode(y=start,func=system_equation,parms = fit,times = times)
}

#' @title The ODE function for the Lotka-Volterra system
system_equation=function(t,y,params){
self_growth_rate=param[,1]
interaction_matrix = params[,-1]
# The right hand side as it would be if devided by the concentration of
# each OTU beforehand
right_side =self_growth_rate+rowSums(interaction_matrix%*%y)
# Finally, multiply with the concentraions
dxdt=y*right_side
}
