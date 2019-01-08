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
#' @param ... Additional parameters to the function \code{\link{ode}} used to solve the system
#' @details By default, the \code{lsoda} is the method used
#' @return An object of type \code{\link{OTU_time_series}}, where the
#' table is a matrix of type \code{\link{deSolve}}
#' @seealso \link{ridge_fit}, \link{deSolve-package}
#' @importFrom deSolve ode
#'
#'
#'
#'
#' @export
predict.LV=function(fit,start,times,diagnostics= TRUE,...){
sol=ode(y=start,func=system_equation,parms = fit,times = times)
OTU_time_series(table=sol[,-1]%>% as.data.frame,time_points = sol[,1])
}

#' @title The ODE function for the Lotka-Volterra system
system_equation=function(t,y,parms){
self_growth_rate=parms[,1]
interaction_matrix = parms[,-1]
# The right hand side as it would be if devided by the concentration of
# each OTU beforehand
right_side =self_growth_rate+rowSums(interaction_matrix%*%y)
# Finally, multiply with the concentraions
dxdt=y*right_side
list(dxdt)
}
