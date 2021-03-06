% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solve_system.R
\name{predict.LV}
\alias{predict.LV}
\title{Predict trajectory of a Lotka-Volterra system from estimated
coefficients}
\usage{
\method{predict}{LV}(object, start, times, ...)
}
\arguments{
\item{object}{An object of class \code{LV_fit} returned from a suitable
fitting function such as \code{\link{ridge_fit}}, containing the Lotka-Volterra coefficients}

\item{start}{The vector of the microbial community at the start
of the simulation. The OTUs must be in the same order as in the
\code{fit} argument.}

\item{times}{The time points to include in the simulation}

\item{...}{Additional parameters to the function \code{\link{ode}} used to solve the system}
}
\value{
An object of type \code{\link{OTU_time_series}}, where the
table is a matrix of type \code{\link{deSolve}}
}
\description{
This function utilizes the \code{\link{deSolve}} package in order to
simulate a Lotka-Volterra system with fitted coefficients. Hence, it can be used
for predictions.
}
\details{
By default, the \code{lsoda} is the method used, but this can be changed by passing an \code{method} argument
}
\examples{
library(micInt)
library(phyloseq)
data("seawater")
subsetted_seawater <- subset_samples(seawater, Reactor == 2)
systems <- integralSystem(OTU_time_series(subsetted_seawater,"Week"), kind = "log_integral")
fit_1 <- ridge_fit(systems,weights = c(self = 10,interaction = 1))
fit_2 <- ridge_fit(systems,weights = c(self = 1,interaction = 1))
prediction_1 <- predict(fit_1,start = rep(1,nrow(fit_1)),times = c(1,2,3,4,5))
prediction_2 <- predict(fit_2,start = rep(1,nrow(fit_2)),times = c(1,2,3,4,5))
plot_trajectory(list(fit_1=prediction_1,fit_2 =prediction_2))


}
\seealso{
\link{ridge_fit}, \link{deSolve-package}
}
