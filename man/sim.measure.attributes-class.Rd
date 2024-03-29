% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_measure_class.R
\docType{class}
\name{sim.measure.attributes-class}
\alias{sim.measure.attributes-class}
\alias{sim.measure.attributes}
\alias{sim_measure_function}
\alias{sim.measure.attributes,sim.measure-method}
\alias{sim_measure_function,sim.measure-method}
\alias{sim.measure.attributes,character-method}
\title{Attributes of similarity measures}
\usage{
sim.measure.attributes(x, ...)

sim_measure_function(x, ...)

\S4method{sim.measure.attributes}{sim.measure}(x, ...)

\S4method{sim_measure_function}{sim.measure}(x, ...)

\S4method{sim.measure.attributes}{character}(x, type_measure, signed)
}
\arguments{
\item{x}{An \link{sim.measure}}

\item{...}{Ignored}

\item{type_measure}{Character, the type of similarity measure which should be repesented}

\item{signed}{Logical, is the similarity measure signed?}
}
\description{
Accessor and construction functions for dealing with the
\code{sim.measure} attributes of similarity measures. The \code{sim.measure.attributes}
extracts or constructs an object showing the attributes, whereas \code{sim_measure_function} extracts
the similarity function itself.
}
\section{Slots}{

\describe{
\item{\code{string}}{Character, the similarity score's human readable name}

\item{\code{type_measure}}{Character, what type to similarity measure is this (equvialent to \code{type} field
in \link{sim.measure})}

\item{\code{signed}}{Logical, does this similarity measure return signed results?}
}}

\examples{
# This measures how well the arguments agree with the signs
# This primary function accept two vectors only
my_measure_simple <- function(x,y){
 (sum((x > 0) & (y > 0)) + sum((x < 0) & (y < 0))) /
  (sum(x > 0) + sum(x < 0))
}
# We extend this function into acception a matrix argument when y is now supplied
my_measure <- function(x,y=NULL){
if(is.null(y)){
n_vectors <- ncol(x)
res <- matrix(NA_real_,n_vectors,n_vectors)
for (i in seq_len(n_vectors)){
    for(j in seq_len(n_vectors)){
    res[i,j] <- my_measure_simple(x[,i],x[,j])
    }
}
    res

}
else{
my_measure_simple(x,y)
}
}

my_sim_score <-sim.measure(my_measure,string= "My measure", categorical = FALSE,
 mean_scaleable = FALSE, signed = FALSE, type = "non-paramtric")
sim_fun <- sim_measure_function(my_sim_score)
sim_fun(c(1,0,-1),c(2,-1,-1))
sim.measure.attributes(my_sim_score)
}
