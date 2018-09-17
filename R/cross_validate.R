#' @title Measure fit of Lotka-Volterra coefficients
#'
#' @description When using cross-validation to find the Lotka-Volterra
#' coefficients, this function is used to determine the error of the cross-validation. It does
#' so by applying the
#'
#' @param test_equations
#' The list of list of equations to test, should be independent of the \code{solution_matrix}.
#' For more information, see the description of \code{test_equation} in \link{ridge_fit}
#'
#' @param solution_matrix
#' The proposed solution of which the fitness is found, such as the ones returned from
#' \link{ridge_fit}
#'
#'
#' @details
#' The errors are combined and calculated among all differences between the actual right
#' sides of the equations and the predicted ones
#'
#' @return
#' A list of two:
#'  \itemize{
#' \item \code{'RMSE'} The root mean square error of the right sides
#' \item \code{'MAE'} The mean absolute error of the right sides
#' }
#'
#'

test_LV_fit=function(test_equations,solution_matrix){
n_OTUs = length(test_equations)
errors = lapply(1:n_OTUs, function(i){
equation = test_equations[[i]]
equation$A%*%solution_matrix[i,]-equation$b
}
)
cat_errors = do.call(c,errors)
list(
  RMSE = sqrt(mean(cat_errors^2)),
MAE = mean(abs(cat_errors))
)
}


#' @title
#' Cross-validate fit for Lotka-Volterra coefficients
#'
#' @description
#' The function preforms k-fold cross-validation on a grid of regularization parameters
#'
#' @param equations
#'
#' A list of equation systems returned from \code{\link{integralSystem}}
#'
#' @param n_folds
#'
#' The number of folds to apply
#'
#' @param weights
#'
#' A n times 2 data.frame containing the weights to be tried out. The two colmns are: \itemize{
#' \item \code{'self'} The regularization for the maximal growth rate
#' \item \code{'interaction'} The regulartization for the interactions between OTUs
#' }
#' If missing, the cross validation is carried out over a quadratic grid from 0 to 10 in steps of 0.1
#'
#'
#' @return
#'
#' A data.frame containing the parameter values, and the two columns:
#'  \itemize{
#' \item \code{'RMSE'} The root mean square error of the right sides
#' \item \code{'MAE'} The mean absolute error of the right sides
#' }
#'
#'
#'
#' @export
cv.LV_fit=function(equations,n_folds,weigths = expand.grid(self=0.1*0:100,interaction=0.1*0:100)){
# We first find the number of equation in each equation system
n_equations = lapply(equations, function(x) nrow(x$A))
# the total number of systems
n_systems = length(equations)
# and the number of parameter combinations to test
n_combinations = nrow(weights)
list_weights = lapply(1:n_combinations,function(i) weigths[i,])
errors = lapply(list_weights,FUN = function(weights){
# At this level, the weights are fixed and we assign the systems
# into different folds
names(weigths) = c('self','interaction')
fold = sample(1:n_folds,size = n_systems,replace = TRUE)
number_in_fold = sapply(1:n_folds, function(i) sum(fold==i))
lapply(1:n_folds,function(i){
train = equations[fold != i]
test = equations[fold == i]
fit = ridge_fit(train,weigths)
CV_res = as.data.frame(test_LV_fit(test_equations = test,solution_matrix = fit))
})
summary_statistics = do.call(rbind,errors)
RMSE = summary_statistics$RMSE^2*number_in_fold %>% mean %>% sqrt
MAE = summary_statistics$MAE*number_in_fold %>% mean
as.data.frame(list(RMSE=RMSE,MAE=MAE))
}
)
results = do.call(rbind, errors)
return(cbind(weigths,results))
}
