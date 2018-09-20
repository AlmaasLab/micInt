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
# If the solution is not available, it makes no sense to calculate the statistics
if(is.na(solution_matrix)){
return(list(RMSE = NA_real_, MAE = NA_real_))
}
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
#' @param time_series
#'
#' A list of \code{\link{OTU_time_series}} to be cross-validated
#'
#' @param kind Charachter, one of \code{c('integral','log_integral')} Choose whether use the integral or log-integral approach
#' described in the article cited.
#'
#' @references
#' P. H. Kloppers and J. C. Greeff. ``Lotka-Volterra model parameter
#' estimation using experiential data''. In: \emph{Appl. Math. Comput. 224}
#' (Nov. 2013), pp. 817–825. ISSN: 0096-3003. DOI: \link{10.1016/j.amc.2013.08.093}
#'
#' @param n_folds
#'
#' The number of folds to apply. If missing, leave-one-out cross-validation is preformed.
#'
#' @param weights
#'
#' A n times 2 data.frame containing the weights to be tried out. The two colmns are: \itemize{
#' \item \code{'self'} The regularization for the maximal growth rate
#' \item \code{'interaction'} The regulartization for the interactions between OTUs
#' }
#' If missing, the cross validation is carried out over a quadratic grid from 0 to 10 in steps of 0.1
#'
#' @details
#'
#' In this setting, each time series is a ``data point'' in the cross-validation procedure, meaning that each
#' fit is preformed on some of the time series and the other are used for validation.
#'
#'
#'
#'
#' @import magrittr
#'
#' @return
#'
#' A data.frame containing the parameter values, and the two columns of cross-validation error:
#'  \itemize{
#' \item \code{'RMSE'} The root mean square error of the right sides
#' \item \code{'MAE'} The mean absolute error of the right sides
#' }
#'
#'
#'
#' @export
cv.LV_fit=function(time_series,n_folds=length(time_series),kind='integral',
                   weights = expand.grid(self=0.1*0:100,interaction=0.1*0:100)){
# We first find the number of time series
# the total number of systems
n_time_series = length(time_series)
# and the number of parameter combinations to test
n_combinations = nrow(weights)
list_weights = lapply(1:n_combinations,function(i) weights[i,]%>% as.numeric)
# We cache the systems in order to avoid re-calculating them
systems= lapply(time_series,function(x) integralSystem(x,kind=kind))
errors = lapply(list_weights,FUN = function(weights){
# At this level, the weights are fixed and we assign the time series
# into different folds
names(weights) = c('self','interaction')
fold = sample(rep(1:n_folds,length.out=n_time_series))
number_in_fold = sapply(1:n_folds, function(i) sum(fold==i))
fold_errors=lapply(1:n_folds,function(i){
train_equations = systems[fold != i] %>% stack_equations
test_equations = systems[fold == i] %>% stack_equations
# We have to consider the cases where at least one of the systems are
# singular
fit = tryCatch(ridge_fit(train_equations,weights),error = function(e) {
  NA_real_
  }
)
CV_res = as.data.frame(test_LV_fit(test_equations = test_equations,solution_matrix = fit))
}
)
summary_statistics = do.call(rbind,fold_errors)
# Note the parentesis the next two lines. Without them, the expression
# is not evaluated correctly as the multiplicator operator
# has lower precedence than the piping operator
RMSE = (summary_statistics$RMSE^2*number_in_fold) %>% mean %>% sqrt
MAE = (summary_statistics$MAE*number_in_fold) %>% mean
as.data.frame(list(RMSE=RMSE,MAE=MAE))
}
)
results = do.call(rbind, errors)
return(cbind(weights,results))
}
