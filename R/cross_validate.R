# @title Measure fit of Lotka-Volterra coefficients
#
# @description When using cross-validation to find the Lotka-Volterra
# coefficients, this function is used to determine the error of the cross-validation. It does
# so by applying the
#
# @param test_equations
# The list of list of equations to test, should be independent of the \code{solution_matrix}.
# For more information, see the description of \code{test_equation} in \link{ridge_fit}
#
# @param solution_matrix
# The proposed solution of which the fitness is found, such as the ones returned from
# \link{ridge_fit}
#
#
# @details
# The errors are combined and calculated among all differences between the actual right
# sides of the equations and the predicted ones
#
# @return
# A list of two:
#  \itemize{
# \item \code{'RMSE'} The root mean square error of the right sides
# \item \code{'MAE'} The mean absolute error of the right sides
# }
#
#

test_LV_fit <- function(test_equations, solution_matrix) {
  # If the solution is not available, it makes no sense to calculate the statistics
  if (is.na(solution_matrix) %>% any()) {
    return(list(RMSE = NA_real_, MAE = NA_real_))
  }
  n_OTUs <- length(test_equations)
  errors <- lapply(1:n_OTUs, function(i) {
    equation <- test_equations[[i]]
    equation$A %*% solution_matrix[i, ] - equation$b
  })
  cat_errors <- do.call(c, errors)
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
#' (Nov. 2013), pp. 817–825. ISSN: 0096-3003. DOI: \url{https://doi.org/10.1016/j.amc.2013.08.093}
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
#' @param show_progress
#'
#' Logical, if \code{TRUE}, the progress is printed for each combination of weights is is to cross-validate, ignored
#' if \code{parallel=TRUE}
#'
#' @inheritParams runAnalysis
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
#' An object of class \code{cvLV}, inheiriting \code{data.frame} containing the parameter values, and the two columns of cross-validation error:
#'  \itemize{
#' \item \code{'RMSE'} The root mean square error of the right sides
#' \item \code{'MAE'} The mean absolute error of the right sides
#' }
#'
#' @examples
#' library(micInt)
#' library(phyloseq)
#' library(magrittr)
#' data("seawater")
#' physeq_list <- subdivide_by_environment(seawater,"Reactor")
#' time_series <- lapply(physeq_list$phyloseq,OTU_time_series,
#' time_points ="Week")
#' systems <- integralSystem(time_series,kind = "integral")
#' cv_res <- cv.LV(time_series,n_folds = 3, kind = "integral",
#' weights = expand.grid(self= c(1,2),
#' interaction = c(1,2)
#' )
#' )
#' best_parameters <- cv_res[which.min(cv_res$RMSE),c("self","interaction")] %>% unlist()
#' fit <- ridge_fit(systems,best_parameters)
#' predict(fit,start = rep(1,nrow(fit)),times = c(0,1,4,6,10,20))
#'
#'
#' @export
cv.LV <- function(time_series, n_folds = length(time_series), kind = "integral",
                  weights = expand.grid(self = 0.1 * 0:100, interaction = 0.1 * 0:100),
                  show_progress = TRUE, parallel = FALSE,
                  ncpus = getOption("micInt.ncpus", 1L), cl = NULL) {
  # We first find the number of time series
  # the total number of systems
  n_time_series <- length(time_series)
  # and the number of parameter combinations to test
  n_combinations <- nrow(weights)
  list_weights <- lapply(1:n_combinations, function(i) weights[i, ] %>% as.numeric())
  # We cache the systems in order to avoid re-calculating them
  systems <- lapply(time_series, function(x) integralSystem(x, kind = kind))
  n_weights <- length(list_weights)
  this_index <- 1
  errorFUN <- function(weights) {
    # At this level, the weights are fixed and we assign the time series
    # into different folds
    names(weights) <- c("self", "interaction")
    if(show_progress){
      message(glue::glue("Cross-validating with weight combination {this_index} of {n_weights}
                         self: {weights['self']}  interaction: {weights['interaction']}"))

    }
    fold <- sample(rep(1:n_folds, length.out = n_time_series))
    number_in_fold <- vapply(1:n_folds, function(i) sum(fold == i),
                             numeric(1))
    fold_errors <- lapply(1:n_folds, function(i) {
      train_equations <- systems[fold != i] %>% stack_equations()
      test_equations <- systems[fold == i] %>% stack_equations()
      # We have to consider the cases where at least one of the systems are
      # singular
      fit <- tryCatch(ridge_fit(train_equations, weights), error = function(e) {
        NA_real_
      })
      CV_res <- as.data.frame(test_LV_fit(test_equations = test_equations, solution_matrix = fit))
    })
    summary_statistics <- do.call(rbind, fold_errors)
    # Note the parentesis the next two lines. Without them, the expression
    # is not evaluated correctly as the multiplicator operator
    # has lower precedence than the piping operator
    RMSE <- (summary_statistics$RMSE^2 * number_in_fold) %>% mean() %>% sqrt()
    MAE <- (summary_statistics$MAE * number_in_fold) %>% mean()
    this_index <<- this_index + 1
    return(as.data.frame(list(RMSE = RMSE, MAE = MAE)))
  }
  if(!parallel){
  errors <- lapply(list_weights, FUN = errorFUN)
  }
  else{
    n_cores <- min(detectCores(), ncpus)
    if(!is.null(cl)){
      cluster <- cl
    }
    else{
      cluster <- makeCluster(n_cores)
    }
    clusterEvalQ(cluster, {
      require(micInt)
      require(magrittr)
      RhpcBLASctl::blas_set_num_threads(1L)
      }
      )
    clusterExport(cl = cluster,
                  varlist = c("systems","this_index","show_progress","n_weights","n_folds",
                              "n_time_series"), envir = environment())
      errors <-
        tryCatch({
          parLapply(
            cl = cluster, X = list_weights,
            fun = errorFUN)
          },
          finally = {
            # Makes sure the cluster shuts down even though an error has occured
            stopCluster(cluster)
          }
        )
    }
  results <- do.call(rbind, errors)
  result_frame <- cbind(weights, results)
  class(result_frame) <- c("cvLV",class(result_frame))
  return(result_frame)
  }



#' @title Create cross-validation colorplot
#'
#' @description This function views the cross-validation error in a colorplot as a function
#' of the regularization weights. Hence, this approach is suitable to detect whether the cross-validation
#' procedure contains a reasonable optimum.
#'
#' @param object A \code{cvLV} object returned from \code{\link{cv.LV}}
#' @param target The cross-validation target to pick, one of \code{'RMSE'} or \code{'MAE'}
#' @param ... other arguments passed to methods
#'
#' @importFrom rlang .data
#' @return A \link{ggplot} colorplots cross-validation errors for the different weights combination
#' @examples
#' library(micInt)
#' library(phyloseq)
#' data("seawater")
#' physeq_list <- subdivide_by_environment(seawater,"Reactor")
#' time_series <- lapply(physeq_list$phyloseq,OTU_time_series,
#' time_points ="Week")
#' cv_res <- cv.LV(time_series,n_folds = 3,
#' kind = "integral",
#' weights = expand.grid(self= c(1,2),
#' interaction = c(1,2)
#' )
#' )
#' autoplot(cv_res, target = "MAE")
#' @export
autoplot.cvLV <- function(object,target = 'RMSE',...){
  ggplot(object,mapping=aes(x=.data$self,y=.data$interaction,z=!! rlang::sym(target)))+
    geom_raster(aes(fill=!! rlang::sym(target)))+scale_fill_gradientn(colours= viridis::viridis(10))
}
