#*****************************************************************************
# Analyses the data with ccrepe using different similarity measurements
# The output is by default a cvs file displaying the different significant
# OTU correlations, their correlations, p-value, q-value and taxonomy.
#****************************************************************************
#' @title
#'Conduct ccrepe analysis
#'
#' @description
#' A wrapper around the \code{ccrepe} function, provides parallel analysis
#'
#' @param ccrepe_job A list of jobs to be passed to \code{ccrepe}. The lists themselves are named lists with the arguments
#' being passed to \code{\link{ccrepe}}
#'
#' @param commonargs \code{ccrepe} arguments common for all jobs
#'
#' @param parallel Should the jobs be run in parallel?
#'
#' @param verbose Should the function display how much time it spent?
#'
#' @return  A list of the results of the various jobs.
#' Each element of this list containings the ccrepe results
#' @seealso
#' \link{ccrepe}
#' @examples
#' library(micInt)
#' data("seawater")
#' sim.scores <- similarity_measures(subset= c("spearman","pearson"))
#' sim_funs <- lapply(sim.scores,sim_measure_function)
#' refined_table <- refine_data(seawater)
#' ccrepe_commonargs <- list(x = refined_table, iterations = 100,
#' memory.optimize = TRUE
#' , min.subj = 5)
#' ccrepe_job <- list(spearman=list(sim.score = sim_funs[["spearman"]]),
#' pearson = list(sim.score = sim_funs[["pearson"]]))
#' ccrepe_analysis(ccrepe_job,ccrepe_commonargs, parallel = FALSE)
#' @import parallel
#' @importFrom utils modifyList
#' @export
ccrepe_analysis <- function(ccrepe_job,commonargs,
                            parallel = TRUE, verbose = TRUE) {
  if (parallel) {
    n_cores <- detectCores()
    cluster <- makeCluster(n_cores)
    clusterEvalQ(cluster,
                 {require(micInt)
                   RhpcBLASctl::blas_set_num_threads(1L)
                   }
      )
    clusterSetRNGStream(cl=cluster)
  }
  else {
    n_cores <- 1
  }

  ccrepe_res <- list()
  start <- Sys.time()
  # ccrepe_res=mclapply(ccrepe_job,
  #                  function(x)list(res=do.call(ccrepe,x$ccrepe_args)),
  #                  mc.cores = n_cores
  #                  )
  if (parallel) {
    ccrepe_res <-
      tryCatch(
        parLapply(
          cl = cluster, X = ccrepe_job,
          fun = function(x) do.call(ccrepe, c(x,commonargs))
        ),
        finally = {
          # Makes sure the cluster shuts down even though an error has occured
          stopCluster(cluster)
        }
      )
  }
  else {
    ccrepe_res <- lapply(
      X = ccrepe_job,
      FUN = function(x) {
        # message(x$string)
        do.call(ccrepe, c(x,commonargs))
      }
    )
  }
  stop <- Sys.time()
  if (verbose) {
    message("Time to execute the ccrepe analysis")
    message(stop - start)
  }

  return(ccrepe_res)
}
