#*****************************************************************************
# Analyses the data with ccrepe using different similarity measurements
# The output is by default a cvs file displaying the different significant
# OTU correlations, their correlations, p-value, q-value and taxonomy.
#****************************************************************************
#' @title
#' ccrepe_analysis
#'
#' @description
#' A wrapper around the \code{ccrepe} function, provides parallel analysis
#'
#' @param ccrepe_job A list of jobs to be passed to \code{ccrepe}. The list themselves are named lists with the arguments
#' being passed to \link{\code{ccrepe}}
#'
#' @param parallel Should the jobs be run in parallel?
#'
#' @param verbose Should the function display how much time it spent?
#'
#' @return  A list of the results of the various jobs.
#' Each element of this list containings the ccrepe results in addition to the outputargs passed as the jobs
#'
#' @seealso
#'
#' \link{ccrepe}
#'
#' @import ccrepe
#' @import parallel
#' @importFrom utils modifyList
#' @export
ccrepe_analysis=function(ccrepe_job,
                         parallel=TRUE,verbose=TRUE){
if(parallel){
n_cores=detectCores()
cluster=makeCluster(n_cores)
clusterEvalQ(cluster,require(ccrepe))
}
else{
    n_cores=1
}
ccrepe_res=list()
start=Sys.time()
# ccrepe_res=mclapply(ccrepe_job,
#                  function(x)list(res=do.call(ccrepe,x$ccrepe_args)),
#                  mc.cores = n_cores
#                  )
if(parallel){ccrepe_res=
  tryCatch(parLapply(cl = cluster,X=ccrepe_job,
                                fun=function(x)list(res=do.call(ccrepe,x$ccrepe_args))
                               )
           ,finally = {
             # Makes sure the cluster shuts down even though an error has occured
             stopCluster(cluster)
           }
           )
}
else{
  ccrepe_res=lapply(X=ccrepe_job,
                       FUN=function(x){
                         print(x$string)
                         list(res=do.call(ccrepe,x$ccrepe_args))
                         }
  )
}
stop=Sys.time()
if(verbose){
  print('Time to execute the ccrepe analysis')
  print(stop-start)
}
for (i in 1:length(ccrepe_job)){
  ccrepe_res[[i]]=modifyList(ccrepe_res[[i]],ccrepe_job[[i]]$output_args)
}
return(ccrepe_res)
}
