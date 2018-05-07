#*****************************************************************************
# Analyses the data with ccrepe using different similarity measurements
# The output is by default a cvs file displaying the different significant
# OTU correlations, their correlations, p-value, q-value and taxonomy.
#****************************************************************************
#' ccrepe_analysis
#'
#' @description
#' A wrapper around the \code{ccrepe} function
#'
#' @import ccrepe
#' @import tictoc
#' @import parallel
#' @importFrom utils modifyList
#' @export
ccrepe_analysis=function(ccrepe_job,
                         parallel=TRUE,verbose=TRUE){
if(parallel){
n_cores=detectCores()
cluster=makeCluster(n_cores)
clusterExport(cluster,c('ccrepe'))
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
  print(toc.outmsg(start,stop,'Time to execute the ccrepe analysis'))
}
for (i in 1:length(similarity_measures)){
  ccrepe_res[[i]]=modifyList(ccrepe_res[[i]],ccrepe_job[[i]]$output_args)
}
return(ccrepe_res)
}
