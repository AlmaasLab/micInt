#' @title runAnalysis
#'
#' @description
#' Runs an automized processing of the OTU table, passes the jobs to \code{ccrepe} and saves the results
#'
#' @param OTU_table The raw OTU table to be treated
#'
#' @param abundance_cutoff The mean abundance cutoff for the OTUs
#'
#' @param parallel Should the analysis be run in parallel?
#'
#' @param sim.scores The similarity measures of class \link{sim.measure}
#' to use. If it is \code{NULL}, all measures available in the package
#' will be used (recommanded for most purposes).
#'
#' @param subset Character, the subset of similarity measures to use, denoted by
#' the its name in the list (not necessarly its string) returned from \link{similarity_measures} or similarity
#' measure modiftying function such as \link{noisify}
#'If \code{NULL}, all available measures will be used
#'
#' @param file Should the tables of significant interactions be printed to a file?
#' @param returnVariables Which variables should the function return (character vector)?
#' Available options are: \itemize{
#' \item \code{similarity_measures_significance}:
#' The \code{interactions_table} of significant interactions
#' \item \code{refined_table}: The processed  OTU table
#' \item \code{min_dataset}: The smallest non-zero entity in the refined table
#' }
#' In addition, all paramerters for this function are available
#' If \code{NULL}, all variables in the function will be returned
#'
#' @param magnitude_factor When making noisified functions, the magnitude of the noise
#' will be this number multiplied with \code{min_dataset}
#'
#' @param prefix The prefix of the file names being written. Ignored if \code{file=FALSE}.
#'
#' @param metadataCols The names (character vector) or position (integer) of the
#' metadata columns to remove from the table before analyzing it
#'
#' @details
#' If the function is told to output a file and no prefix is given, the csv-files will all share a common prefix of the form:
#' \code{q_crit=(critical q-value)_cutoff=(the mean abundance cutoff)_magfac=(the magnitude factor)},
#' where all numbers are in scientific notation. Then the sim.score name follows, then the postfix and finally the csv
#' extention.
#' The postfix is by default empty.
#'
#' In order for an OTU-table to be valid, the following criteria must hold:
#'
#' \itemize{
#' \item
#' The data points (sample) are in columns, the abundances for each
#' OTU is in rows.
#' \item
#' The rows may only hold OTU abundances
#' \item
#' There may be as many metadata colums as preferable. However, the all
#' need to be declared in the \code{metadataCols} argument and the column
#' \code{taxonomy} has be there in order for the output file to contain the
#' taxonomy.
#' \item The row names of the table are the OTU names and the column names are the
#' sample names
#' }
#'
#'
#'
#' @import stringr
#' @export
runAnalysis=function(OTU_table,abundance_cutoff=1e-04,q_crit=0.05,parallel=TRUE,
                    returnVariables=NULL,subset=NULL,sim.scores=NULL,file=FALSE,magnitude_factor=10,prefix=NULL,
                    metadataCols=c('OTU Id','taxonomy'),
                    postfix=""){
if(is.null(prefix))
prefix = create_prefix(q_crit = q_crit, cutoff = abundance_cutoff, magfac = magnitude_factor)
refined_table=refine_data(OTU_table,abundance_cutoff=abundance_cutoff,metadataCols = metadataCols)
# The smallest value in the data set
min_dataset=min(apply(refined_table,MARGIN = 2,function(x) min(x[x>0])))
magnitude=magnitude_factor*min_dataset
if(is.null(sim.scores)){
  sim.scores =noisify(magnitude = magnitude)
  ccrepe_job=create_ccrepe_jobs(data=refined_table,sim.scores = sim.scores,
                                prefix=prefix,postfix=paste0(postfix,".csv"))
}
if(!is.null(subset))
{
  ccrepe_job=ccrepe_job[subset]
}
stringlist=lapply(ccrepe_job, function(x) list(string=x$string))
ccrepe_res=ccrepe_analysis(ccrepe_job,parallel = parallel)
outputargs=add_outputargs(ccrepe_res,OTU_table=OTU_table,file=file,
                          threshold.value=q_crit,
                          return.value =TRUE)
similarity_measures_significance=lapply(outputargs,
       function(x)(do.call(what=micInt::output_ccrepe_data,
                          args=x))
)
if(is.null(returnVariables)){
  return(mget(ls()))
}
else{
  return(mget(returnVariables))
}
}
#' @title significanceDiganostics
#'
#' @description
#' Makes interaction density plots and abundance
#' product plots for significant interations
#'
#' @param similarity_measures_significance A list over interaction tables, see \link{output_ccrepe_data}.
#'
#' @param refined_table The refined table returned from \link{}
#'
#' @inheritParams runAnalysis
#'
#' @inheritParams output_ccrepe_data
#'
#'
#'
#'
#' @importFrom graphics par plot
significanceDiganostics=function(similarity_measures_significance,OTU_table,refined_table,metadataCols = NULL,
                                 threshold.type='q'){
if(is.null(metadataCols)){
metadataCols = which(unlist(lapply(OTU_table, FUN = function(x) !is.numeric(x)))
)
}
OTU_stat=OTU_stats(OTU_table,metadataCols = metadataCols)
# We want to plot the results, but want to exclude the ones being filtered
# out because of low abundance
OTU_stat=OTU_stat[OTU_stat$ID %in% names(refined_table),]
par(ask=TRUE)
# Plotting number of significant interactions versus the relative abundaces of the
# OTUs
for(i in 1:length(similarity_measures_significance)){
plot(similarity_measures_significance[[i]],OTU_stat,type = 'num_int',abundance_type ='mean',main=names(similarity_measures_significance)[i]
     )
  plot(similarity_measures_significance[[i]],OTU_stat,type = 'num_int',abundance_type ='median',main=names(similarity_measures_significance)[i]
  )
  plot(similarity_measures_significance[[i]],OTU_stat,type = 'num_int',abundance_type ='max',main=names(similarity_measures_significance)[i]
  )
}
# Plotting the q-values/p-value of the significant interactions against the product
# of the mean abundances of the interacting OTUs
  for(i in 1:length(similarity_measures_significance)){
    if(length(similarity_measures_significance[[i]][[valueColumn]])==0){
      next
    }
  plot(similarity_measures_significance[[i]], OTU_stat,type='ab_prod',cutoff_type = type,
  main=names(similarity_measures_significance)[i]
  )
}
}


#' @title autoplot.interaction_table
#'
#' @param table An \code{interaction_table}
#'
#' @param OTU_stat OTU statistics obtained from the function \link{OTU_stats}
#'
#' @param type One of the following: \itemize{
#' \item \code{'num_int'} Plots the number of interactions for each OTU against
#' its abundance
#' \item \code{'ab_prod'} Plots the product of the abundances of each significant pairs of
#' OTUs against the q- or p-value of the interaction
#' }
#'
#' @param cutoff_type
#' One of the following:
#' \itemize{
#' \item \code{'q'} Plots according to the q-values
#' \item \code{'p'} Plots according to the p-values
#' }
#' This parameter is ignored if \code{type='num_int'}
#'
#' @param abundance_type The type of abundance to use in the plots. One of \code{c('mean','median','max')}.
#'
#'
#'
#' @description Makes a diagnostic plot of an \code{interaction_table}
#'
#' @import ggplot2
#'
#' @return A \link{ggplot} object showing the desired diagnostic plot
#'
#'
#' @export
autoplot.interaction_table = function(table,OTU_stat,type = 'num_int',cutoff_type='q',
                                  abundance_type ='mean'){
 plot_after = switch (abundance_type,
   'mean' = OTU_stat$meanAbundance,
   'max' = OTU_stat$maxAbundance,
   'median' = OTU_stat$medianAbundance
 )
if(type == 'num_int'){
num_int= countInteractions(OTU_stat$ID,table)
xlab= switch (abundance_type,
              'mean' = 'Mean abundance',
              'max' = 'Max abundance',
              'median' = 'Median abundance'
)
return(ggplot()+geom_point(aes_string(x="plot_after",y="num_int"))+scale_x_continuous(name= xlab,trans = 'log10')+
  scale_y_continuous(name='Number of interactions'))
}
else{
  if(cutoff_type=='q')
  {
    description='q-value'
    valueColumn='q.value'
  }
  else{
    description='p-value'
    valueColumn='p.value'
  }
abundance_product = abundanceProduct(table,OTU_stat,type=abundance_type)
y = abundance_product
x = table[[valueColumn]]
return(ggplot()+geom_point(aes_string(x="x",y="y"))+xlab(description)+ylab('Abundance product')+
  scale_x_continuous(trans='log10')+scale_y_continuous(trans='log10'))
}
}

#' @title ratio_shared_interactions
#'
#' @description
#' Finds the Jaccard index of number of shared interactions between the similarity measures
#'
#' @param similarity_measures_significance A list objects of the class \code{interaction_table}
#'
#' @return A matrix showing the Jaccard indecies of the number of significant interactions.
#' Measures with no significant interactions are ignored.
#'
#' @export
ratio_shared_interactions=function(similarity_measures_significance){
  # Ignores tables with zero columns
  to_include=unlist(lapply(similarity_measures_significance,function(x) nrow(x)!= 0))
  OTU_pairs=lapply(similarity_measures_significance[to_include],extract_OTUs)
  n_sim_measures=length(OTU_pairs)
  res=matrix(data=NA, nrow=n_sim_measures,ncol=n_sim_measures)
  dimnames(res)=list(names(OTU_pairs),names(OTU_pairs))
  for (i in 1:(n_sim_measures-1)){
    for (j in (i+1):n_sim_measures){
      res[i,j]=count_matches(OTU_pairs[[i]],OTU_pairs[[j]])/min(length(OTU_pairs[[i]]),length(OTU_pairs[[j]]))
      res[j,i]=res[i,j]
    }
  }
  return(res)
}
# Ensures that each pair of OTUs are ordered by the lowest index first
extract_OTUs=function(sim_result){
  OTUs=sim_result[c('OTU_1','OTU_2')]
  OTUs_orded=as.data.frame(apply(OTUs,1,function(x) c(min(x),max(x))))
  OTUs_list=as.list(OTUs_orded)
  return(OTUs_list)
}
# Count the number of identical (not necessarly contigious) rows
count_matches=function(OTUs_1,OTUs_2){
 length(intersect(OTUs_1,OTUs_2))
}

