source('ccrepe_analysis.R')
source('data_refinement.R')
source('OTU_stats.R')
source('outputargs.R')
source('similarity_measures.R')
runAnalysis=function(OTU_table,abundance_cutoff=1e-04,q_crit=0.05,parallel=TRUE,
                    returnVariables=NULL,subset=NULL){
library(stringr)
prefix=paste('q_crit=',format(q_crit,scientific = TRUE),'_cutoff=',format(abundance_cutoff,
                                                          scientific = TRUE),sep = '')
refined_table=refine_data(OTU_table,abundance_cutoff=abundance_cutoff)
# The smallest value in the data set
min_dataset=min(apply(refined_table,MARGIN = 2,function(x) min(x[x>0])))
magnitude=10*min_dataset
ccrepe_job=create_ccrepe_jobs(data=refined_table,sim.scores =noisify(magnitude = magnitude),
                               prefix=prefix)
if(!is.null(subset))
{ccrepe_job=ccrepe_job[subset]

  }
stringlist=lapply(ccrepe_job, function(x) list(string=x$string))
# Removes the nc.score jobs as the tend to be unreliable on Cruncher
ccrepe_job=ccrepe_job[!str_detect(names(ccrepe_job),'nc.score')]
ccrepe_res=ccrepe_analysis(ccrepe_job,parallel = parallel)
outputargs=add_outputargs(ccrepe_res,OTU_table=OTU_table,file=FALSE,
                          threshold.value=q_crit,
                          return.value =TRUE)
similarity_measures_significance=lapply(outputargs,
       function(x)do.call(output_ccrepe_data,
                          x))
if(is.null(returnVariables)){
save(OTU_table,abundance_cutoff,magnitude,min_dataset,q_crit,
            similarity_measures_significance,refined_table,
     file=paste(prefix,'.RData',sep='')
)
}
else{
  return(mget(returnVariables))
}
}
significanceDiganostics=function(similarity_measures_significance,refined_table,OTU_table,
                                 type='q',score.name='pearson'){
OTU_stat=OTU_stats(OTU_table)
# We want to plot the results, but want to exclude the ones being filtered
# out because of low abundance
OTU_stat=OTU_stat[OTU_stat$ID %in% names(refined_table),]
numInteractions=lapply(similarity_measures_significance,
                       function(x) countInteractions(OTU_stat$ID,x)
                       )
par(ask=TRUE)
# Plotting number of significant interactions versus the relative abundaces of the
# OTUs
for(i in 1:length(numInteractions)){
plot(OTU_stat$meanAbundance,numInteractions[[i]],xlab = 'Mean abundance',
     ylab='Number of interactions',main=names(numInteractions)[i],
     log='x')
}
for(i in 1:length(numInteractions)){
  plot(OTU_stat$medianAbundance,numInteractions[[i]],xlab = 'Median abundance',
       ylab='Number of interactions',main=names(numInteractions)[i],
       log='x')
}
for(i in 1:length(numInteractions)){
  plot(OTU_stat$maxAbundance,numInteractions[[i]],xlab = 'Max abundance',
       ylab='Number of interactions',main=names(numInteractions)[i],
       log='x')
}
# Plotting the q-values/p-value of the significant interactions against the product
# of the mean abundances of the interacting OTUs
if(type=='q')
{
  description='q-value'
  valueColumn='q.value'
}
else{
  description='p-value'
  valueColumn='p.value'
}
abundance_product=lapply(similarity_measures_significance,
                         function(X) abundanceProduct(X,OTU_stat))
  for(i in 1:length(numInteractions)){
    if(length(similarity_measures_significance[[i]][[valueColumn]])==0){
      next
    }
  plot(similarity_measures_significance[[i]][[valueColumn]],abundance_product[[i]],
  main=names(numInteractions)[i],log='xy',xlab=description,ylab='Abundance product'
  )
}
}
# Finds the Jaccard index of number of shared interactions between the similarity measures
ratio_shared_interactions=function(similarity_measures_significance){
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
# Count the number identical (not necessarly contigious) rows
count_matches=function(OTUs_1,OTUs_2){
 length(intersect(OTUs_1,OTUs_2))
}

