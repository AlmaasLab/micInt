#' @import matrixStats
OTU_stats=function(OTU_table){
ID=OTU_table$`OTU Id`
refined_data=as.matrix(remove_metadata(OTU_table))
meanAbundance=colMeans(refined_data)
medianAbundance=colMedians(refined_data)
maxAbundance=colMaxs(refined_data)
numberNonZero=colSums(refined_data != 0)
propotionNonZero=colMeans(refined_data != 0)
taxonomy=OTU_table$taxonomy
res=data.frame(ID,meanAbundance,medianAbundance,maxAbundance,
           numberNonZero,propotionNonZero,taxonomy)
}
# Counts the number of interactions for each OTU
countInteractions=function(IDs,interactions_table){
  dataset_OTUs=c(interactions_table$`OTU_1`,interactions_table$`OTU_2`)
  counts=vapply(IDs,function(ID){
    sum(dataset_OTUs==ID)
  },FUN.VALUE = 1)
  names(counts)=IDs
  return(counts)
}
# Computes the abunance product for at set of interactions
abundanceProduct=function(interactions_table,OTU_stat,type='mean'){
abundances=switch(type,
  'mean' = OTU_stat$meanAbundance,
  'median'= OTU_stat$medianAbundance,
  'max'= OTU_stat$maxAbundance
)
names(abundances)=OTU_stat$ID
abundance_1=abundances[interactions_table$OTU_1]
abundance_2=abundances[interactions_table$OTU_2]
abundance_1*abundance_2
}

#' @title overall_stats
#' @description Compiles a summary of significant interactions for the given similarity measures
#' @param similarity_measures_significance List of \code{interaction_table} objects returned from
#' \link{create_interaction_table}
#' @param outputargs The corresponding outputargs list sendt to \link{create_interaction_table}
#' @ A \code{data.frame} consisting of the following fields: \itemize{
#' \item \code{name} The name of the similarity measure
#' \item \code{type} The type of similarity measure, equal to the one in the \code{type}
#' slot in it \code{sim.measure} object
#' \item \code{num_significant} The number of significant interaction
#'
#' }
overall_stats=function(similarity_measures_significance){
overall=data.frame()

}
