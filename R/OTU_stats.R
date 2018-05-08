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
