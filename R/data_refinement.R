#*********************************************************
# ********************************************************
#' @title remove_metadata
#' @description
#' Removes metadata columns from dataset
remove_metadata=function(OTU_table,metadataCols=c('OTU Id','taxonomy')){
  options(stringsAsFactors = FALSE)
  metadata=which(metadataCols == names(OTU_table))
  # Removes the metadata for the data set
  refined_table=as.data.frame(t(OTU_table[,-metadata]))
  colnames(refined_table)=OTU_table$`OTU Id`
  return(refined_table)
}
#' @title cut_abundances
#'
#' @description
#' Removes metadata columns from dataset
#'
#' @import matrixStats
cut_abundances=function(refined_table,abundance_cutoff=0,type='mean',renormalize=TRUE){
  m_refined_table=as.matrix(refined_table)
  abundances=switch(type,
         mean=colMeans(refined_table),
         median=colMedians(m_refined_table),
         max=colMaxs(m_refined_table),
         numberNonZero=colSums(m_refined_table != 0),
         proprotionNonZero=colMeans(m_refined_table != 0)
         )
  refined_table=refined_table[,abundances>abundance_cutoff]
  if(renormalize){
    # Renormalizes the table in order for the rows to sum to 1
    refined_table=renormalize(refined_table)
  }
  return(refined_table)
}
#' @title refine_data
#' @description
#' Removes metadata from OTU table and cuts off the least abundant
#' species, defined by the cutoff parameter
#' @export
refine_data=function(OTU_table,abundance_cutoff=0,cutoff_type='mean',renormalize=TRUE)
  {
  refined_table=remove_metadata(OTU_table)
  # Cuts away the least abundant species
  cut_abundances(refined_table,abundance_cutoff,type=cutoff_type,renormalize=renormalize)
}
renormalize=function(table)
{
  matrix=as.matrix(table)
  res=diag(1/rowSums(matrix))%*%matrix
  table[,]=res
  return(table)
}
#' @title output_ccrepe_data
#' @description
#' Takes input from ccrepe and transforms it into a convenient table
#'
#' @param data The results from ccrepe
#'
#' @importFrom utils modifyList write.csv write.csv2
#'
#' @export
output_ccrepe_data=function(data,OTU_table=NULL,threshold.type='q',threshold.value=0.05,output.file=FALSE,filename=NULL,
                   return.value=TRUE,csv_option='2',removeDuplicates=TRUE){
                    options(stringsAsFactors = FALSE)
                    p.values=as.data.frame(data$p.values)
                    z.stat=as.data.frame(data$z.stat)
                    sim.score=as.data.frame(data$sim.score)
                    q.values=as.data.frame(data$q.values)
                    if (threshold.type =='q'){
                      threshold_matrix=q.values
                    }
                    else if(threshold.type =='p'){
                      threshold_matrix=p.values
                    }
                    else{
                      stop("no valid theshold method given, must be 'p' (local p-value)
                           or 'q' (familywise false discovery rate)")
                    }
                    significant_pairs=which(threshold_matrix<threshold.value,arr.ind = TRUE)
                    significant_interactions = as.data.frame(matrix(colnames(threshold_matrix)[significant_pairs],ncol=2))
                    colnames(significant_interactions)=c('OTU_1','OTU_2')
                    if(nrow(significant_interactions)==0){
                      significant_interactions$sim.score=numeric()
                      significant_interactions$p.value=numeric()
                      significant_interactions$q.value=numeric()
                      significant_interactions$z.stat=numeric()
                      if (!is.null(OTU_table)){
                        # Add taxonomy if available
                        significant_interactions$taxonomy_1=character()
                        significant_interactions$taxonomy_2=character()
                      }
                    }
                    else{
                    # Removing duplicates
                    if(removeDuplicates){
                    significant_interactions=significant_interactions[significant_interactions$OTU_1>significant_interactions$OTU_2, ]
                    }
                    if(!is.null(sim.score) && nrow(sim.score)!=0){
                    significant_interactions$sim.score=apply(significant_interactions,1,function(x){sim.score[x[1],x[2]]})
                    }
                    if(!is.null(p.values)&& nrow(p.values)!=0){
                    significant_interactions$p.value=apply(significant_interactions,1,function(x){p.values[x[1],x[2]]})
                    }
                    if(!is.null(q.values)&& nrow(q.values)!=0)
                    significant_interactions$q.value=apply(significant_interactions,1,function(x){q.values[x[1],x[2]]})
                    if(!is.null(z.stat)&& nrow(z.stat)!=0){
                    significant_interactions$z.stat=apply(significant_interactions,1,function(x){z.stat[x[1],x[2]]})
                    }
                    if (!is.null(OTU_table)){
                      # Add taxonomy if available
                      significant_interactions$taxonomy_1=OTU_table[match(significant_interactions$OTU_1,as.character(OTU_table$`OTU Id`)),]$taxonomy
                      significant_interactions$taxonomy_2=OTU_table[match(significant_interactions$OTU_2,as.character(OTU_table$`OTU Id`)),]$taxonomy
                    }
                    # Sorts by dersired significance column
                    if (threshold.type=='q')
                      significant_interactions=significant_interactions[order(significant_interactions$q.value),]
                    else
                      significant_interactions=significant_interactions[order(significant_interactions$p.value),]
                    }
                    if(output.file){
                      if (csv_option=='2'){
                      write.csv2(significant_interactions,file =filename)
                        }
                      else{
                        write.csv(significant_interactions,file =filename)
                      }
                    }
                    if(return.value){
                      return(significant_interactions)
                    }
                    else{
                      return(NULL)
                    }
}
