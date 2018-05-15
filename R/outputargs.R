add_outputargs=function(ccrepe_res,OTU_table=NULL,file=TRUE,return.value=FALSE,
                        threshold.type='q',threshold.value=0.05,
                        csv_option='2',removeDuplicates=TRUE,stringlist=NULL){
ccrepe_res=lapply(ccrepe_res,function(i){
  n=names(i)
  if ('res' %in% n){
  names(n)=n
  n['res']='data'
  names(i)<-n
  return(i)
  }
}
)
output_defaultargs=list(OTU_table = OTU_table,output.file = file,
                        return.value = return.value,threshold.value=threshold.value,
                        threshold.type=threshold.type,
                        csv_option='2',removeDuplicates=removeDuplicates)
outputargs = lapply(ccrepe_res,function(x) modifyList(x,output_defaultargs))
}
