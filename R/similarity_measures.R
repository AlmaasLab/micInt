#*************************************************************
# Outputs the file name specifications for the similarity measures
# This must be updated as the similarity scores are added
#*************************************************************
score_names=function(){
  # The sim score names without noise
  basic=list(pearson='pearson',spearman='spearman',bray_curtis='bray_curtis',
       kendall='kendall',bray_curtis='bray_curtis',jaccard='jaccard_index',
       gen_jaccard='generalized_jaccard_index',nc.score='nc_score',
       mutual_information='mutual_information')
}
#*************************************************************
# Similarity measures passed as functions to ccrepe and
# output_ccrepe_data
#*************************************************************
similarity_measures=function(subset=NULL){
measures=list()
score_names=score_names()
#************************************************************
# Pearson linear correlation
#***********************************************************
pearson_cor=function(x,y=NULL){cor(x=x,y=y,method='pearson')}
measures$pearson$FUN=pearson_cor
measures$pearson$string='pearson'
#************************************************************
# Spearman correlation
#***********************************************************
spearman_cor=function(x,y=NULL){cor(x=x,y=y,method='spearman')}
measures$spearman$FUN=spearman_cor
measures$spearman$string='spearman'
#************************************************************
# Kendall's tau
#***********************************************************
kendall_cor=function(x,y=NULL){cor(x=x,y=y,method='kendall')}
measures$kendall$FUN=kendall_cor
measures$kendall$string='kendall'
#************************************************************
# Bray-Curtis
#*************************************************************
# We must transform the disimilarities into similarities
bray_curtis_cor=function(x,y=NULL){
  library(vegan)
  if (is.null(y)){
    ones=rep(1,dim(x)[2])
    res=ones%*%t(ones)-as.matrix(vegdist(t(x),method='bray',upper=TRUE))
  }
  else{
    # The vegdist-function accepts matrix inputs only, but this can be fixed
    # the following way:
    res=bray_curtis_cor(cbind(x,y))[1,2]
  }
  return(res)
}
measures$bray_curtis$FUN=bray_curtis_cor
measures$bray_curtis$string='bray_curtis'
#************************************************************
# Jaccard index
#*************************************************************
jaccard_cor=function(x,y=NULL){
  library(vegan)
  if (is.null(y)){
    ones=rep(1,dim(x)[2])
    res=ones%*%t(ones)-as.matrix(vegdist(t(x),method='jaccard',binary = TRUE,upper=TRUE))
  }
  else{
    res=jaccard_cor(cbind(x,y))[1,2]
  }
  return(res)
}
measures$jaccard$FUN=jaccard_cor
measures$jaccard$string='jaccard_index'
#************************************************************
# Generalized Jaccard
#*************************************************************
gen_jaccard_cor=function(x,y=NULL){
  library(vegan)
  if (is.null(y)){
    ones=rep(1,dim(x)[2])
    res=ones%*%t(ones)-as.matrix(designdist(t(x),method='1-J/(A+B-J)',terms='minimum'))
  }
else{
  res=gen_jaccard_cor(cbind(x,y))[1,2]
}
return(res)
}
measures$gen_jaccard$FUN=gen_jaccard_cor
measures$gen_jaccard$string='generalized_jaccard_index'
#************************************************************
# Mutual information
#*************************************************************
mutual_information=function(x,y=NULL){
  library(infotheo)
  if(is.null(y)){
    disc=discretize(x)
    res=mutinformation(disc)
  }
  else{
    res=mutual_information(cbind(x,y))[1,2]
  }
  return(res)
}
measures$mutual_information$FUN=mutual_information
measures$mutual_information$string='mutual_information'
#*******************************************************************************
# nc.score
#******************************************************************************
measures$nc.score$FUN=function(x,y=NULL){
                                                library(ccrepe)
                                                nc.score(x,y)
}
measures$nc.score$string='nc_score'
score_names=lapply(measures,function(x)x$string)
if(!is.null(subset)){
  measures=measures[score_names %in% subset]
}
return(measures)
}

#*************************************************************
# Creates ccrepe jobs out of sim.score functions
#*************************************************************
create_ccrepe_jobs=function(data=NULL,sim.scores=similarity_measures(),ccrepe_defaultargs=NULL,prefix='significant_interactions'
                            ,postfix='.csv')
{
  # Note: At least one of data and ccrepe_defaultargs must be different
  # from NULL
  if(is.null(ccrepe_defaultargs))
  {
    ccrepe_defaultargs=list(x=data,min.subj = 10,verbose = TRUE)
  }
  jobs=lapply(sim.scores, function(sim.score) list(ccrepe_args=c(ccrepe_defaultargs,
                                                            sim.score=sim.score$FUN),
                                              output_args=list(filename =paste(prefix,'_',sim.score$string,postfix,sep=''))
                                              ,string=sim.score$string
  )
  )
}
#*************************************************************
# Creates similarity measure functions which adds random noise
# at certain magnitude to the data before calculating
# the similarity score
#*************************************************************
noisify=function(sim.scores=similarity_measures(),magnitude=1e-5,noise=c('none','uniform','normal')){
  res=list()
  if('none' %in% noise){
    res=c(res,sim.scores)
  }
  if('uniform' %in% noise){
    noiseFUN=function(n){
      runif(n,min=-magnitude,max=mangnitude)
    }
    uniform_functions=lapply(sim.scores,
                             function(sim.score) list(
                               FUN=noisificationTemplate(sim.score$FUN,noiseFUN),
                               string=paste0(sim.score$string,'_uniform')
                               )
    )
                               
    names(uniform_functions)=sapply(names(sim.scores),
                                    function(name) paste0(name,'_uniform'))
    res=c(res,uniform_functions)                      
                             }
  if('normal' %in% noise){
    noiseFUN=function(n){
      rnorm(n,mean=0,sd=magnitude)
    }
    normal_functions=lapply(sim.scores,
                                      function(sim.score) list(
                                        FUN=noisificationTemplate(sim.score$FUN,noiseFUN),
                                        string=paste0(sim.score$string,'_normal')
                                      )
    )
    names(normal_functions)=sapply(names(sim.scores),
                                  function(name) paste0(name,'_normal'))
    res=c(res,normal_functions)
  }
  return(res)
  }
  
#*************************************************************
# Helper function for nosify, creates the acutual noised
# functions
#*************************************************************
noisificationTemplate=function(FUN,noiseFUN)
  {
  addNoise=function(x){
    abs(x+noiseFUN(n=length(x)))
  }
    noisifiedFunction=function(x,y=NULL)
      if(is.null(y)){
        noisedFrame=apply(x,MARGIN = 2,
                    FUN = addNoise)
        return(FUN(noisedFrame))
      }
    else{
      noisedX=addNoise(x)
      noisedY=addNoise(y)
      return(FUN(noisedX,noisedY))
    }
}
