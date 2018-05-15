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

#*************************************************************
#' @title similarity_measures
#'
#' @description
#' Similarity measures passed as functions to ccrepe and
#' output_ccrepe_data
#'
#' @param subset The subset of the similarity measures to be returned. \
#' By default, all similarity measures are returned
#'
#' @return A list of the similarity measures with the function ifself (\code{FUN}), and its human reable name (\code{string})
#'
#' @import vegan
#' @import ccrepe
#' @import infotheo
#' @importFrom stats cor
#' @export
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
                                                nc.score(x,y)
}
measures$nc.score$string='nc_score'
#*******************************************************************************
# Euclidian distance
#******************************************************************************
euclidean_similarity=function(x,y=NULL){
  if (is.null(y)){
    ones=rep(1,dim(x)[2])
    res=ones%*%t(ones)-as.matrix(designdist(t(x),method='(A+B-2*J)/(1+A+B-2*J)',terms='quadratic'))
  }
  else{
    res=euclidean_similarity(cbind(x,y))[1,2]
  }
  return(res)
}
measures$euclidean$FUN=euclidean_similarity
measures$euclidean$string='squared_euclidean_similarity'

#*******************************************************************************
# Cosine distance
#******************************************************************************
cosine_similarity=function(x,y=NULL){
  if (is.null(y)){
    ones=rep(1,dim(x)[2])
    res=as.matrix(designdist(t(x),method='J/sqrt(A*B)',terms='quadratic'))
  }
  else{
    res=cosine_similarity(cbind(x,y))[1,2]
  }
  return(res)
}
measures$cosine$FUN=cosine_similarity
measures$cosine$string='cosine_similarity'

score_names=lapply(measures,function(x)x$string)
if(!is.null(subset)){
  measures=measures[score_names %in% subset]
}
return(measures)
}

#' @title create_ccrepe_jobs
#'
#' Creates ccrepe jobs out of sim.score functions
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
#' @title Noisify
#'
#' @description
#' Creates similarity measure functions which adds random noise
#' at certain magnitude to the data before calculating
#' the similarity score
#'
#' @param sim.scores The list of similarity measures to nosify
#'
#' @param magnitude The magnitude of the noise to add. This is: The standard deviation for normal noise and the
#' radius for the interval for uniform noise
#'
#' @param noise The type of noise to add.
#' \code{'none'} just returns the similarity measure without doing anything.
#' \code{'uniform'} adds uniformly distributed noise
#' \code{'normal'} adds normaly distributed  noise
#' Any combination of the options can be passes as a character vector.
#' In this case, sim.score functions are returned for the selected types of noise.
#'
#' @return A list of \code{sim.score} functions corresponding to the input where noise have been
#' added
#'
#' @importFrom stats rnorm runif
#'
#' @export
noisify=function(sim.scores=similarity_measures(),magnitude=1e-5,noise=c('none','uniform','normal')){
  res=list()
  if('none' %in% noise){
    res=c(res,sim.scores)
  }
  if('uniform' %in% noise){
    noiseFUN=function(n){
      runif(n,min=-magnitude,max=magnitude)
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
