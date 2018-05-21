#' @name sim.measure
#' @title sim.measure
#' @description
#' A similarity measure structure with the function ifself, and metadata fields
#'
#' @slot FUN The similarity measure function to use, of the form \code{f(x,y)}.
#' The function must satisfy the following: \itemize{
#' \item If \code{x} and \code{y} are both given, they are treated as vectors and the
#' similarity score between them is returned
#' \item If only \code{x} is given, it is treated as a matrix (or data frame) where the
#' features are in columns and sample in the rows. The function must then return
#' the matrix of pairwise similarities
#' \item The function must return a similarity score in the interval \eqn{[0,1]}, where 1 is
#' perfect similarity and 0 is perfect dissimilarity. Alternativly, signed similarity scores
#' in the interval \eqn{[-1,1]} can also be used.
#' }
#' @slot string The similarity score's human readable name.
#'
#'
#' @export
setClass(Class="sim.measure",slots=c(FUN="function",string="character"))

#' @name sim.measure
#'@rdname sim.measure
sim.measure=function(FUN,string){
 new("sim.measure",FUN=FUN,string=string)
}


#' @title similarity_measures
#'
#' @description
#' Similarity measures passed as functions to ccrepe and
#' output_ccrepe_data
#'
#' @param subset The subset of the similarity measures to be returned. \
#' By default, all similarity measures are returned
#'
#' @return A list of the \link{sim.measure} objects defined by the function
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
measures$pearson=sim.measure(FUN=pearson_cor,string='pearson')
#************************************************************
# Spearman correlation
#***********************************************************
spearman_cor=function(x,y=NULL){cor(x=x,y=y,method='spearman')}
measures$spearman=sim.measure(FUN=spearman_cor,string='spearman')
#************************************************************
# Kendall's tau
#***********************************************************
kendall_cor=function(x,y=NULL){cor(x=x,y=y,method='kendall')}
measures$kendall=sim.measure(FUN=kendall_cor,string='kendall')
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
measures$bray_curtis=sim.measure(FUN=bray_curtis_cor,string='bray_curtis')
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
measures$jaccard=sim.measure(FUN=jaccard_cor,string='jaccard_index')
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
measures$gen_jaccard(FUN=gen_jaccard_cor,string='generalized_jaccard_index')
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
measures$mutual_information=sim.measure(FUN=mutual_information,string='mutual_information')
#*******************************************************************************
# nc.score
#******************************************************************************
measures$nc.score=sim.measure(FUN=function(x,y=NULL){
  nc.score(x,y)
},string='nc_score')
#*******************************************************************************
# Euclidean distance
#******************************************************************************
euclidean_similarity=function(x,y=NULL){
  if (is.null(y)){
    ones=rep(1,dim(x)[2])
    res=ones%*%t(ones)-as.matrix(designdist(t(x),method='((A+B-2*J)/P)/(1+(A+B-2*J)/P)',terms='quadratic'))
  }
  else{
    res=euclidean_similarity(cbind(x,y))[1,2]
  }
  return(res)
}
measures$euclidean=sim.measure(FUN=euclidean_similarity,string='squared_euclidean_similarity')
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
measures$cosine=sim.measure(FUN=cosine_similarity,string='cosine_similarity')

score_names=lapply(measures,function(x)x@string)
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
                                                            sim.score=sim.score@FUN),
                                              output_args=list(filename =paste0(prefix,'_',sim.score@string,postfix))
                                              ,string=sim.score@string
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
#' @param sim.scores The list of objects of class \link{sim.measure} to nosify
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
#' @return A list of \link{sim.measure} objects corresponding to the input where noise have been
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
                               FUN=noisificationTemplate(sim.score@FUN,noiseFUN),
                               string=paste0(sim.score@string,'_uniform')
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
                                        FUN=noisificationTemplate(sim.score@FUN,noiseFUN),
                                        string=paste0(sim.score@string,'_normal')
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
