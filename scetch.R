temp_1=apply(a,1,names)
temp_2=apply(a,2,names)
library(caret)
createFolds(y=rep(1,100))
data(oil)
createResample(1:100, times = 10, list = TRUE)
createDataPartition(1:100)
groupKFold(1:100,k=10)
createTimeSlices(1:100,10)
t= list(data.frame(list(a=2,b=5)),data.frame(list(a=6,b=-2)))

d=do.call(rbind,t)
a=paste(c('a','b','c'))
length(a)
do.call(paste,list('a','b','c',sep='_'))
paste(x=c('a','b','c'),sep='_',collapse = NULL)
temp=c(a=1,b=2,c=3)
temp['d']
library(magrittr)
filtered=temp[c('a','d','b')] %>% na.omit
class(filtered)
filtered
