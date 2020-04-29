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
ds=data.frame(a=c('a','b'),b=c('c','d'))
rownames(ds)=c('Rad 1','Rad 2')
apply(ds,MARGIN=1,FUN=paste,collapse=' ')
time = OTU_time_series(phyloseq_object,time_points = phyloseq_object@sam_data$Week)
time_2 = OTU_time_series(phyloseq_object,time_points = 'Week')
comp_phyloseq = phyloseq_object %>% microbiome::transform(transform = 'compositional')
new_phyloseq_object <- phyloseq::prune_taxa(taxa = comp_phyloseq %>%
                                              phyloseq::taxa_sums(.)/phyloseq::nsamples(comp_phyloseq) >
                                              10^{-3},comp_phyloseq) %>%  phyloseq::subset_samples(Source=='Water')

reactors = phyloseq::sample_data(new_phyloseq_object)[['Reactor']] %>% unique

time_series = lapply(reactors , function(reactor){
                    subsetted_phyloseq = phyloseq::prune_samples(x=new_phyloseq_object,
                                                                 samples=phyloseq::sample_data(new_phyloseq_object)[['Reactor']]==reactor)
                     OTU_time_series(subsetted_phyloseq, 'Week')
  }

                     )
names(time_series) = reactors
time_series
systems = micInt::integralSystem(time_series = time_series)
sol = micInt::ridge_fit(systems,weights = c(self=0.1,interaction=0.01))
start = time_series[[1]]@table[1,,drop=TRUE] %>% unlist()
pred <- predict(sol,start = start,times = c(time_series[[1]]@time_points))
micInt::plot_trajectory(time_series_list = list(predicted=pred, reference=time_series[[1]]),label = TRUE)
x <- seq_len(100)
y <- sqrt(x)
label <- rep('Hallo!',100)
frame <- data.frame(x,y,label)
library(ggplot2)
ggplot(frame)+aes_string(x="x",y="y")+geom_text(aes_string(label="label"),cex=1)
soilrep
sample_data(soilrep)
table <- list(a=c(1,4),b=c(2,5),c=c(3,7)) %>% as_tibble()
table
liste_2 <- list(a=c(1,4),b=c(4,1),c=c(3,7)) %>% as_tibble()
pmap(liste,sum)
pmap_int(names(liste_2))
debug_pipe()
library(pipecleaner)
debug_pipeline(subdivide_by_environment(soilrep,variables=c('Treatment','warmed','clipped')))
b <- subdivide_by_environment(soilrep,variables=c('Treatment','warmed','clipped'),keep_variables = FALSE,keep_empty = FALSE)
table %>% rowwise() %>% mutate(d=lapply(list(.data), FUN=function(variable){
variable + 1
}) %>% median()
)
table %>% rowwise() %>% do({lapply(.data, FUN=function(variable){
  variable + 1
}) %>% unlist() %>% median() %>% tibble
  })
d = numeric(nrow(liste))
for(i in 1:nrow(liste)){
row <- liste[i,]
d[i]=lapply(row,FUN=function(variable){
  variable + 1
}) %>% unlist() %>% median()
}
s <- c(rep(0,50),rep(100,6))
scaled_sample <- micInt::scale_by_column(soilrep,s)
new <- scale(otu_table(soilrep),center = FALSE,scale = s)
otu_table(scaled_sample)
test_phyloseq = phyloseq()
x <- c(1,1)
y <- c(1,1)
cos_sim <- similarity_measures()[['cosine']]
cos_sim@FUN(x,y)
