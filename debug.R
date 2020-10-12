#OTU_time_series(subsetted_phyloseq,'Week')
# testset <- readRDS('tests/testthat/toy_phyloseq.rds')
# t <- runAnalysis(testset,parallel = FALSE,subset = c('pearson','spearman'))
# t$similarity_measures_significance$spearman
library(phyloseq)
data("soilrep")
subdivide_by_environment(soilrep,variables=c('Treatment','warmed','clipped'),keep_variables=FALSE,keep_empty=FALSE)
