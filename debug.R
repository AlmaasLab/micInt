#OTU_time_series(subsetted_phyloseq,'Week')
testset <- readRDS('tests/testthat/toy_phyloseq.rds')
t <- runAnalysis(testset,parallel = FALSE,subset = c('pearson','spearman'))
t$similarity_measures_significance$spearman
