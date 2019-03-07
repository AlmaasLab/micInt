context('Lotka-Volterra')
library(micInt)
library(magrittr)
library(phyloseq)
new_phyloseq_object <- readRDS('toy_phyloseq.rds')
reactors = phyloseq::sample_data(new_phyloseq_object)[['Reactor']] %>% unique
test_that('Construction of OTU time series from phyloseq objects works',
          {
time_series <-  lapply(reactors , function(reactor){
  subsetted_phyloseq <-  phyloseq::prune_samples(x=new_phyloseq_object,
                                               samples=phyloseq::sample_data(new_phyloseq_object)[['Reactor']]==reactor)
  OTU_time_series(subsetted_phyloseq, 'Week')
})
names(time_series) <-  names(new_phyloseq_object)
assign(x='time_series',value = time_series, pos= 1)
expect_s4_class(time_series[[1]],'OTU_time_series')
}
)
test_that('Fitting procedure for Lotka-Volterra approach works',{
systems <-  micInt::integralSystem(time_series = time_series)
sol <-  micInt::ridge_fit(systems,weights = c(self=0.1,interaction=0.01))
start <- time_series[[1]]@table[1,,drop=TRUE] %>% unlist()
pred <- predict(sol,start = start,times = c(time_series[[1]]@time_points))
expect_s3_class(micInt::plot_trajectory(time_series_list = list(predicted=pred, reference=time_series[[1]]),label = TRUE),
                'gg')
}
)

