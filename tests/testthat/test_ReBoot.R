context('ReBoot')
library(micInt)
test_that('The ReBoot pipeline works', {
testdata<-matrix(c(0.29787234, 0.2978723, 0.2553191, 0.1489362,
                   0.17073171, 0.3170732, 0.2682927, 0.2439024,
                   0.09302326, 0.3255814, 0.2558140, 0.3255814,
                   0.32352941, 0.3235294, 0.1470588, 0.2058824,
                   0.17241379, 0.1724138, 0.4137931, 0.2413793,
                   0.29729730, 0.2162162, 0.2702703, 0.2162162,
                   0.22500000, 0.3250000, 0.2000000, 0.2500000,
                   0.12820513, 0.3589744, 0.2307692, 0.2820513,
                   0.20000000, 0.2250000, 0.2250000, 0.3500000,
                   0.10256410, 0.3076923, 0.1794872, 0.4102564
),nrow=10,ncol=4,byrow = TRUE)
dimnames(testdata) <-  list(
  c("Subject.1", "Subject.2","Subject.3","Subject.4","Subject.5","Subject.6","Subject.7","Subject.8","Subject.9","Subject.10"),
  c("bug.1", "bug.2", "bug.3","bug.4")) # column names
expect_error(runAnalysis(t(testdata),parallel = FALSE,subset = c('pearson','spearman')),
             regexp = NA)
}
)
