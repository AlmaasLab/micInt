context('Data manipulation utilities')
library(micInt)
test_that("Scaling of abundances in phyloseq objects works properly",{
  testset <- phyloseq::phyloseq(phyloseq::otu_table(matrix(c(1,2,3,4),nrow = 2,ncol = 2,byrow = FALSE,
                                                           dimnames = list(c("1","2"),c("sa1","sa2"))),
                                                    taxa_are_rows = TRUE),
                                phyloseq::sample_data(data.frame(scale=c(1,10),row.names = c("1",'2'))))
  scaled_testset <- micInt::scale_by_column(testset,'scale')
  expect_equal(phyloseq::otu_table(scaled_testset),expected = phyloseq::otu_table(matrix(c(1,2,30,40),nrow = 2,ncol = 2,byrow = FALSE,
                                                                     dimnames = list(c("1","2"),c("sa1","sa2"))),taxa_are_rows = TRUE),check.attributes = FALSE)
}
)
