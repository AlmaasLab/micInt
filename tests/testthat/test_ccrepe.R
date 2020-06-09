library(testthat)
find_expr <- function(name, env = parent.frame()) {
      subs <- do.call("substitute", list(as.name(name), env))
        paste0(deparse(subs, width.cutoff = 500), collapse = "\n")
  }

is_approximately <- function(expected,tol=10e-7,label=NULL)
{
      if (is.null(label)) {
              label <- find_expr("expected")
          } else if (!is.character(label) || length(label) != 1) {
                  label <- deparse(label)
              }
      
    function(actual)
        {
            same <- all.equal.numeric(as.vector(actual),as.vector(expected),tol=tol)

            expectation(
                 "success",
                identical(same,TRUE),
                paste0("not equal to ", label, " within tolerance ",tol,"\n", same)
                )
        }
}

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
		

                                        #***********************************************************
                                        #* ccrepe tests against testdata                           *
                                        #***********************************************************
context("CCREPE")
ccrepe.results <-  ccrepe  (x=testdata,   min.subj=10,iterations = 10000)
tol = 0.05

p.values.results <- as.matrix(read.table('../p.values.results.txt'))
     
q.values.results <- as.matrix(read.table('../q.values.results.txt'))

sim.score.results <- as.matrix(read.table('../sim.score.results.txt'))
							
z.stat.results <- as.matrix(read.table('../z.stat.results.txt'))


expect_equal(p.values.results,ccrepe.results$p.values,tolerance=tol)
expect_equal(q.values.results,ccrepe.results$q.values,tolerance=tol)
expect_equal(sim.score.results,ccrepe.results$sim.score, tolerance=tol)
expect_equal(z.stat.results,ccrepe.results$z.stat, tolerance=tol)

context("Absolute abundances")

testdata_absolute <- diag(exp(seq_len(10))) %*% testdata
rownames(testdata_absolute) <- rownames(testdata)
ccrepe.results_absolute <-  ccrepe  (x=testdata_absolute,   min.subj=10,renormalize = FALSE,iterations = 10000)
p.values.results_absolute <- as.matrix(read.table('../p.values.results_absolute.txt'))

q.values.results_absolute <- as.matrix(read.table('../q.values.results_absolute.txt'))

sim.score.results_absolute <- as.matrix(read.table('../sim.score.results_absolute.txt'))

z.stat.results_absolute <- as.matrix(read.table('../z.stat.results_absolute.txt'))
expect_equal(p.values.results_absolute,ccrepe.results_absolute$p.values, tolerance=tol)
expect_equal(q.values.results_absolute,ccrepe.results_absolute$q.values,tolerance=tol)
expect_equal(sim.score.results_absolute,ccrepe.results_absolute$sim.score, tolerance=tol)
expect_equal(z.stat.results_absolute,ccrepe.results_absolute$z.stat, tolerance=tol)
testdata_absolute_wrong <- abs(testdata_absolute + matrix(rnorm(n=10*4,mean=0, sd=5e+5),nrow = 10, ncol = 4))
ccrepe.results_absolute_wrong <- ccrepe  (x=testdata_absolute_wrong,   min.subj=10,renormalize = FALSE,iterations = 10000)

# Should result in failures
expect_false(isTRUE(all.equal(p.values.results_absolute,ccrepe.results_absolute_wrong$p.values,tolerance=tol)))
expect_false(isTRUE(all.equal(q.values.results_absolute,ccrepe.results_absolute_wrong$q.values,tolerance=tol)))
expect_false(isTRUE(all.equal(sim.score.results_absolute,ccrepe.results_absolute_wrong$sim.score,tolerance=tol)))
expect_false(isTRUE(all.equal(z.stat.results_absolute,ccrepe.results_absolute_wrong$sim.score,tolerance=tol)))

context('Difference between absolute and relative data')
expect_false(isTRUE(all.equal(sim.score.results_absolute,sim.score.results,tolerance=tol)))

context("Memory optimization")
ccrepe.results.memory.optimized <-  ccrepe  (x=testdata,   min.subj=10,iterations = 10000,memory.optimize = TRUE)
expect_equal(p.values.results,ccrepe.results.memory.optimized$p.values,tolerance=tol)
expect_equal(q.values.results,ccrepe.results.memory.optimized$q.values,tolerance=tol)
expect_equal(sim.score.results,ccrepe.results.memory.optimized$sim.score, tolerance=tol)
expect_equal(z.stat.results,ccrepe.results.memory.optimized$z.stat, tolerance=tol)


# context("Permutation p-values")

## test_that("Renormalization gives normalized data without missing", {
##     data <- matrix(1:10,nrow=2,byrow=TRUE)
##     data.norm <- ccrepe_norm(data)

##     expect_that( is.matrix(data.norm), is_true())
##     expect_that( apply(data.norm,1,sum), equals(c(1,1)) )
##     expect_that( data.norm, equals( matrix(c(1/15,2/15,3/15,4/15,5/15,6/40,7/40,8/40,9/40,0.25),
##                                            byrow=TRUE,
##                                            nrow=2) ) )
## })

## test_that("Renormalization works with missing data", {
##     data <- matrix(1:10,nrow=2,byrow=TRUE)
##     data[1,2] = NA
##     data.norm <- ccrepe_norm(data)

##     expect_that( data.norm[1,], equals(rep(0,ncol(data))) )
## })


## test_that("Permutation of one matrix permutes the data by columns", {
##     data <- matrix(1:12, ncol=3)
##     permute.id.matrix <- matrix(c(4,3,2,1,2,1,4,3,2,3,1,4),ncol=3)
##     data.permute <- permute(data,permute.id.matrix)

##     expect_that( is.matrix(data.permute), is_true() )
##     expect_that( data.permute, equals(matrix(c(4,3,2,1,6,5,8,7,10,11,9,12),ncol=3)) )
## })

## v_dist.na   <- c(3.4,2,NA,2,4.6,0,10,2,8,NA,-2,0,NA)
## v_dist      <- c(3.4,2,2,4.6,0,10,2,8,-2,0)
## v_dist.null <- NA
## obs.value1 <- 3.4             # high p-value, in dist
## obs.value2 <- 10              # low p-value, in dist
## obs.value3 <- 4               # high p-value, not in dist
## obs.value4 <- -5              # low p-value, not in dist
## dist.value1 <- 2              # repeated value
## dist.value2 <- 0              # repeated value

## test_that("get_count returns a correct value with no missing", {

##     count1 <- get_count(dist.value1, v_dist)
##     count2 <- get_count(dist.value2, v_dist)
##     count3 <- get_count(obs.value3, v_dist)
##     count4 <- get_count(obs.value1, v_dist)

##     expect_that( count1, equals(3) )
##     expect_that( count2, equals(2) )
##     expect_that( count3, equals(0) )
##     expect_that( count4, equals(1) )
## })

## test_that("get_perm_p.value works with no missing", {

##     p.value1 <- get_perm_p.value(v_dist,obs.value1)
##     p.value2 <- get_perm_p.value(v_dist,obs.value2)
##     p.value3 <- get_perm_p.value(v_dist,obs.value3)
##     p.value4 <- get_perm_p.value(v_dist,obs.value4)

##     expect_that( p.value1, equals(1) )
##     expect_that( p.value2, equals(.1) )
##     expect_that( p.value3, equals(.9) )
##     expect_that( p.value4, equals(0) )
## })

## test_that("get_perm_p.value works with missing", {

##     p.value1 <- get_perm_p.value(v_dist.na,obs.value1)
##     p.value2 <- get_perm_p.value(v_dist.na,obs.value2)
##     p.value3 <- get_perm_p.value(v_dist.na,obs.value3)
##     p.value4 <- get_perm_p.value(v_dist.na,obs.value4)

##     expect_that( p.value1, equals(1) )
##     expect_that( p.value2, equals(.1) )
##     expect_that( p.value3, equals(.9) )
##     expect_that( p.value4, equals(0) )
## })

## test_that("get_perm_p.value returns missing if only missing in v_dist", {

##     p.value <- get_perm_p.value(v_dist.null,obs.value1)

##     expect_that( is.na(p.value), is_true() )
## })

## test_that("get_perm_p.value returns missing if obs.value is missing", {

##     p.value <- get_perm_p.value(v_dist, NA)

##     expect_that( is.na(p.value), is_true() )
## })
