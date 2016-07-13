context("Algorithm runs properly")

data(toyGenotypes)
data(toyResults)

test_that("algorithm runs on toy data", {
    expect_error(run_stego(NA),"genotypes must be matrix-like object")
    expect_error(run_stego(rnorm(100)),"genotypes must be matrix-like object")
    expect_error(run_stego(matrix(rnorm(1000),ncol=10)),"Non-binary values for phased data")
    expect_error(run_stego(matrix(rbinom(10000,1,.5),ncol=50), phased=F),"Unphased data needs")
    capture.output(res <- run_stego(toyGenotypes, simFun=cor))
    
    expect_equal(res$s_matrix_dip, toyResults$toyAllTogether$s_matrix_dip)
    expect_equal(res$s_matrix_hap, toyResults$toyAllTogether$s_matrix_hap)
    expect_equal(res$pkweightsMean, toyResults$toyAllTogether$pkweightsMean)
    expect_equal(res$var_s_dip, toyResults$toyAllTogether$var_s_dip)
    expect_equal(res$var_s_hap, toyResults$toyAllTogether$var_s_hap)
    expect_equal(res$simMat, toyResults$toyAllTogether$simMat)
    expect_equal(res$kinships, toyResults$toyAllTogether$kinships)
    expect_equal(res$structurePValue, toyResults$toyAllTogether$structurePValue)
    expect_equal(res$crypticPValue, toyResults$toyAllTogether$crypticPValue)
    expect_equal(res$analysisType, toyResults$toyAllTogether$analysisType)
    
    
    
    capture.output(res <- run_stego(toyGenotypes[,rep(c(T,F),100),,with=F]+toyGenotypes[,rep(c(F,T),100),with=F], groups="all.together", phased=F, simFun=cor))
    expect_equal(res, toyResults$toyUnphasedAllTogether)
    
    labels <- paste("Group",c(LETTERS[rep(1:5,20)]))
    expect_warning(capture.output(res <- run_stego(toyGenotypes, groups="each.separately", labels=labels, simFun=cor)))
    expect_equal(res, toyResults$toyEachSeparately)
    
    labels <- paste("Group",c(LETTERS[rep(1:5,10)],LETTERS[rep(6:10,10)]))
    super <- c(rep("Super A",50), rep("Super B",50))
    expect_warning(capture.output(res <- run_stego(toyGenotypes, groups="pairwise.within.superpop", labels=labels, super=super, simFun=cor)))

    expect_equal(res, toyResults$toyPairwiseWithinSuperpop)
    
})

test_that("algorithm throws error on bad data", {
    labels <- paste("Group",c(LETTERS[rep("A",40)]))
    expect_error(capture.output(res <- run_stego(toyGenotypes, groups="each.separately", labels=labels)))
    expect_error(capture.output(res <- run_stego(NA, groups="each.separately")))
})

test_that("data is pruned", {
    expect_equal(pruneGenotypes(toyGenotypes, 1),toyGenotypes)
    expect_equal(pruneGenotypes(toyGenotypes, NA),toyGenotypes)
    expect_equal(dim(pruneGenotypes(toyGenotypes, 2)),c(2500,200))
    expect_equal(dim(pruneGenotypes(toyGenotypes, 100)),c(50,200))
    expect_error(pruneGenotypes(toyGenotypes, 0),"Block size must be a positive integer less than the number of rows in the genotype data.")
    expect_error(pruneGenotypes(toyGenotypes, 10000),"Block size must be a positive integer less than the number of rows in the genotype data.")
    expect_error(pruneGenotypes(toyGenotypes, 2.2),"Block size must be a positive integer less than the number of rows in the genotype data.")
})
