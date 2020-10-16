library(testthat)
library(slingshot)
library(SingleCellExperiment)

data(list = 'slingshotExample', package = "slingshot")
if (!"cl" %in% ls()) {
  rd <- slingshotExample$rd
  cl <- slingshotExample$cl  
}
condition <- factor(rep(c('A','B'), length.out = nrow(rd)))
condition[110:139] <- 'A'
sds <- slingshot(rd, cl)

test_that("Sds does merge correctly",{
    sds_merged <- merge_sds(sds, sds, 
        mapping = matrix(c(1, 1, 2, 2), nrow = 2, byrow = TRUE))
    expect_equal(ncol(reducedDim(sds_merged)), ncol(reducedDim(sds)))
    expect_equal(nrow(reducedDim(sds_merged)), 2 * nrow(reducedDim(sds)))
    expect_equal(slingClusterLabels(sds_merged)[rownames(reducedDim(sds)), ],
                 slingClusterLabels(sds))
    expect_equal(quantile(slingPseudotime(sds_merged, na = FALSE)[, 1]),
                 quantile(slingPseudotime(sds, na = FALSE)[, 1]))
    expect_equal(quantile(slingPseudotime(sds_merged, na = FALSE)[, 2]),
                 quantile(slingPseudotime(sds, na = FALSE)[, 2]))
    
})