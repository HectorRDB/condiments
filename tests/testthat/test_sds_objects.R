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
    sds_merged <- merge_sds(sds, sds, scale = FALSE,
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

test_that("sds is fitted ok", {
  sdss <- slingshot_conditions(sds, condition)
  expect_is(sdss, "list")
  expect_is(sdss[[1]], "SlingshotDataSet")
  expect_is(sdss[[2]], "SlingshotDataSet")
  sdss <- slingshot_conditions(sds, rep(1, nrow(rd)))
  expect_equal(sdss[[1]], sds)
  # Input SingleCellExperiment
  pd <- DataFrame(cond = condition)
  rownames(pd) <- colnames(sds)
  sce <- SingleCellExperiment(assay = list(counts = t(reducedDim(sds))),
                              colData = pd)
  sce@int_metadata$slingshot <- sds
  set.seed(23)
  sdss1 <- slingshot_conditions(sds, condition)
  sdss2 <- slingshot_conditions(sce, condition)
  expect_equal(sdss1, sdss2)
  sdss2 <- slingshot_conditions(sce, "cond")
  expect_equal(sdss1, sdss2)
})
