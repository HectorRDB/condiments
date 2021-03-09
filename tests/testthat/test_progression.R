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


test_that("The progressionTest work on expected inputs",{
  # Input SlingshotDataSet
  set.seed(23)
  test <- progressionTest(pseudotime = sds, conditions = condition, rep = 2)
  expect_is(test, "data.frame")
  expect_equal(dim(test), c(1, 3))
  expect_equal(colnames(test),  c("lineage", "statistic", "p.value"))
  set.seed(23)
  test_all <- progressionTest(pseudotime = sds, conditions = condition,
                              rep = 2, lineage = TRUE)
  expect_equal(nrow(test_all), length(slingCurves(sds)) + 1)
  expect_equal(test[, 2:3], test_all[1, 2:3])
  set.seed(23)
  test_lineages <- progressionTest(pseudotime = sds, conditions = condition,
                                   rep = 2, lineage = TRUE, global = FALSE)
  test_lineages <- as.data.frame(test_lineages)
  rownames(test_lineages) <- 2:3
  expect_equal(nrow(test_lineages), length(slingCurves(sds)))
  expect_equal(test_lineages[1, 2:3], test_all[2, 2:3])
  # Input SingleCellExperiment
  pd <- DataFrame(cond = condition)
  rownames(pd) <- colnames(sds)
  sce <- SingleCellExperiment(assay = list(counts = t(reducedDim(sds))),
                              colData = pd)
  sce@int_metadata$slingshot <- sds
  set.seed(12)
  test_sce <- progressionTest(pseudotime = sce, conditions = "cond", rep = 2)
  expect_identical(test_sce, test)
  # Input matrix
  pst <- slingPseudotime(sds, na = FALSE)
  ws <- slingCurveWeights(sds)
  set.seed(12)
  test_mat <- progressionTest(pseudotime = pst, cellWeights = ws, conditions = condition, rep = 2)
  expect_identical(test_mat, test)
})

test_that("The progressionTest work on all tests", {
  set.seed(23)
  test <- progressionTest(pseudotime = sds, conditions = condition, rep = 2,
                          method = "Classifier")
  expect_is(test, "data.frame")
  expect_equal(dim(test), c(1, 3))
  expect_equal(colnames(test),  c("lineage", "statistic", "p.value"))
  set.seed(23)
  test <- progressionTest(pseudotime = sds, conditions = condition, rep = 2,
                          method = "Permutation")
  expect_is(test, "data.frame")
  expect_equal(dim(test), c(1, 3))
  expect_equal(colnames(test),  c("lineage", "statistic", "p.value"))
  set.seed(23)
  test <- progressionTest(pseudotime = sds, conditions = condition, rep = 2,
                          method = "mmd")
  expect_is(test, "data.frame")
  expect_equal(dim(test), c(1, 3))
  expect_equal(colnames(test),  c("lineage", "statistic", "p.value"))
  set.seed(23)
  condition <- factor(rep(c('A','B'), length.out = nrow(rd)))
  test <- progressionTest(pseudotime = sds, conditions = condition, rep = 2,
                          method = "wasserstein_permutation")
  expect_is(test, "data.frame")
  expect_equal(dim(test), c(1, 3))
  expect_equal(colnames(test),  c("lineage", "statistic", "p.value"))
})

test_that("The progressionTest error when it should", {
  # Error when wrong condition and sds
  expect_error(progressionTest(pseudotime = sds, conditions = condition, rep = 2,
                               method = "bulls**t"))
  # Error when wrong condition and matrix
  pst <- slingPseudotime(sds, na = FALSE)
  ws <- slingCurveWeights(sds)
  set.seed(12)
  expect_error(progressionTest(pseudotime = pst, cellWeights = ws,
                               conditions = condition, rep = 2,
                               method = "bulls**t"))
  # Error when missing sds
  pd <- DataFrame(cond = condition)
  rownames(pd) <- colnames(sds)
  sce <- SingleCellExperiment(assay = list(counts = t(reducedDim(sds))),
                              colData = pd)
  expect_error(progressionTest(pseudotime = sce, conditions = "cond", rep = 2))
})
