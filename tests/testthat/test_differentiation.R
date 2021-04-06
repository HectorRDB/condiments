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

test_that("The differentiationTest work on expected inputs",{
  # Input SlingshotDataSet
  set.seed(23)
  test <- differentiationTest(cellWeights = sds, conditions = condition)
  expect_is(test, "data.frame")
  expect_equal(dim(test), c(1, 3))
  expect_equal(colnames(test),  c("pair", "statistic", "p.value"))
  set.seed(23)
  test_all <- differentiationTest(cellWeights = sds, conditions = condition,
                                  pairwise = TRUE)
  expect_equal(nrow(test_all), choose(nLineages(sds), 2))
  expect_equal(test, test_all[1,])
  set.seed(23)
  test_pairs <- differentiationTest(cellWeights = sds, conditions = condition,
                                        pairwise = TRUE, global = FALSE)
  test_pairs <- as.data.frame(test_pairs)
  expect_equal(nrow(test_pairs), choose(nLineages(sds), 2))
  expect_equivalent(test_pairs, test_all)
  # Input SingleCellExperiment
  pd <- DataFrame(cond = condition)
  rownames(pd) <- rownames(sds)
  sce <- SingleCellExperiment(
    assay = list(counts = t(slingReducedDim(sds))),
    colData = pd
  )
  sce@int_metadata$slingshot <- sds
  set.seed(23)
  test_sce <- differentiationTest(cellWeights = sce, conditions = "cond")
  expect_identical(test_sce, test)
  set.seed(23)
  test_mat <- differentiationTest(cellWeights = slingCurveWeights(sds),
                                  conditions = condition)
  expect_identical(test_mat, test)
})

test_that("The differentiationTest work on all tests",{
  # Input SlingshotDataSet
  set.seed(23)
  test <- differentiationTest(cellWeights = sds, conditions = condition,
                              method = "Classifier")
  expect_is(test, "data.frame")
  expect_equal(dim(test), c(1, 3))
  expect_equal(colnames(test),  c("pair", "statistic", "p.value"))
  test <- differentiationTest(cellWeights = sds, conditions = condition,
                              method = "mmd")
  expect_is(test, "data.frame")
  expect_equal(dim(test), c(1, 3))
  expect_equal(colnames(test),  c("pair", "statistic", "p.value"))
  test <- differentiationTest(cellWeights = sds, conditions = condition,
                              method = "wasserstein_permutation")
  expect_is(test, "data.frame")
  expect_equal(dim(test), c(1, 3))
  expect_equal(colnames(test),  c("pair", "statistic", "p.value"))
})


test_that("The differentiationTest error when it should", {
  pd <- DataFrame(cond = condition)
  rownames(pd) <- rownames(sds)
  sce <- SingleCellExperiment(
    assay = list(counts = t(slingReducedDim(sds))),
    colData = pd
  )
  expect_error(differentiationTest(cellWeights = sce, conditions = "cond"))
})

test_that("The function wroks when reassign is false",{
  sds <- slingshot(rd, cl, reassign = FALSE, reweight = FALSE)
  set.seed(23)
  test <- differentiationTest(cellWeights = sds, conditions = condition,
                              method = "Classifier")
  expect_is(test, "data.frame")
  expect_equal(dim(test), c(1, 3))
  expect_equal(colnames(test),  c("pair", "statistic", "p.value"))
})

