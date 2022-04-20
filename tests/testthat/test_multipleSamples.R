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

test_that("All tests work multiple samples", {
  set.seed(22)
  Samples <- sample(1:2, 140, replace = TRUE)
  test <- progressionTest_multipleSamples(pseudotime = slingPseudotime(sds),
                                          cellWeights = slingCurveWeights(sds),
                                          conditions = condition, Samples = Samples)
  expect_is(test, "data.frame")
  test <- fateSelectionTest_multipleSamples(cellWeights = slingCurveWeights(sds),
                                            conditions = condition, Samples = Samples)
  expect_is(test, "data.frame")
  test <- topologyTest_multipleSamples(sds = sds, conditions = condition,
                                       Samples = Samples, rep = 10)
  expect_is(test, "data.frame")
})


test_that("All tests work with expected inputs with multiple samples",{
  set.seed(23)
  Samples <- sample(1:2, 140, replace = TRUE)
  # Input SlingshotDataSet
  set.seed(23)
  test <- progressionTest_multipleSamples(sds, conditions = condition, Samples = Samples)
  expect_is(test, "data.frame")
  test <- fateSelectionTest_multipleSamples(sds, conditions = condition, Samples = Samples)
  expect_is(test, "data.frame")
  test <- topologyTest_multipleSamples(sds, conditions = condition,
                                       Samples = Samples, rep = 10)
  expect_is(test, "data.frame")
  # Input SingleCellExperiment
  pd <- DataFrame(cond = condition)
  rownames(pd) <- rownames(sds)
  sce <- SingleCellExperiment(
    assay = list(counts = t(slingReducedDim(sds))),
    colData = pd
  )
  sce@int_metadata$slingshot <- sds
  set.seed(23)
  test <- progressionTest_multipleSamples(sce, conditions = condition, Samples = Samples)
  expect_is(test, "data.frame")
  test <- fateSelectionTest_multipleSamples(sce, conditions = condition, Samples = Samples)
  expect_is(test, "data.frame")
  test <- topologyTest_multipleSamples(sce, conditions = condition,
                                       Samples = Samples, rep = 10)
  expect_is(test, "data.frame")
})
