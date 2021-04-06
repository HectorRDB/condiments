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

test_that("The topologyTest work on expected inputs",{
  # Input SlingshotDataSet
  test <- topologyTest(sds = sds, conditions = condition, rep = 2)
  expect_is(test, "data.frame")
  expect_true(test$statistic >= 0)
  expect_true(test$p.value >= 0 & test$p.value <= 1)
  set.seed(12)
  test <- topologyTest(sds = sds, conditions = condition, rep = 2, threshs = 0)
  expect_true(test$statistic > 0)
  expect_is(test, "data.frame")
  # Input SingleCellExperiment
  pd <- DataFrame(cond = condition)
  rownames(pd) <- rownames(sds)
  sce <- SingleCellExperiment(
    assay = list(counts = t(slingReducedDim(sds))),
    colData = pd
  )
  sce@int_metadata$slingshot <- sds
  set.seed(12)
  test_sce <- topologyTest(sds = sce, conditions = "cond", rep = 2, threshs = 0)
  expect_identical(test_sce, test)
})

test_that("The topologyTest work on expected tests",{
  # Input SlingshotDataSet
  set.seed(21)
  test <- topologyTest(sds = sds, conditions = condition, rep = 2, methods =
                         c("KS_all", "Classifier", "mmd", "wasserstein_permutation"),
                       threshs = c(0, .01))
  expect_is(test, "data.frame")
  expect_true(all(test$statistic[1:4] >= 0))
  expect_true(all(test$p.value >= 0) & all(test$p.value <= 1))
})
