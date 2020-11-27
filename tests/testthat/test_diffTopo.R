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

test_that("The diffTopoTest work on expected inputs",{
  # Input SlingshotDataSet
  test <- diffTopoTest(sds = sds, conditions = condition, rep = 2)
  expect_is(test, "list")
  expect_true(test$statistic >= 0)
  expect_true(test$p.value >= 0 & test$p.value <= 1)
  set.seed(12)
  test <- diffTopoTest(sds = sds, conditions = condition, rep = 2, thresh = 0)
  expect_true(test$statistic > 0)
  expect_is(test, "list")
  # Input SingleCellExperiment
  pd <- DataFrame(cond = condition)
  rownames(pd) <- colnames(sds)
  sce <- SingleCellExperiment(assay = list(counts = t(reducedDim(sds))),
                              colData = pd)
  sce@int_metadata$slingshot <- sds
  set.seed(12)
  test_sce <- diffTopoTest(sds = sce, conditions = "cond", rep = 2, thresh = 0)
  expect_identical(test_sce, test)
})

test_that("The diffTopoTest work on expected tests",{
  # Input SlingshotDataSet
  set.seed(21)
  test <- diffTopoTest(sds = sds, conditions = condition, rep = 2, method = "KS_all")
  expect_is(test, "list")
  expect_true(test$statistic >= 0)
  expect_true(test$p.value >= 0 & test$p.value <= 1)
  test <- diffTopoTest(sds = sds, conditions = condition, rep = 2, method = "Classifier")
  expect_is(test, "list")
  expect_true(test$statistic >= 0)
  expect_true(test$p.value >= 0 & test$p.value <= 1)
  test <- diffTopoTest(sds = sds, conditions = condition, rep = 2, method = "mmd")
  expect_is(test, "list")
  expect_true(test$p.value >= 0 & test$p.value <= 1)
})
