library(testthat)
library(slingshot)
library(SingleCellExperiment)

data(list = 'slingshotExample', package = "slingshot")
if (!"cl" %in% ls()) {
  rd <- slingshotExample$rd
  cl <- slingshotExample$cl
}

# Create fake data
set.seed(3)
condition <- factor(rep(c('A','B'), length.out = nrow(rd)))
condition[110:139] <- 'A'
sds <- slingshot(rd, cl)
keep <- slingCurveWeights(sds)[,1] > 0
rd <- rd[keep, ]
cl <- cl[keep]
condition <- condition[keep]
sds <- slingshot(rd, cl)

test_that("All tests work", {
  set.seed(22)
  test <- topologyTest(sds, conditions = condition, rep = 3)
  expect_is(test, "data.frame")
  set.seed(07)
  test <- progressionTest(sds, conditions = condition)
  expect_equal(nrow(test), 1)
  set.seed(07)
  test_full <- progressionTest(sds, conditions = condition, lineages = TRUE)
  expect_equal(test, test_full[1,])
  set.seed(07)
  expect_error(test <- differentiationTest(sds, conditions = condition))
  expect_error(differentiationTest(cellWeights = slingCurveWeights(sds),
                                   conditions = condition))
})
