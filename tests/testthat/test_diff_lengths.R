library(testthat)
library(slingshot)
library(SingleCellExperiment)

# 1 and 2 lineages ----
data(list = 'slingshotExample', package = "slingshot")
if (!"cl" %in% ls()) {
  rd <- slingshotExample$rd
  cl <- slingshotExample$cl
}
set.seed(3)
condition <- factor(rep(c('A','B'), length.out = nrow(rd)))
condition[110:139] <- 'A'
condition[cl == 4] <- 'A'
sds <- slingshot(rd, cl)

test_that("All tests work with one or two lineages", {
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
  test_1 <- fateSelectionTest(sds, conditions = condition)
  set.seed(07)
  test_2 <- fateSelectionTest(cellWeights = slingCurveWeights(sds),
                              conditions = condition)
  expect_equal(test_1, test_2)
})
