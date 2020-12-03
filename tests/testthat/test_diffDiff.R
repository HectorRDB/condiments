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
  test <- differentiationTest(sds = sds, conditions = condition)
  expect_is(test, "data.frame")
  expect_equal(dim(test), c(1, 3))
  expect_equal(colnames(test),  c("pair", "statistic", "p.value"))
  set.seed(23)
  test_all <- differentiationTest(sds = sds, conditions = condition, pair = TRUE)
  expect_equal(nrow(test_all), choose(nLineages(sds), 2) + 1)
  expect_equal(test, test_all[1,])
  set.seed(23)
  test_pairs <- differentiationTest(sds = sds, conditions = condition,
                                        pair = TRUE, global = FALSE)
  test_pairs <- as.data.frame(test_pairs)
  rownames(test_pairs) <- "2"
  expect_equal(nrow(test_pairs), choose(nLineages(sds), 2))
  expect_equivalent(test_pairs[1, ], test_all[2, ])
  # Input SingleCellExperiment
  pd <- DataFrame(cond = condition)
  rownames(pd) <- colnames(sds)
  sce <- SingleCellExperiment(assay = list(counts = t(reducedDim(sds))),
                              colData = pd)
  sce@int_metadata$slingshot <- sds
  set.seed(23)
  test_sce <- differentiationTest(sds = sce, conditions = "cond")
  expect_identical(test_sce, test)
})
test_that("The differentiationTest work on all tests",{
  # Input SlingshotDataSet
  set.seed(23)
  test <- differentiationTest(sds = sds, conditions = condition, method = "Classifier")
  expect_is(test, "data.frame")
  expect_equal(dim(test), c(1, 3))
  expect_equal(colnames(test),  c("pair", "statistic", "p.value"))
  test <- differentiationTest(sds = sds, conditions = condition, method = "mmd")
  expect_is(test, "data.frame")
  expect_equal(dim(test), c(1, 3))
  expect_equal(colnames(test),  c("pair", "statistic", "p.value"))
})
