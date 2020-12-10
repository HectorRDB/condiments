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


test_that("nLineages work", {
  pd <- DataFrame(cond = condition)
  rownames(pd) <- colnames(sds)
  sce <- SingleCellExperiment(assay = list(counts = t(reducedDim(sds))),
                              colData = pd)
  sce@int_metadata$slingshot <- sds
  expect_equal(nLineages(sds), 2)
  expect_equal(nLineages(sce), 2)
})

test_that("imbalance score works", {
  score <- imbalance_score(rd, conditions = condition)
  expect_equal(length(score$scores), nrow(rd))
  expect_equal(length(score$scaled_scores), nrow(rd))
  expect_equal(names(score$scaled_scores), rownames(rd))
})
