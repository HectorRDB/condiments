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
  rownames(pd) <- rownames(sds)
  sce <- SingleCellExperiment(
    assay = list(counts = t(slingReducedDim(sds))),
    colData = pd
  )
  expect_error(nLineages(sce))
  sce@int_metadata$slingshot <- sds
  expect_equal(nLineages(sds), 2)
  expect_equal(nLineages(sce), 2)
})

test_that("imbalance score works", {
  score <- imbalance_score(rd, conditions = condition)
  expect_equal(length(score$scores), nrow(rd))
  expect_equal(length(score$scaled_scores), nrow(rd))
  expect_equal(names(score$scaled_scores), rownames(rd))
  pd <- DataFrame(cond = condition)
  rownames(pd) <- rownames(sds)
  sce <- SingleCellExperiment(
    assay = list(counts = t(slingReducedDim(sds))),
    colData = pd
  )
  reducedDim(sce, "rd") <- rd
  score1 <- imbalance_score(sce, conditions = "cond")
  score2 <- imbalance_score(sce, conditions = condition)
  expect_equal(score1, score2)
  expect_true(all(score$scores == score2$scores$scores))
  expect_true(all(score$scaled_scores == score2$scores$scaled_scores))
  expect_error(imbalance_score(rd, conditions = rep(1, nrow(rd))))
})

test_that("example work", {
  test <- create_differential_topology()
  expect_equal(dim(test$sd), c(200, 4))
  expect_equal(dim(test$mst), c(12, 4))
})


