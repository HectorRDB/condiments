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

test_that("Sds does merge correctly",{
    sds_merged <- merge_sds(sds, sds, scale = FALSE,
        mapping = matrix(c(1, 1, 2, 2), nrow = 2, byrow = TRUE))
    expect_equal(ncol(slingReducedDim(sds_merged)), ncol(slingReducedDim(sds)))
    expect_equal(nrow(slingReducedDim(sds_merged)), 2 * nrow(slingReducedDim(sds)))
    expect_equal(slingClusterLabels(sds_merged)[rownames(sds), ],
                 slingClusterLabels(sds))
    expect_equal(quantile(slingPseudotime(sds_merged, na = FALSE)[, 1]),
                 quantile(slingPseudotime(sds, na = FALSE)[, 1]))
    expect_equal(quantile(slingPseudotime(sds_merged, na = FALSE)[, 2]),
                 quantile(slingPseudotime(sds, na = FALSE)[, 2]))
})

test_that("Sds fail when it should", {
  expect_error(merge_sds(sds, slingReducedDim(sds), scale = FALSE,
                          mapping = matrix(c(1, 1, 2, 2), nrow = 2,
                                           byrow = TRUE)))
  expect_error(merge_sds(sds, sds, scale = FALSE, condition_id = 1:2,
                         mapping = matrix(c(1, 1, 2, 2), nrow = 1,
                                          byrow = TRUE)))
  sds2 <- sds
  sds2@metadata$lineages$Lineage2 <- NULL
  expect_error(merge_sds(sds, sds2, scale = FALSE,
                         mapping = matrix(c(1, 1, 2, 2), nrow = 1,
                                          byrow = TRUE)))
  expect_error(merge_sds(sds, sds, scale = FALSE,
                         mapping = matrix(c(1, 1, 3, 2), nrow = 2,
                                          byrow = TRUE)))
})

test_that("sds is fitted ok", {
  sdss <- slingshot_conditions(sds, condition)
  expect_is(sdss, "list")
  expect_is(sdss[[1]], "PseudotimeOrdering")
  expect_is(sdss[[2]], "PseudotimeOrdering")
  sdss <- slingshot_conditions(sds, rep(1, nrow(rd)))
  expect_equal(sdss[[1]], sds)
  # Input SingleCellExperiment
  pd <- DataFrame(cond = condition)
  rownames(pd) <- rownames(sds)
  sce <- SingleCellExperiment(
    assay = list(counts = t(slingReducedDim(sds))),
    colData = pd
  )
  sce@int_metadata$slingshot <- sds
  set.seed(23)
  sdss1 <- slingshot_conditions(sds, condition)
  sdss2 <- slingshot_conditions(sce, condition)
  expect_equal(sdss1, sdss2)
  sdss2 <- slingshot_conditions(sce, "cond")
  expect_equal(sdss1, sdss2)
})

test_that("missing cluster", {
  cl2 <- cl
  cl2[condition == "B" & cl == 4] <- 3
  sds <- slingshot(rd, cl2)
  sdss <- slingshot_conditions(sds, condition)
  expect_is(sdss, "list")
  expect_is(sdss[[1]], "PseudotimeOrdering")
  expect_is(sdss[[2]], "PseudotimeOrdering")
})
