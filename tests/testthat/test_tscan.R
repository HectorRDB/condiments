library(testthat)
library(slingshot)
library(SingleCellExperiment)
library(TSCAN)

data(list = 'slingshotExample', package = "slingshot")
if (!"cl" %in% ls()) {
  rd <- slingshotExample$rd
  cl <- slingshotExample$cl
}
condition <- factor(rep(c('A','B'), length.out = nrow(rd)))
condition[110:139] <- 'A'
mst <- createClusterMST(rd, cl)
mapping <- try(mapCellsToEdges(rd, mst, cl))
if (class(mapping) != "try-error") {
  ordering <- pathStat(orderCells(mapping, mst, start = 1))

  test_that("Weights can be extracted from the pseudotime",{
    ws <- weights_from_pst(ordering)
    expect_equal(dim(ws), dim(ordering))
    expect_true(all(ws >= 0))
    expect_true(all(ws <= 1))
    expect_equal(is(ws), is(ordering))
    ws <- weights_from_pst(as.data.frame(ordering))
    expect_equal(dim(ws), dim(ordering))
    expect_true(all(ws >= 0))
    expect_true(all(ws <= 1))
    expect_equal(is(ws), is(ordering))
  })


  test_that("Differential progression and differentiation work",{
    ws <- weights_from_pst(ordering)
    expect_is(progressionTest(ordering, ws, condition),
              "data.frame")
    expect_is(fateSelectionTest(ws, condition),
              "data.frame")
  })
}
