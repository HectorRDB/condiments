library(testthat)

test_that("K_S test works as expected, compared to the default ks test",{
  for (i in 1:10){
    x <- runif(1000, 0, 1)
    y <- runif(1000, 0, 1)
    ks_og <- ks.test(x, y, alternative = "two.sided")
    new_ks <- ks_test(x, y, thresh = 0)
    expect_equivalent(ks_og$statistic, new_ks$statistic)
    expect_equivalent(ks_og$p.value, new_ks$p.value, tolerance = 1e-5)
    expect_equal(ks_og$data.name, new_ks$data.name)
  }
})

test_that("K_S with equal weights equal ks with no weights",{
  for (i in 1:10){
    x <- runif(1000, 0, 1)
    y <- runif(1000, 0, 1)
    ks <- ks_test(x, y)
    ks_weights <- ks_test(x, y, w_x = rep(.5, 1000), w_y = rep(906, 1000))
    expect_equivalent(ks$statistic, ks_weights$statistic)
    expect_equivalent(ks$p.value, ks_weights$p.value)
    expect_equal(ks$data.name, ks_weights$data.name)
  }
})
