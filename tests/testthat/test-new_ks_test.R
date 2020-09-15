library(testthat)

test_that("K_S test works as expected",{
  for (i in 1:10){
    x <- runif(1000, 0, 1)
    y <- runif(1000, 0, 1)
    ks_og <- ks.test(x, y, alternative = "two.sided")
    new_ks <- ks_test(x, y, thresh = 0)
    expect_equivalent(ks_og$statistic, new_ks$statistic)
    expect_equivalent(ks_og$p.value, new_ks$p.value)
    expect_equal(ks_og$data.name, new_ks$data.name)
  }
})
