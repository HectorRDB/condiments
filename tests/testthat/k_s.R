library(testthat)

test_that("K_S test works as expected",{
  x <- runif(1000, 0, 1)
  y <- runif(1000, 0, 1)
  ks_og <- ks.test(x, y, alternative = "two.sided")
  new_ks <- ks_test(x, y, thresh = 0)
  expect_equivalent(ks_og$statistic, new_ks$statistic)
  expect_equivalent(ks_og$p.value, new_ks$p.value)
  expect_equal(ks_og$data.name, new_ks$data.name)
})



function (x, y, ..., alternative = c("two.sided", "less", "greater"),
          exact = NULL)
{
  alternative <- match.arg(alternative)
  DNAME <- deparse1(substitute(x))
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 1L)
    stop("not enough 'x' data")
  PVAL <- NULL
  if (is.numeric(y)) {
    DNAME <- paste(DNAME, "and", deparse1(substitute(y)))
    y <- y[!is.na(y)]
    n.x <- as.double(n)
    n.y <- length(y)
    if (n.y < 1L)
      stop("not enough 'y' data")
    if (is.null(exact))
      exact <- (n.x * n.y < 10000)
    METHOD <- "Two-sample Kolmogorov-Smirnov test"
    TIES <- FALSE
    n <- n.x * n.y/(n.x + n.y)
    w <- c(x, y)
    z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))
    if (length(unique(w)) < (n.x + n.y)) {
      if (exact) {
        warning("cannot compute exact p-value with ties")
        exact <- FALSE
      }
      else warning("p-value will be approximate in the presence of ties")
      z <- z[c(which(diff(sort(w)) != 0), n.x + n.y)]
      TIES <- TRUE
    }
    STATISTIC <- switch(alternative, two.sided = max(abs(z)),
                        greater = max(z), less = -min(z))
    nm_alternative <- switch(alternative, two.sided = "two-sided",
                             less = "the CDF of x lies below that of y", greater = "the CDF of x lies above that of y")
    if (exact && (alternative == "two.sided") && !TIES)
      PVAL <- 1 - .Call(C_pSmirnov2x, STATISTIC, n.x, n.y)
  }
  names(STATISTIC) <- switch(alternative, two.sided = "D",
                             greater = "D^+", less = "D^-")
  if (is.null(PVAL)) {
    pkstwo <- function(x, tol = 1e-06) {
      if (is.numeric(x))
        x <- as.double(x)
      else stop("argument 'x' must be numeric")
      p <- rep(0, length(x))
      p[is.na(x)] <- NA
      IND <- which(!is.na(x) & (x > 0))
      if (length(IND))
        p[IND] <- .Call(C_pKS2, p = x[IND], tol)
      p
    }
    PVAL <- if (alternative == "two.sided")
      1 - pkstwo(sqrt(n) * STATISTIC)
    else exp(-2 * n * STATISTIC^2)
  }
  PVAL <- min(1, max(0, PVAL))
  RVAL <- list(statistic = STATISTIC, p.value = PVAL, alternative = nm_alternative,
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}
