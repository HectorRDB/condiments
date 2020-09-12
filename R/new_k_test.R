#' New KS Test
#'
#' @description Kolmogorov-Smirnov Two-Sample Test, with threshold for differences
#'
#' @param x Vector of values sampled from the first distribution
#' @param y Vector of values sampled from the second distribution
#' @param thresh The threshold needed to clear between the two cumulative distributions
#' @examples
#'  x <- runif(100)
#'  y <- runif(100, min = .001, max = .001)
#'  .ks.test(x, y, thresh = .001)
#' @export
ks_test <-  function (x, y, thresh = .05) {
  DNAME <- deparse1(substitute(x))
  DNAME <- paste(DNAME, "and", deparse1(substitute(y)))
  x <- x[!is.na(x)]
  n <- length(x)
  n.x <- as.double(n)
  if (n < 1L)  stop("not enough 'x' data")
  y <- y[!is.na(y)]
  n.y <- length(y)
  if (n.y < 1L) stop("not enough 'y' data")
  w <- n.x * n.y / (n.x + n.y)
  xcdf <- ecdf(x)
  ycdf <- ecdf(y)
  dist <- lapply(c(x, y), function(a) {
    stat <- abs(xcdf(a) - ycdf(a)) - thresh
    return(max(stat, 0))
  })
  STATISTIC <- max(unlist(dist))
  PVAL <- 1 - pkstwo(sqrt(w) * STATISTIC)

  PVAL <- min(1, max(0, PVAL))
  METHOD <- paste0("Two-sample Kolmogorov-Smirnov test with threshold ", thresh)
  RVAL <- list(statistic = STATISTIC, p.value = PVAL, alternative = "two-sided",
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}
