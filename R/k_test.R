#' New KS Test
#'
#' @description Weight Kolmogorov-Smirnov Two-Sample Test,
#' with threshold for differences
#'
#' @param x Vector of values sampled from the first distribution
#' @param y Vector of values sampled from the second distribution
#' @param thresh The threshold needed to clear between the two cumulative distributions
#' @param w_x The observation weights for x
#' @param w_y The observation weights for y
#' @importFrom spatstat ewcdf
#' @examples
#'  x <- runif(100)
#'  y <- runif(100, min = .001, max = .001)
#'  .ks.test(x, y, thresh = .001)
#' @export
ks_test <-  function (x, y, thresh = .05, w_x = rep(1, length(x)),
                      w_y = rep(1, length(y))) {
  # Inspired by https://stackoverflow.com/a/55664242/13768995
  DNAME <- deparse1(substitute(x))
  DNAME <- paste(DNAME, "and", deparse1(substitute(y)))
  w_x <- w_x[!is.na(x)]
  x <- x[!is.na(x)]
  x <- x[w_x > 0]
  w_x <- w_x[w_x > 0]
  w_x <- w_x / sum(w_x)
  n.x <- length(x)
  n.x <- as.double(n.x)
  if (n.x < 1L)  stop("not enough 'x' data")
  w_y <- w_y[!is.na(y)]
  y <- y[!is.na(y)]
  y <- y[w_y > 0]
  w_y <- w_y[w_y > 0]
  w_y <- w_y / sum(w_y)
  n.y <- length(y)
  if (n.y < 1L) stop("not enough 'y' data")
  n.x <- sum(w_x)^2/sum(w_x^2)
  n.y <- sum(w_y)^2/sum(w_y^2)
  w <- n.x * n.y / (n.x + n.y)
  x_ewcdf <- spatstat::ewcdf(x, weights = w_x)
  y_ewcdf <- spatstat::ewcdf(y, weights = w_y)
  xy <- c(x, y)
  STATISTIC <- max(abs(x_ewcdf(xy) - y_ewcdf(xy)))
  pkstwo <- function(x, tol = 1e-06) {
    if (is.numeric(x))
      x <- as.double(x)
    else stop("argument 'x' must be numeric")
    p <- rep(0, length(x))
    p[is.na(x)] <- NA
    IND <- which(!is.na(x) & (x > 0))
    if (length(IND))
      p[IND] <- .Call(stats:::C_pKS2, p = x[IND], tol)
    p
  }
  PVAL <- 1 - pkstwo(sqrt(w) * STATISTIC)

  PVAL <- min(1, max(0, PVAL))
  METHOD <- paste0("Two-sample Weighted Kolmogorov-Smirnov test with threshold ", thresh)
  RVAL <- list(statistic = STATISTIC, p.value = PVAL, alternative = "two-sided",
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}
