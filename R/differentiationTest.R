#' Differential differentiation
#'
#' @description Test whether or not the cell repartition between lineages is
#' independent of the conditions
#'
#' @param ... See the \code{\link{fateSelectionTest}}

#' @return See the \code{\link{fateSelectionTest}}
#' @examples
#' data('slingshotExample', package = "slingshot")
#' rd <- slingshotExample$rd
#' cl <- slingshotExample$cl
#' condition <- factor(rep(c('A','B'), length.out = nrow(rd)))
#' condition[110:139] <- 'A'
#' sds <- slingshot::slingshot(rd, cl)
#' differentiationTest(sds, condition)
#' @export

differentiationTest <- function(...) {
  return(fateSelectionTest(...))
}
