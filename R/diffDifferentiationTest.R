.diffDifferentiationTest <- function(sds, conditions, global = TRUE,
                                     pairwise = FALSE, method = "Permutation") {
  pairs <- utils::combn(length(slingshot::slingLineages(sds)), 2)
  df <- apply(pairs, 1, function(pair) {
    ws <- slingshot::slingCurveWeights(sds, as.probs = TRUE)[, pair]

    return()
  })
  return()
}


#' Differential Differentiation Test
#'
#' @description Test whether or not the cell repartition between lineages is
#' independent of the conditions
#'
#' @param sds The final object after running slingshot. Can be either a
#' \code{\link{SlingshotDataSet}} or a \code{\link{SingleCellExperiment}} object.
#' @param conditions Either the vector of conditions, or a character indicating which
#' column of the metadata contains this vector
#' @param global If TRUE, test for all pairs simultaneously.
#' @param pairwise If TRUE, test for all pairs independently.
#' @param method For now, only "Permutation" is accepted.
#' @importFrom slingshot slingshot SlingshotDataSet
#' @importFrom utils combn
#' @return A data frame with 3 columns:
#' \itemize{
#'   \item *Pval* the pvalue for the test at the global or lineage level
#'   \item *Statistic* for individual lineages, either the modified KS statistic
#'   if \code{method = "KS"}, or the weighted difference of means, if
#'   \code{method = "Permutation"}. For the global test, the combined Z-score.
#'   \item *Pair* for individual pairs, the lineages numbers. For global,
#'   \code{"All"}.
#' }
#' @examples
#' data('slingshotExample', package = "slingshot")
#' rd <- slingshotExample$rd
#' cl <- slingshotExample$cl
#' condition <- factor(rep(c('A','B'), length.out = nrow(rd)))
#' condition[110:139] <- 'A'
#' sds <- slingshot::slingshot(rd, cl)
#' diffDifferentiationTest(sds, condition)
#' @rdname diffDifferentiationTest
setMethod(f = "diffDifferentiationTest",
          signature = c(sds = "SlingshotDataSet"),
          definition = function(sds, conditions, global = TRUE, pairwise = FALSE,
                                method = "Permutation"){
            res <- .diffDifferentiationTest(sds = sds,
                                            conditions = conditions,
                                            global = global,
                                            pairwise = pairwise,
                                            method = method)
            return(res)
          }
)


#' @rdname diffDifferentiationTest
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData
setMethod(f = "diffDifferentiationTest",
          signature = c(sds = "SingleCellExperiment"),
          definition = function(sds, conditions,  global = TRUE,
                                pairwise = FALSE, method = "Permutation"){
            if (is.null(sds@int_metadata$slingshot)) {
              stop("For now this only works downstream of slingshot")
            }
            if (length(conditions == 1)) {
              if (conditions %in% colnames(SummarizedExperiment::colData(sds))) {
                conditions <- SummarizedExperiment::colData(sds)[, conditions]
              } else {
                stop("conditions is not a column of colData(sds)")
              }
            }
            return(diffDifferentiationTest(slingshot::SlingshotDataSet(sds),
                                           conditions = conditions,
                                           global = global,
                                           pairwise = pairwise,
                                           method = method))
          }
)
