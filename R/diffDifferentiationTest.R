.diffDifferentiationTest <- function(sds, conditions, global = TRUE,
                                     lineages = FALSE, method = "KS",
                                     thresh = 0.05) {
}


#' Differential Differentiation Test
#'
#' @description Test whether or not the cell repartition between lineages is
#' independent of the conditions
#'
#' @param sds The final object after running slingshot. Can be either a
#' \code{\link{slingshotDataset}} or a \code{\link{SingleCellExperiment}} object.
#' @param conditions Either the vector of conditions, or a character indicating which
#' column of the metadata contains this vector
#' @param global If TRUE, test for all pairs simultaneously.
#' @param pairwise If TRUE, test for all pairs independently.
#' @param method Either "KS" for the weighted Kolmogorov-Smirnov or "Permutation"
#' for a permutation. See details. Default to KS.
#' @param thresh the threshold for the KS test. See \code{\link{ks_test}}.
#' @import slingshot
#' @importFrom dplyr n_distinct bind_rows mutate
#' @importFrom magrittr %>%
#' @export
#' @rdname diffDifferentiationTest
setMethod(f = "diffDifferentiationTest",
          signature = c(sds = "SlingshotDataSet"),
          definition = function(sds, conditions, global = TRUE, lineages = FALSE,
                                method = "KS", thresh = .05){
            if (n_distinct(conditions) != 2) {
              stop("For now, this test only works with two conditions")
            }
            res <- .diffProgressionTest(sds = sds,
                                        conditions = conditions,
                                        global = global,
                                        lineages = lineages,
                                        method = method,
                                        thresh = thresh)
            return(res)
          }
)


#' @export
#' @rdname diffTopoTest
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment colData
setMethod(f = "diffProgressionTest",
          signature = c(sds = "SingleCellExperiment"),
          definition = function(sds, conditions,  global = TRUE,
                                lineages = FALSE, method = "KS", thresh = .05){
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
            return(diffTopoTest(slingshot::SlingshotDataSet(sds),
                                conditions = conditions,
                                global = global,
                                lineages = lineages,
                                method = method,
                                thresh = thresh))
          }
)
