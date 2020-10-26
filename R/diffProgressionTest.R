.diffProgressionTest <- function(sds, conditions, global = TRUE, lineages = FALSE,
                                 method = "KS", thresh = 0.05, rep = 1e4) {
  A <- unique(conditions)[1]
  B <- unique(conditions)[2]
  pst <- slingshot::slingPseudotime(sds, na = TRUE)
  w <- slingshot::slingCurveWeights(sds, as.probs = TRUE)
  colnames(pst) <- colnames(w) <-
    paste0("Lineage", seq_len(ncol(pst)))
  lineages_test <- lapply(colnames(pst), function(l){
    w_l <- w[, l]
    pst_l <- pst[, l]
    if (method == "KS") {
      test_l <- ks_test(x = pst_l[conditions == A],
                        w_x = w_l[conditions == A],
                        y = pst_l[conditions == B],
                        w_y = w_l[conditions == B],
                        thresh = thresh)
      return(c("Pval" = test_l$p.value, "Statistic" = test_l$statistic))
    } else if (method == "Permutation") {
      d_l <- stats::weighted.mean(pst_l[conditions== A], w_l[conditions == A]) -
        stats::weighted.mean(pst_l[conditions== B], w_l[conditions == B])
      d_il <- replicate(rep, {
        conditions_i <- sample(conditions)
        return(stats::weighted.mean(pst_l[conditions_i== A], w_l[conditions_i == A]) -
                 stats::weighted.mean(pst_l[conditions_i== B], w_l[conditions_i == B]))
      })
      return(c("Pval" = max(mean(abs(d_l) > abs(d_il)), 1 / rep),
               "Statistic" = d_l))
    } else {
      stop("Method must be one of KS or permutation")
    }
  }) %>%
    dplyr::bind_rows(.id = "Lineage") %>%
    dplyr::mutate(Lineage = as.character(Lineage)) %>%
    dplyr::select(Lineage, Pval, Statistic)
  glob_test <- .Stouffer(pvals = lineages_test$Pval,
                         weights = colSums(w))
  glob_test <- data.frame("Lineage" = "All",
                          "Pval" = glob_test$Pval,
                          "Statistic" = glob_test$Statistic)
  if (global == TRUE & lineages == FALSE) return(glob_test)
  if (global == FALSE & lineages == TRUE) return(lineages_test)
  if (global == TRUE & lineages == TRUE) {
    return(dplyr::bind_rows(glob_test, lineages_test))
  }
}


#' Differential Progression Test
#'
#' @description Test whether or not the pseudotime distribution are identical
#' within lineages between conditions
#'
#' @param sds The final object after running slingshot. Can be either a
#' \code{\link{SlingshotDataSet}} or a \code{\link{SingleCellExperiment}} object.
#' @param conditions Either the vector of conditions, or a character indicating
#' which column of the metadata contains this vector.
#' @param global If TRUE, test for all lineages simultaneously.
#' @param lineages If TRUE, test for all lineages independently.
#' @param method Either "KS" for the weighted Kolmogorov-Smirnov or "Permutation"
#' for a permutation. See details. Default to KS.
#' @param thresh The threshold for the KS test. See \code{\link{ks_test}}.
#' Ignored if \code{method = "Permutation"}. Default to .05.
#' @param rep Number of permutations to run. Ignored if \code{method = "KS"}.
#' Default to \code{1e4}.
#' @importFrom slingshot slingshot SlingshotDataSet slingPseudotime slingCurveWeights
#' @importFrom stats weighted.mean
#' @importFrom dplyr n_distinct bind_rows mutate select
#' @details
#' For every lineage, we compare the pseudotimes of the cells from either
#' conditions, using the lineage weights as observations weights.
#' \itemize{
#'   \item If \code{method = "KS"}, this uses the updated KS test,
#'   see \code{\link{ks_test}} for details.
#'   \item If \code{method = "Permutation"}, the difference of weighted mean
#'   pseudotime between condition is computed, and a p-value is found by
#'   permuting the condition labels.
#' }
#' Lineage-levels p-values are combined using Stouffer's Z-score method, using
#' the sum of cellweights per lineage to weight each p-value.
#' @references  Stouffer, S.A.; Suchman, E.A.; DeVinney, L.C.; Star, S.A.;
#' Williams, R.M. Jr. (1949).
#' *The American Soldier, Vol.1: Adjustment during Army Life.*
#' Princeton University Press, Princeton.
#' @md
#' @return A data frame with 3 columns:
#' \itemize{
#'   \item *Pval* the pvalue for the test at the global or lineage level
#'   \item *Statistic* for individual lineages, either the modified KS statistic
#'   if \code{method = "KS"}, or the weighted difference of means, if
#'   \code{method = "Permutation"}. For the global test, the combined Z-score.
#'   \item *Lineage* for individual lineages, the lineage number. For global,
#'   \code{"All"}.
#' }
#' @examples
#' data('slingshotExample', package = "slingshot")
#' rd <- slingshotExample$rd
#' cl <- slingshotExample$cl
#' condition <- factor(rep(c('A','B'), length.out = nrow(rd)))
#' condition[110:139] <- 'A'
#' sds <- slingshot::slingshot(rd, cl)
#' diffProgressionTest(sds, condition)
#' @export
#' @rdname diffProgressionTest
setMethod(f = "diffProgressionTest",
          signature = c(sds = "SlingshotDataSet"),
          definition = function(sds, conditions, global = TRUE, lineages = FALSE,
                                method = "KS", thresh = .05, rep = 1e4){
            if (n_distinct(conditions) != 2) {
              stop("For now, this test only works with two conditions")
            }
            res <- .diffProgressionTest(sds = sds,
                                        conditions = conditions,
                                        global = global,
                                        lineages = lineages,
                                        method = method,
                                        thresh = thresh,
                                        rep = rep)
            return(res)
          }
)


#' @export
#' @rdname diffProgressionTest
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData
setMethod(f = "diffProgressionTest",
          signature = c(sds = "SingleCellExperiment"),
          definition = function(sds, conditions,  global = TRUE, lineages = FALSE,
                                method = "KS", thresh = .05, rep = 1e4){
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
            return(diffProgressionTest(slingshot::SlingshotDataSet(sds),
                                       conditions = conditions,
                                       global = global,
                                       lineages = lineages,
                                       method = method,
                                       thresh = thresh,
                                       rep = rep))
          }
)
