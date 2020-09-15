.diffDifferentiationTest <- function(sds, cd, global = TRUE, lineages = FALSE,
                                     method = "KS", thresh = 0.05) {
  A <- unique(cd)[1]
  B <- unique(cd)[2]
  pst <- slingshot::slingPseudotime(sds, na.rm = TRUE)
  w <- slingshot::slingCurveWeights(sds)
  w <- sweep(w, 1, FUN = "/", STATS = apply(w, 1, sum))
  lineages_test <- lapply(seq_len(ncol(pst)), function(l){
    w_l <- w[, l]
    pst_l <- pst[, l]
    if (method == "KS") {
      test_l <- ks_test(x = pst_l[cd == A], w_x = w_l[cd == A],
                        y = pst_l[cd == B], w_y = w_l[cd == B],
                        thresh = thresh)
      return(c("Pval" = test_l$p.value, "Statistic" = test_l$statistic))
    } else if (method == "Permutation") {
      d_l <- weighted.mean(pst_l[cd== A], w_l[cd == A]) -
        weighted.mean(pst_l[cd== B], w_l[cd == B])
      d_il <- replicate(1e4, {
        cd_i <- sample(cd)
        return(weighted.mean(pst_l[cd_i== A], w_l[cd_i == A]) -
                 weighted.mean(pst_l[cd_i== B], w_l[cd_i == B]))
      })
      return(c("Pval" = mean(abs(d_l) > abs(d_il)),
               "Statistic" = d_l))
    } else {
      stop("Method must be one of KS or permutation")
    }
  }) %>%
    dplyr::bind_rows(.id = "Lineage") %>%
    dplyr::mutate(Lineage = as.character(Lineage))
  glob_test <- .Stouffer(pvals = lineages_test$Pval,
                         weights = colSums(w))
  glob_test <- data.frame("Pval" = glob_test$Pval,
                          "Statistic" = glob_test$Statistic,
                          "Lineage" = "All")
  if (global == TRUE & lineages == FALSE) return(glob_test)
  if (global == FALSE & lineages == TRUE) return(lineages_test)
  if (global == TRUE & lineages == TRUE) {
    return(dplyr::bind_rows(glob_test, lineages_test))
  }
}


#' Differential Differentiation Test
#'
#' @description Test whether or not the cell repartition between lineages is
#' independent of the conditions
#'
#' @param sds The final object after running slingshot. Can be either a
#' \code{\link{slingshotDataset}} or a \code{\link{SingleCellExperiment}} object.
#' @param cd Either the vector of conditions, or a character indicating which
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
          definition = function(sds, cd, global = TRUE, lineages = FALSE,
                                method = "KS", thresh = .05){
            if (n_distinct(cd) != 2) {
              stop("For now, this test only works with two conditions")
            }
            res <- .diffProgressionTest(sds = sds,
                                        cd = cd,
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
          definition = function(sds, cd,  global = TRUE, lineages = FALSE,
                                method = "KS", thresh = .05){
            if (is.null(sds@int_metadata$slingshot)) {
              stop("For now this only works downstream of slingshot")
            }
            if (length(cd == 1)) {
              if (cd %in% colnames(SummarizedExperiment::colData(sds))) {
                conditions <- SummarizedExperiment::colData(sds)[, cd]
              } else {
                stop("cd is not a column of colData(sds)")
              }
            }
            return(diffTopoTest(slingshot::SlingshotDataSet(sds),
                                cd = conditions,
                                global = global,
                                lineages = lineages,
                                method = method,
                                thresh = thresh))
          }
)
