.Stouffer <- function(pvals, weights) {
  Zs <- lapply(pvals, function(pval) {
    return(qnorm(pval / 2, lower.tail = FALSE))
  }) %>%
    unlist()
  Z <- sum(Zs * weights) / sqrt(sum(weights^2))
  return(c("Pval" = pnorm(Z), "Statistic" = Z))
}


.diffProgressionTest <- function(sds, cl, thresh = 0.05, global = TRUE,
                                 lineages = FALSE) {
  A <- unique(cl)[1]
  B <- unique(cl)[2]
  pst <- slingshot::slingPseudotime(sds, na.rm = TRUE)
  w <- slingshot::slingCurveWeights(sds)
  w <- sweep(w, 1, FUN = "/", STATS = apply(w, 1, sum))
  lineages_test <- lapply(seq_len(ncol(pst)), function(l){
    w_l <- w[, l]
    pst_l <- pst[, l]
    test_l <- ks_test(x = pst_l[cl == A], w_x = w_l[cl == A],
                      y = pst_l[cl == B], w_y = w_l[cl == B],
                      thresh = thresh)
    return(c("Pval" = test_l$p.value, "Statistic" = test_l$statistic))
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


#' Differential Topology Test
#'
#' @description Test whether or not slingshot should be fitted independently
#' for different conditions or not.
#'
#' @param sds The final object after running slingshot. Can be either a
#' \code\link{{slingshotDataset}} or a \code\link{{SingleCellExperiment}} object.
#' @param cl Either the vector of conditions, or a character indicating which
#' column of the metadata contains this vector
#' @param thresh the threshold for the KS test. See \code{\link{ks_test}}.
#' @param global If TRUE, test for all lineages simultaneously.
#' @param lineages If TRUE, test for all lineages independently.
#' @import slingshot
#' @importFrom dplyr n_distinct bind_rows mutate
#' @importFrom magrittr %>%
#' @examples
#' sd <- create_differential_topology(n_cells = 200, shift = 0,
#'                                    unbalance_level = 1)
#' @export
#' @rdname diffTopoTest
setMethod(f = "diffProgressionTest",
          signature = c(sds = "SlingshotDataSet"),
          definition = function(sds, cl, thresh = .05, global = TRUE, lineages = FALSE){
            if (n_distinct(cl) != 2) {
              stop("For now, this test only works with two conditions")
            }
            res <- .diffProgressionTest(sds = sds,
                                        cl = cl,
                                        thresh = thresh,
                                        global = global,
                                        lineages = lineages)
            return(res)
          }
)


#' @export
#' @rdname diffTopoTest
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment colData
setMethod(f = "diffProgressionTest",
          signature = c(sds = "SingleCellExperiment"),
          definition = function(sds, cl, rep = 200, thresh = .05){
            if (is.null(sds@int_metadata$slingshot)) {
              stop("For now this only works downstream of slingshot")
            }
            if (length(cl == 1)) {
              if (cl %in% colnames(SummarizedExperiment::colData(sds))) {
                conditions <- SummarizedExperiment::colData(sds)[, cl]
              } else {
                stop("cl is not a column of colData(sds)")
              }
            }
            return(diffTopoTest(slingshot::SlingshotDataSet(sds),
                                cl = conditions,
                                thresh = thresh,
                                global = global,
                                lineages = lineages))
          }
)
