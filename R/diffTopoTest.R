.condition_sling <- function(sds, cd) {
  psts <- NULL
  ws <- NULL
  for (cond in unique(cd)) {
    sds_cond <- sds
    sds_cond@reducedDim <- sds_cond@reducedDim[cd == cond, ]
    sds_cond@clusterLabels <- sds_cond@clusterLabels[cd == cond, ]
    # TODO: handle cases where one cluster copntains only one condition!!
    # missing_clusters <- which(colSums(sds_cond@clusterLabels) == 0)
    # sds_cond@clusterLabels <- sds_cond@clusterLabels[, -missing_clusters]
    # sds_cond@adjacency <- sds_cond@adjacency[-missing_clusters, -missing_clusters]
    # sds_cond@lineages <- lapply(sds_cond@lineages, function(l) {
    #   return(l[!l %in% names(missing_clusters)])
    # })
    sds_cond <- slingshot::getCurves(sds_cond)
    pst_cond <- slingshot::slingPseudotime(sds_cond) %>% as.vector()
    psts <- c(psts, pst_cond)
    w_conds <- slingshot::slingCurveWeights(sds_cond)
    w_conds <- sweep(w_conds, 1, FUN = "/", STATS = apply(w_conds, 1, sum)) %>%
      as.vector()
    ws <- c(ws, w_conds)
  }
  return(list("psts" = psts, "ws" = ws))
}


.diffTopoTest <- function(sds, cd, rep = 200, thresh = .05) {
  og <- .condition_sling(sds, cd)
  permutations <- lapply(seq_len(rep), function(i) {
    condition <- sample(cd)
    return(.condition_sling(sds, condition))
  })
  psts <- lapply(permutations, '[', 1) %>% unlist()
  ws <- lapply(permutations, '[', 2) %>% unlist()
  return(ks_test(x = og$psts, w_x = og$ws,
                 y = psts, w_y = ws,
                 thresh = thresh))
}


#' Differential Topology Test
#'
#' @description Test whether or not slingshot should be fitted independently
#' for different conditions or not.
#'
#' @param sds A slingshot object already run on the full dataset. Can be either a
#' \code{\link{slingshotDataset}} or a \code{\link{SingleCellExperiment}} object.
#' @param cd Either the vector of conditions, or a character indicating which
#' column of the metadata contains this vector.
#' @param rep How many permutations to run. Default to 200.
#' @param thresh the threshold for the KS test. See \code{\link{ks_test}}.
#' @import slingshot
#' @examples
#' data('slingshotExample')
#' rd <- slingshotExample$rd
#' cl <- slingshotExample$cl
#' condition <- factor(rep(c('A','B'), length.out = nrow(rd)))
#' condition[110:139] <- 'A'
#' sds <- getLineages(rd, cl)
#' diffTopoTest(sds, condition, rep = 20)
#' @export
#' @rdname diffTopoTest
setMethod(f = "diffTopoTest",
          signature = c(sds = "SlingshotDataSet"),
          definition = function(sds, cd, rep = 200, thresh = .05){
            res <- .diffTopoTest(sds = sds, cd = cd, rep = rep, thresh = thresh)
            return(res)
          }
)


#' @export
#' @rdname diffTopoTest
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment colData
setMethod(f = "diffTopoTest",
          signature = c(sds = "SingleCellExperiment"),
          definition = function(sds, cd = cd, rep = 200, thresh = .05){
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
            return(diffTopoTest(slingshot::SlingshotDataSet(sds), cd = cd,
                                rep = rep, thresh = thresh))
          }
)
