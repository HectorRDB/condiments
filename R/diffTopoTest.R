.condition_sling <- function(sds, conditions) {
  psts <- NULL
  ws <- NULL
  for (cond in unique(conditions)) {
    sds_cond <- sds
    sds_cond@reducedDim <- sds_cond@reducedDim[conditions == cond, ]
    sds_cond@clusterLabels <- sds_cond@clusterLabels[conditions == cond, ]
    # TODO: handle cases where one cluster contains only one condition!!
    # missing_clusters <- which(colSums(sds_cond@clusterLabels) == 0)
    # sds_cond@clusterLabels <- sds_cond@clusterLabels[, -missing_clusters]
    # sds_cond@adjacency <- sds_cond@adjacency[-missing_clusters, -missing_clusters]
    # sds_cond@lineages <- lapply(sds_cond@lineages, function(l) {
    #   return(l[!l %in% names(missing_clusters)])
    # })
    sds_cond <- slingshot::getCurves(sds_cond, approx_points = 100)
    pst_cond <- slingshot::slingPseudotime(sds_cond)
    psts <- rbind(psts, pst_cond)
    w_conds <- slingshot::slingCurveWeights(sds_cond)
    ws <- rbind(ws, w_conds)
  }
  ws <- sweep(ws, 1, FUN = "/", STATS = apply(ws, 1, sum))
  return(list("psts" = psts, "ws" = ws))
}


.diffTopoTest <- function(sds, conditions, rep = 200, thresh = .05) {
  og <- .condition_sling(sds, conditions)
  permutations <- lapply(seq_len(rep), function(i) {
    condition <- sample(conditions)
    return(.condition_sling(sds, condition))
  })
  psts <- lapply(permutations, '[', 1) %>% unlist()
  ws <- lapply(permutations, '[', 2) %>% unlist()
  res <- ks_test(x = og$psts, w_x = og$ws,
                 y = psts, w_y = ws,
                 thresh = thresh)
  res$alternative <- "Infer separate trajectories for each condition"
  return(res)
}


#' Differential Topology Test
#'
#' @description Test whether or not slingshot should be fitted independently
#' for different conditions or not.
#'
#' @param sds A slingshot object already run on the full dataset. Can be either a
#' \code{\link{SlingshotDataSet}} or a \code{\link{SingleCellExperiment}} object.
#' @param conditions Either the vector of conditions, or a character indicating which
#' column of the metadata contains this vector.
#' @param rep How many permutations to run. Default to 200.
#' @param thresh the threshold for the KS test. See \code{\link{ks_test}}.
#' @return
#' The test output. An object of class \code{htest}
#' @importFrom slingshot SlingshotDataSet getCurves slingPseudotime slingCurveWeights
#' @examples
#' data('slingshotExample', package = "slingshot")
#' rd <- slingshotExample$rd
#' cl <- slingshotExample$cl
#' condition <- factor(rep(c('A','B'), length.out = nrow(rd)))
#' condition[110:139] <- 'A'
#' sds <- slingshot::getLineages(rd, cl)
#' diffTopoTest(sds, condition, rep = 20)
#' @export
#' @rdname diffTopoTest
setMethod(f = "diffTopoTest",
          signature = c(sds = "SlingshotDataSet"),
          definition = function(sds, conditions, rep = 200, thresh = .05){
            res <- .diffTopoTest(sds = sds, conditions = conditions, rep = rep,
                                 thresh = thresh)
            return(res)
          }
)


#' @export
#' @rdname diffTopoTest
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData
setMethod(f = "diffTopoTest",
          signature = c(sds = "SingleCellExperiment"),
          definition = function(sds, conditions = conditions, rep = 200,
                                thresh = .05){
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
                                rep = rep,
                                thresh = thresh))
          }
)
