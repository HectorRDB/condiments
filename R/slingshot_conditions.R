.clean_mst <- function(sds_cond, cluss) {
  # Lineage List
  sds_cond@lineages <- lapply(sds_cond@lineages, function(lin) {
    return(lin[!lin %in% cluss])
  })
  # Adjency matrix
  mat <- sds_cond@adjacency
  mat <- mat[!rownames(mat) %in% cluss, ]
  mat <- mat[, !colnames(mat) %in% cluss]
  sds_cond@adjacency <- mat
  # Cluster Labels
  sds_cond@clusterLabels <- sds_cond@clusterLabels[,
        !colnames(sds_cond@clusterLabels) %in% cluss]
  # Clean lineage list
  inside <- sapply(sds_cond@lineages, function(lin1) {
    sapply(sds_cond@lineages, function(lin2){
      return(all(lin1 %in% lin2))
    })
  })
  diag(inside) <- FALSE
  while (any(inside)) {
    lin <- which(rowSums(inside) > 0)[1]
    sds_cond@lineages <- sds_cond@lineages[-lin]
  }
  names(sds_cond@lineages) <- paste0("Lineage", seq_along(sds_cond@lineages))
  # Params
  sds_cond@slingParams$end.clus <- sapply(sds_cond@lineages, tail, n = 1)
  dist_mat <- sds_cond@slingParams$dist
  dist_mat <- dist_mat[!rownames(dist_mat) %in% cluss, ]
  dist_mat <- dist_mat[, !colnames(dist_mat) %in% cluss]
  sds_cond@slingParams$dist <- dist_mat
  return(sds_cond)
}

.slingshot_conditions <- function(sds, conditions, approx_points = 100, ...) {
  if (n_distinct(conditions) == 1) {
    cond <- conditions[1]
    return(list(cond = sds))
  }
  sdss <- list()
  for (cond in unique(conditions)) {
    sds_cond <- sds
    sds_cond@reducedDim <- sds_cond@reducedDim[conditions == cond, ]
    sds_cond@clusterLabels <- sds_cond@clusterLabels[conditions == cond, ]
    if (any(colSums(sds_cond@clusterLabels) == 0)) {
      cluss <- colnames(sds_cond@clusterLabels)[
        colSums(sds_cond@clusterLabels) == 0] 
      clus <- cluss[1]
      message(paste0("Cluster ", clus, " contains no cells condition", cond, ". ",
                     "This means you should either lower the clustering resolution before ",
                     "running trajectory inference or fit one trajectory per condition"))
      sds_cond <- .clean_mst(sds_cond, cluss)
    }
    sdss[[cond]] <- slingshot::getCurves(sds_cond, approx_points = approx_points,
                                         ...)
  }
  return(sdss)
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
#' @param approx_points Passed to \code{\link[slingshot]{getCurves}}
#' @param ... Other arguments passed to \code{\link[slingshot]{getCurves}}
#' @return
#' A list of \code{\link[slingshot]{SlingshotDataSet}}, one per condition.
#' @examples
#' data('slingshotExample', package = "slingshot")
#' rd <- slingshotExample$rd
#' cl <- slingshotExample$cl
#' condition <- factor(rep(c('A','B'), length.out = nrow(rd)))
#' condition[110:139] <- 'A'
#' sds <- slingshot::slingshot(rd, cl)
#' sdss <- slingshot_conditions(sds, condition)
#' @export
#' @importFrom slingshot SlingshotDataSet getCurves
#' @importFrom dplyr n_distinct
#' @rdname slingshot_conditions
setMethod(f = "slingshot_conditions",
          signature = c(sds = "SlingshotDataSet"),
          definition = function(sds, conditions, approx_points = 100, ...) {
            sdss <- .slingshot_conditions(sds = sds, conditions = conditions,
                                          approx_points = approx_points, ...)
            return(sdss)
          }
)

#' @export
#' @rdname slingshot_conditions
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData
setMethod(f = "slingshot_conditions",
          signature = c(sds = "SingleCellExperiment"),
          definition = function(sds, conditions, approx_points = 100, ...) {
            if (is.null(sds@int_metadata$slingshot)) {
              stop("For now this only works downstream of slingshot")
            }
            if (length(conditions) == 1) {
              if (conditions %in% colnames(SummarizedExperiment::colData(sds))) {
                conditions <- SummarizedExperiment::colData(sds)[, conditions]
              } else {
                stop("conditions is not a column of colData(sds)")
              }
            }
            return(slingshot_conditions(slingshot::SlingshotDataSet(sds),
                                        conditions = conditions,
                                        approx_points = approx_points, ...))
          }
)
