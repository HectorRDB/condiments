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
      clus <- colnames(sds_cond@clusterLabels)[
        colSums(sds_cond@clusterLabels) == 0][1]
      stop(paste0("Cluster ", clus, " contains no cells condition", cond, ". ",
                  "This means you should either lower the clustering resolution before ",
                  "running trajectory inference or fit one trajectory per condition"))
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
