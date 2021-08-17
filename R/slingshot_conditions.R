.message_missing_clus <- function(clus, cond) {
  message(paste0("Cluster ", clus, " contains no cells of condition ", cond,
                 ". This means you should either lower the clustering ",
                 "resolution before running trajectory inference or fit",
                 " one trajectory per condition"))
  message("Proceeding through")
  return(NULL)
}

.clean_mst <- function(sds_cond, cluss) {
  # Lineage List
  sds_cond@metadata$lineages <- lapply(slingLineages(sds_cond),
                                       function(lin) {
    return(lin[!lin %in% cluss])
  })
  # mst
  for (clus in cluss) {
    sds_cond@metadata$mst <- igraph::delete_vertices(slingMST(sds_cond), clus)
  }
  # Cluster Labels
  sds_cond@elementMetadata$clusterLabels <- slingClusterLabels(sds_cond)[,
                          !colnames(slingClusterLabels(sds_cond)) %in% cluss]
  # Clean lineage list
  inside <- sapply(slingLineages(sds_cond), function(lin1) {
    sapply(slingLineages(sds_cond), function(lin2){
      return(all(lin1 %in% lin2))
    })
  })
  diag(inside) <- FALSE
  while (any(inside)) {
    lin <- which(colSums(inside) > 0)[1]
    sds_cond@metadata$lineages <- slingLineages(sds_cond)[-lin]
    inside <- sapply(slingLineages(sds_cond), function(lin1) {
      sapply(slingLineages(sds_cond), function(lin2){
        return(all(lin1 %in% lin2))
      })
    })
    inside <- as.matrix(inside)
    diag(inside) <- FALSE
  }
  # names(sds_cond@metadata$lineages) <-
  #   paste0("Lineage", seq_along(slingLineages(sds_cond)))
  # Params
  sds_cond@metadata$slingParams$end.clus <-
    sapply(slingLineages(sds_cond), utils::tail, n = 1)
  # sds_cond <- sds_cond[, paste0("Lineage", seq_along(slingLineages(sds_cond)))]
  sds_cond <- sds_cond[, names(sds_cond@metadata$lineages)]
  return(sds_cond)
}

.sling_cond <- function(sds, conditions, approx_points = 100, verbose = TRUE,
                        ...) {
  if (n_distinct(conditions) == 1) {
    cond <- conditions[1]
    return(list(cond = sds))
  }
  sdss <- list()
  for (cond in unique(conditions)) {
    sds_cond <- sds[conditions == cond, ]
    if (any(colSums(slingClusterLabels(sds_cond)) == 0)) {
      cluss <- colnames(slingClusterLabels(sds_cond))[
        colSums(slingClusterLabels(sds_cond)) == 0]
      clus <- cluss[1]
      if (verbose) {
        .message_missing_clus(clus, cond)
      }
      sds_cond <- .clean_mst(sds_cond, cluss)
    }
    sdss[[cond]] <- slingshot::getCurves(sds_cond, approx_points = approx_points,
                                         ...)
  }
  return(sdss)
}

.recompute_skeleton <- function(sds) {
  mst <- sds@metadata$mst
  rd <- slingReducedDim(sds)
  clusters <- apply(slingClusterLabels(sds), 1, function(r) {which(r == 1)})
  centers <- base::rowsum(rd, group = clusters)
  centers <- apply(centers, 2, function(dim) {dim / as.vector(table(clusters))})
  centers <- as.list(as.data.frame(t(centers)))
  centers <- lapply(centers, function(cts) {
    names(cts) <- colnames(rd)
    return(cts)
  })
  igraph::V(mst)$coordinates <- centers
  sds@metadata$mst <- mst
  return(sds)
}

#' Conditions slingshot
#'
#' @description Based on an original slingshot object, refit one trajectory per
#' condition, using the same skeleton.
#'
#' @param sds A slingshot object already run on the full dataset. Can be either a
#' \code{\link{SlingshotDataSet}} or a \code{\link{SingleCellExperiment}} object.
#' @param conditions Either the vector of conditions, or a character indicating which
#' column of the metadata contains this vector.
#' @param approx_points Passed to \code{\link[slingshot]{getCurves}}
#' @param adjust_skeleton Boolean, default to `TRUE`. Whether to recompute the locations
#' of the nodes after fitting per conditions.
#' @param verbose Boolean, default to `TRUE`. Control whether messages are printed.
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
#' @importFrom slingshot SlingshotDataSet getCurves as.PseudotimeOrdering
#' @importFrom utils tail
#' @importFrom dplyr n_distinct
#' @rdname slingshot_conditions
setMethod(f = "slingshot_conditions",
          signature = c(sds = "SlingshotDataSet"),
          definition = function(sds, conditions, approx_points = 100,
                                adjust_skeleton = TRUE, verbose = TRUE, ...) {
            sdss <- slingshot_conditions(sds = as.PseudotimeOrdering(sds),
                                         conditions = conditions,
                                         approx_points = approx_points,
                                         adjust_skeleton = adjust_skeleton,
                                         verbose = verbose,
                                         ...)
            return(sdss)
          }
)

#' @export
#' @rdname slingshot_conditions
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData
setMethod(f = "slingshot_conditions",
          signature = c(sds = "SingleCellExperiment"),
          definition = function(sds, conditions, approx_points = 100,
                                adjust_skeleton = TRUE, verbose = TRUE, ...) {
            if (is.null(sds@int_metadata$slingshot) & is.null(colData(sds)$slingshot)) {
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
                                        approx_points = approx_points,
                                        adjust_skeleton = adjust_skeleton,
                                        verbose = verbose,
                                        ...))
          }
)


#' @rdname slingshot_conditions
#' @importClassesFrom TrajectoryUtils PseudotimeOrdering
#' @importFrom slingshot slingLineages slingClusterLabels slingMST
#' @importFrom igraph delete_vertices
#' @export
setMethod(f = "slingshot_conditions",
          signature = c(sds = "PseudotimeOrdering"),
          definition = function(sds, conditions, approx_points = 100,
                                adjust_skeleton = TRUE, verbose = TRUE, ...) {
            sdss <- .sling_cond(sds = sds, conditions = conditions,
                                approx_points = approx_points,
                                verbose = verbose,
                                ...)
            if (adjust_skeleton) {
              sdss <- lapply(sdss, .recompute_skeleton)
            }
            return(sdss)
          }
)
