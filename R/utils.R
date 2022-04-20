.weights_from_pst <- function(pseudotime) {
  ws <- lapply(pseudotime, is.na) %>%
    as.data.frame()
  ws <- !ws
  ws <- sweep(ws, 1, FUN = "/", STATS = apply(ws, 1, sum))
  return(ws)
}

#' Create weight matrix for methods without partial weights
#'
#' @description Most trajectory inference methods do not perform soft assignment
#' but instead assign cells to all possible lineages before a branching point,
#' and then to one or another. This function re-creates a weight matrix from those
#' matrices of pseudotime
#'
#' @return A object of the same type and dimensions as the original object, with the
#' weights for each curve and cell.
#' @examples
#' data(list = 'slingshotExample', package = "slingshot")
#' if (!"cl" %in% ls()) {
#'   rd <- slingshotExample$rd
#'   cl <- slingshotExample$cl
#' }
#' sds <- slingshot::slingshot(rd, cl)
#' weights_from_pst(slingshot::slingPseudotime(sds))
#' @export
#' @rdname weights_from_pst
setMethod(f = "weights_from_pst",
          signature = c(pseudotime = "matrix"),
          definition = function(pseudotime){
            return(as.matrix(.weights_from_pst(as.data.frame(pseudotime))))
          }
)

#' @export
#' @rdname weights_from_pst
setMethod(f = "weights_from_pst",
          signature = c(pseudotime = "data.frame"),
          definition = function(pseudotime){
            return(.weights_from_pst(pseudotime))
          }
)
# Distinct help ----
.distinct_inputs <- function(x, distinct_samples, conditions) {
  colData <- data.frame(Samples = distinct_samples,
                        Cluster = 1,
                        "conditions" = conditions)
  sce <- SummarizedExperiment::SummarizedExperiment(
    assays = list("Pseudotime" = matrix(c(x, x), nrow = 2, byrow = TRUE)),
    colData = colData
  )
  design <- colData %>%
    dplyr::select(Samples, conditions) %>%
    dplyr::distinct() %>%
    dplyr::arrange(Samples) %>%
    model.matrix(~conditions, .)
  return(list(sce = sce, design = design))
}

# Sds merge ----

#' Merge slingshots datasets
#'
#' @description If trajectory inference needs to be manually done condition per condition,
#' this allows to merge them into one. It requires manual mapping of lineages.
#'
#' @param ... Slingshot datasets
#' @param mapping a matrix, one column per dataset. Each row amounts to lineage mapping.
#' @param condition_id A vector of condition for each condition. Default to integer values
#' in order of appearance
#' @param scale If TRUE (default), lineages that are mapped are scaled to have the same
#' length.
#' @return A modified slingshot dataset that can be used for downstream steps.
#' @import slingshot
#' @importFrom dplyr bind_rows
#' @importClassesFrom TrajectoryUtils PseudotimeOrdering
#' @details The function assumes that each lineage in a dataset maps to exactly one lineage
#' in another dataset. Anything else needs to be done manually.
#' @examples
#' data(list = 'slingshotExample', package = "slingshot")
#' if (!"cl" %in% ls()) {
#'   rd <- slingshotExample$rd
#'   cl <- slingshotExample$cl
#' }
#' sds <- slingshot::slingshot(rd, cl)
#' merge_sds(sds, sds, mapping = matrix(c(1, 2, 1, 2), nrow = 2))
#' @export
merge_sds <- function(..., mapping, condition_id = seq_len(ncol(mapping)),
                      scale = FALSE) {
  sdss <- list(...)
  names(sdss) <- condition_id
  # Checking inputs ----
  ises <- unlist(lapply(sdss, function(sds){
    any(c("SingleCellExperiment", "SlingshotDataSet", "PseudotimeOrdering") %in% is(sds))
  }))
  if (!all(ises)) {
    stop("The datasets must either be SlingshotDataset or SingleCellExperiment objects")
  }
  sdss <- lapply(sdss, as.PseudotimeOrdering)
  # Add something if one is null
  if (ncol(mapping) != length(sdss)) {
    stop("mapping should have one column per dataset")
  }
  nCurves <- lapply(sdss, function(sds) {
    length(slingCurves(sds))
  })
  if (length(unique(unlist(nCurves))) != 1) {
    stop("All datasets must have the same number of lineages")
  }
  nCurves <- nCurves[[1]]
  if (nrow(mapping) != nCurves) stop("Some lineages are not mapped")
  mapped <- apply(mapping, 2, function(maps) {
    sort(maps) == seq_len(nCurves)
  })
  if (!all(mapped)) stop("Some lineages are not mapped")

  # Merging per say ----
  sds <- Reduce('rbind', sdss)
  mapping <- as.data.frame(t(mapping))
  sds@metadata$curves <- lapply(mapping, function(map) {
    curves_i <- Map(function(n_lin, n_dataset) {
      slingCurves(sdss[[n_dataset]])[[n_lin]]
    }, n_lin = map, n_dataset = seq_along(map))
    curve_i <- curves_i[[1]]
    min_ref <- min(curve_i[["lambda"]])
    length_ref <- length(curve_i[["lambda"]])
    curve_i$s <- do.call('rbind', lapply(curves_i, "[[", "s"))
    curve_i$ord <- do.call('c', lapply(curves_i, "[[", 'ord'))
    curve_i$lambda <- do.call('c',lapply(curves_i, function(crv){
     lambda <- crv[["lambda"]]
     if (scale) {
       lambda <- (lambda - min(lambda)) / (max(lambda) - min(lambda))
       lambda <- lambda * length_ref + min_ref
     }
     return(lambda)
    }))
    curve_i$dist_ind <- do.call('c', lapply(curves_i, "[[", 'dist_ind'))
    curve_i$dist <- do.call('sum', lapply(curves_i, "[[", 'dist'))
    curve_i$w <- do.call('c', lapply(curves_i, "[[", 'w'))
    return(curve_i)
  })
  return(sds)
}


# reassign cells
.sling_reassign <- function(sds) {
  # from slingshot package
  W <- slingCurveWeights(sds)
  D <- vapply(slingCurves(sds), function(p){ p$dist_ind }, rep(0, nrow(sds)))
  ordD <- order(D)
  W.prob <- W / rowSums(W)
  WrnkD <- cumsum(W.prob[ordD]) / sum(W.prob)
  Z <- D
  Z[ordD] <- WrnkD
  Z.prime <- 1 - Z^2
  Z.prime[W == 0] <- NA
  W0 <- W
  W <- Z.prime / matrixStats::rowMaxs(Z.prime,na.rm = TRUE) #rowMins(D) / D
  W[is.nan(W)] <- 1 # handle 0/0
  W[is.na(W)] <- 0
  W[W > 1] <- 1
  W[W < 0] <- 0
  W[W0 == 0] <- 0
  # add if z < .5
  idx <- Z < .5
  W[idx] <- 1 #(rowMins(D) / D)[idx]
  # drop if z > .9 and w < .1
  ridx <- rowMaxs(Z, na.rm = TRUE) > .9 &
    rowMins(W, na.rm = TRUE) < .1
  W0 <- W[ridx, ]
  Z0 <- Z[ridx, ]
  W0[!is.na(Z0) & Z0 > .9 & W0 < .1] <- 0
  W[ridx, ] <- W0
  return(W)
}
