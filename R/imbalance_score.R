# Inspired from the EMT package from Uwe Menzel
# Uwe Menzel (2013). EMT: Exact Multinomial Test: Goodness-of-Fit Test for Discrete
# Multivariate data. R package version 1.1. https://CRAN.R-project.org/package=EMT
# citation(EMT)

.findVectors <- function(groups, size) {
  if (groups == 1) {
    mat <- size
  }
  else {
    mat <- matrix(rep(0, groups - 1), nrow = 1)
    for (i in seq_len(size)) {
      mat <- rbind(mat, .findVectors(groups - 1, i))
    }
    mat <- cbind(mat, size - rowSums(mat))
  }
  return(mat)
}

.multinomial.test <- function(cdMatrix, groups, props) {
  size <- ncol(cdMatrix)
  eventMat <- .findVectors(length(groups), size)
  eventProb <- apply(eventMat, 1, function(x) {
    dmultinom(x,  size = size, prob = props)
  })
  pvalues <- apply(cdMatrix, 1, function(conds){
    real <- as.vector(table(factor(conds, levels = groups)))
    pObs <- stats::dmultinom(real, size, props)
    p.value <- sum(eventProb[eventProb <= pObs])
    return(p.value)
  })
  res <- -stats::qnorm(unlist(pvalues) / 2)
  return(res)
}

.imbalance_score <- function(rd, conditions, k = 10, smooth = k) {
  # Code inspired from the monocle3 package
  # https://github.com/cole-trapnell-lab/monocle3/blob/9becd94f60930c2a9b51770e3818c194dd8201eb/R/cluster_cells.R#L194
  if (length(conditions) != nrow(rd)) {
    stop("The conditions and reduced dimensions do not contain the same cells")
  }

  props <- as.vector(table(conditions) / length(conditions))
  groups <- unique(conditions)
  if (length(groups) == 1) stop("conditions should have at least 2 classes")

  # Get the graph
  # We need to add 1 because by default, nn2 counts each cell as its own
  # neighbour
  tmp <- RANN::nn2(rd, rd, k + 1, searchtype = "standard")
  neighborMatrix <- tmp[[1]]
  cdMatrix <- matrix(factor(conditions)[neighborMatrix], ncol = k + 1)

  # Get the smoothed scores
  scores <- .multinomial.test(cdMatrix, groups, props)
  scores <- unlist(scores)
  names(scores) <- rownames(rd)
  formula <- paste0("scores ~ s(",
                    paste0("rd[, ", seq_len(ncol(rd)), "], ", collapse = ""),
                    "k = smooth)")
  mm <- mgcv::gam(stats::as.formula(formula))
  scaled_scores <- mgcv::predict.gam(mm, type = "response")

  return(list("scores" = scores, "scaled_scores" = scaled_scores))
}


#' Imbalance Score
#'
#' @description Compute a imbalance score to show whether nearby cells have the
#' same condition of not
#'
#' @param Object A \code{\link{SingleCellExperiment}} object or a matrix
#' representing the reduced dimension matrix of the cells.
#' @param dimred A string or integer scalar indicating the reduced dimension
#' result in \code{reducedDims(sce)} to plot. Default to 1.
#' @param conditions Either the vector of conditions, or a character indicating which
#' column of the metadata contains this vector
#' @param k The number of neighbors to consider when computing the score.
#'  Default to 10.
#' @param smooth The smoothing parameter. Default to k. Lower values mean that
#' we smooth more.
#' @importFrom RANN nn2
#' @importFrom mgcv gam predict.gam
#' @importFrom stats dmultinom qnorm as.formula
#' @return Either a list with the \code{scaled_scores} and the \code{scores} for
#'  each cell, if input is a matrix, or the \code{\link{SingleCellExperiment}}
#'  object, wit this list in the \code{\link{colData}}.
#' @examples
#' sd <- create_differential_topology(n_cells = 200, shift = 0,
#'                                    unbalance_level = 1)
#' scores <- imbalance_score(sd$rd, sd$conditions, k = 4)
#' cols <- as.numeric(cut(scores$scaled_scores, 8))
#' plot(sd$rd[, "Dim1"], sd$rd[, "Dim2"], xlab = "Dim1", ylab = "Dim2",
#'  pch = 16, col = RColorBrewer::brewer.pal(8, "Blues")[cols])
#' @export
#' @rdname imbalance_score
setMethod(f = "imbalance_score",
          signature = c(Object = "matrix"),
          definition = function(Object, conditions, k = 10, smooth = 10){
            scores <- .imbalance_score(rd = Object,
                                       conditions = conditions,
                                       k = k, smooth = smooth)
            return(scores)
          }
)


#' @export
#' @rdname imbalance_score
#' @importFrom SummarizedExperiment colData
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SingleCellExperiment reducedDims
setMethod(f = "imbalance_score",
          signature = c(Object = "SingleCellExperiment"),
          definition = function(Object,
                                dimred = 1,
                                conditions, k = 10,
                                smooth = 10){
            if (ncol(Object) == 1) stop("The dataset only has one cell")
            if (length(conditions == 1)) {
              if (conditions %in% colnames(SummarizedExperiment::colData(Object))) {
                conditions <- SummarizedExperiment::colData(Object)[, conditions]
              } else {
                stop("conditions is not a column of colData(Object)")
              }
            }
            if (length(SingleCellExperiment::reducedDims(Object)) == 0) {
              stop("Add a reduced dimension method to Object")
            } else {
              rd <- SingleCellExperiment::reducedDims(Object)[[dimred]]
            }
            Object$scores <- as.data.frame(
              .imbalance_score(rd = rd,
                               conditions = conditions,
                               k = k, smooth = smooth)
            )
            return(Object)
          }
)
