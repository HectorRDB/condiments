
.diffTopoTest <- function(rd, cl, k = 10, smooth = k) {
  # Code inspired from the monocle3 package
  # https://github.com/cole-trapnell-lab/monocle3/blob/9becd94f60930c2a9b51770e3818c194dd8201eb/R/cluster_cells.R#L194
  if (length(cl) != nrow(rd)) {
    stop("The conditions and reduced dimensions do not contain the same cells")
  }

  props <- as.vector(table(cl) / length(cl))
  groups <- unique(cl)
  if (length(groups) == 1) stop("cl should have at least 2 classes")

  # Get the graph
  # We need to add 1 because by default, nn2 counts each cell as its own
  # neighbour
  tmp <- RANN::nn2(rd, rd, k + 1, searchtype = "standard")
  neighborMatrix <- tmp[[1]]
  # Remove each cell from being it's own neighbour
  distMatrix <- tmp[[2]][, -1]
  distMatrix <- 1 / distMatrix
  distMatrix <- distMatrix / rowSums(distMatrix)
  clMatrix <- matrix(factor(cl)[neighborMatrix], ncol = k + 1)
  simMatrix <- clMatrix == clMatrix[, 1]
  # Remove each cell from being it's own neighbour
  simMatrix <- simMatrix[, -1]

  # Get the basic scores
  scores <- rowMeans(simMatrix)

  # Get the smoothed scores
  scaled_scores <- .multinomial.test(clMatrix, groups, props)
  scaled_scores <- unlist(scaled_scores)
  names(scaled_scores) <- rownames(rd)
  formula <- paste0("scaled_scores ~ s(",
                    paste0("rd[, ", seq_len(ncol(rd)), "], ", collapse = ""),
                    "k = smooth)")
  mm <- mgcv::gam(as.formula(formula))
  scaled_scores <- predict(mm, type = "response")

  return(list("scores" = scores, "scaled_scores" = scaled_scores))
}


#' Differential Topology Test
#'
#' @description Test whether or not slingshot should be fitted independently
#' for different conditions or not.
#'
#' @param sds A slingshot object already run on the full dataset. Can be either a
#' \code\link{{slingshotDataset}} or a \code\link{{SingleCellExperiment}} object.
#' @param rep How many permutations to run. Default to 200.
#' @param thresh the threshold for the KS test. See \code{\link{ks_test}}.
#' @import slingshot
#' @examples
#' sd <- create_differential_topology(n_cells = 200, shift = 0,
#'                                    unbalance_level = 1)
#' @export
#' @rdname diffTopoTest
setMethod(f = "diffTopoTest",
          signature = c(sds = "SlingshotDataSet"),
          definition = function(sds, rep = 200, thresh = .05){

            return(res)
          }
)


#' @export
#' @rdname diffTopoTest
#' @import SingleCellExperiment
setMethod(f = "diffTopoTest",
          signature = c(sds = "SingleCellExperiment"),
          definition = function(sds, rep = 200, thresh = .05){
            if (is.null(sds@int_metadata$slingshot)) {
              stop("For now this only works downstream of slingshot")
            }
            return(diffTopoTest(slingshot::SlingshotDataSet(sds),
                                rep = rep, thresh = thresh))
          }
)
