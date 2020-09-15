
.diffTopoTest <- function(sds, rep = 200, thresh = .05) {
  return()
}


#' Differential Topology Test
#'
#' @description Test whether or not slingshot should be fitted independently
#' for different conditions or not.
#'
#' @param sds A slingshot object already run on the full dataset. Can be either a
#' \code{\link{slingshotDataset}} or a \code{\link{SingleCellExperiment}} object.
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
            res <- .diffTopoTest(sds = sds, rep = rep, thresh = thresh)
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
