
.diffProgressionTest <- function(sds, thresh = 0.05, global = TRUE, lineages = FALSE) {

}


#' Differential Topology Test
#'
#' @description Test whether or not slingshot should be fitted independently
#' for different conditions or not.
#'
#' @param sds The final object after running slingshot. Can be either a
#' \code\link{{slingshotDataset}} or a \code\link{{SingleCellExperiment}} object.
#' @param thresh the threshold for the KS test. See \code{\link{ks_test}}.
#' @param global If TRUE, test for all lineages simultaneously.
#' @param lineages If TRUE, test for all lineages independently.
#' @import slingshot
#' @examples
#' sd <- create_differential_topology(n_cells = 200, shift = 0,
#'                                    unbalance_level = 1)
#' @export
#' @rdname diffTopoTest
setMethod(f = "diffProgressionTest",
          signature = c(sds = "SlingshotDataSet"),
          definition = function(sds, thresh = .05, global = TRUE, lineages = FALSE){
            res <- .diffProgressionTest(sds = sds,
                                        thresh = thresh,
                                        global = global,
                                        lineages = lineages)
            return(res)
          }
)


#' @export
#' @rdname diffTopoTest
#' @import SingleCellExperiment
setMethod(f = "diffProgressionTest",
          signature = c(sds = "SingleCellExperiment"),
          definition = function(sds, rep = 200, thresh = .05){
            if (is.null(sds@int_metadata$slingshot)) {
              stop("For now this only works downstream of slingshot")
            }
            return(diffTopoTest(slingshot::SlingshotDataSet(sds),
                                thresh = thresh,
                                global = global,
                                lineages = lineages))
          }
)
