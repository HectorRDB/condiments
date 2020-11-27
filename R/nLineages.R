#' nLineages
#'
#' @description Return the number of lineages for a slingshot object
#' @param sds A slingshot object already run on the full dataset. Can be either a
#' \code{\link[slingshot]{SlingshotDataSet}} or a
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' @export
#' @examples
#' data(list = 'slingshotExample', package = "slingshot")
#' nLineages(sds)
#' @rdname nLineages
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @return The number of lineages in the slingshot object
#' @importFrom slingshot SlingshotDataSet
setMethod(f = "nLineages",
          signature = c(sds = "SingleCellExperiment"),
          definition = function(sds){
            if (is.null(sds@int_metadata$slingshot)) {
              stop("No slingshot object")
            } else {
              return(nLineages(slingshot::SlingshotDataSet(sds)))
            }
          }
)

#' @export
#' @rdname nLineages
#' @importClassesFrom slingshot SlingshotDataSet
setMethod(f = "nLineages",
          signature = c(sds = "SlingshotDataSet"),
          definition = function(sds){
            return(length(slingshot::slingCurves(sds)))
          }
)
