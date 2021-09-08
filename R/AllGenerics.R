#' @export
#' @name imbalance_score
#' @title Imbalance Score
#' @import methods
#' @param ... parameters including:
methods::setGeneric(
  name = "imbalance_score",
  signature = "Object",
  def = function(Object, ...) {
    standardGeneric("imbalance_score")
  }
)

#' @export
#' @name weights_from_pst
#' @title weights_from_pst
#' @import methods
#' @param pseudotime A matrix or data.frame of \[ncells\] by \[nCurves\].
#' @param ... Other parameters including:
methods::setGeneric(
  name = "weights_from_pst",
  signature = "pseudotime",
  def = function(pseudotime, ...) {
    standardGeneric("weights_from_pst")
  }
)


#' @export
#' @name topologyTest
#' @title Differential Topology Test
#' @param ... parameters including:
#' @import methods
methods::setGeneric(
  name = "topologyTest",
  signature = "sds",
  def = function(sds, ...) {
    standardGeneric("topologyTest")
  }
)

#' @export
#' @name progressionTest
#' @title Differential Progression Test
#' @param ... parameters including:
#' @import methods
methods::setGeneric(
  name = "progressionTest",
  signature = "pseudotime",
  def = function(pseudotime, ...) {
    standardGeneric("progressionTest")
  }
)

#' @name fateSelectionTest
#' @title Differential fate selection Test
#' @param ... parameters including:
#' @import methods
methods::setGeneric(
  name = "fateSelectionTest",
  signature = "cellWeights",
  def = function(cellWeights, ...) {
    standardGeneric("fateSelectionTest")
  }
)



#' @name nLineages
#' @title Number of lineages
#' @param ... parameters including:
#' @import methods
methods::setGeneric(
  name = "nLineages",
  signature = "sds",
  def = function(sds, ...) {
    standardGeneric("nLineages")
  }
)

#' @name slingshot_conditions
#' @title Refitting slingshot per condition
#' @param ... parameters including:
#' @import methods
methods::setGeneric(
  name = "slingshot_conditions",
  signature = "sds",
  def = function(sds, ...) {
    standardGeneric("slingshot_conditions")
  }
)
