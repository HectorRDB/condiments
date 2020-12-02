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
  signature = "sds",
  def = function(sds, ...) {
    standardGeneric("progressionTest")
  }
)

#' @name differentiationTest
#' @title Differential Differentiation Test
#' @param ... parameters including:
#' @import methods
methods::setGeneric(

  name = "differentiationTest",
  signature = "sds",
  def = function(sds, ...) {
    standardGeneric("differentiationTest")
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

