#' @export
#' @name proximity_score
#' @title Proximity Score
#' @import methods
#' @param ... parameters including:
methods::setGeneric(
  name = "proximity_score",
  signature = "Object",
  def = function(Object, ...) {
    standardGeneric("proximity_score")
  }
)


#' @export
#' @name diffTopoTest
#' @title Differential Topology Test
#' @param ... parameters including:
#' @import methods
methods::setGeneric(
  name = "diffTopoTest",
  signature = "sds",
  def = function(sds, ...) {
    standardGeneric("diffTopoTest")
  }
)

#' @export
#' @name diffProgressionTest
#' @title Differential Progression Test
#' @param ... parameters including:
#' @import methods
methods::setGeneric(
  name = "diffProgressionTest",
  signature = "sds",
  def = function(sds, ...) {
    standardGeneric("diffProgressionTest")
  }
)

#' @name diffDifferentiationTest
#' @title Differential Differentiation Test
#' @param ... parameters including:
#' @import methods
methods::setGeneric(

  name = "diffDifferentiationTest",
  signature = "sds",
  def = function(sds, ...) {
    standardGeneric("diffDifferentiationTest")
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

