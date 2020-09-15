#' @export
#' @name proximityScore
#' @title Proximity Score
#' @param k The number of neighbours to consider when computing the score.
#'  Default to 10.
#' @param smooth The smoothing parameter. Default to k. Lower values mean that
#' we smooth more.

setGeneric(
  name = "proximity_score",
  signature = "Object",
  def = function(Object, ...) {
    standardGeneric("proximity_score")
  }
)
