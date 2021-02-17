#' A toy dataset used in the vignette and in the examples
#'
#' This example has been created using the `create_differential_topology`
#' function.
#'
#' @format A list with two dataframes
#' \itemize{
#'   \item *sd* A dataframe containing, for 1000 cells, the dimensions in two
#'   coordinates, and cluster, lineage and condition assignment.
#'   \item \code{mst}: a data.frame that contains the skeleton of the trajectories
#' }
#' @source The following code reproduces the object
#' \code{
#' set.seed(21)
#' library(condiments)
#' data <- create_differential_topology(n_cells = 1000, shift = 0)
#' data$sd$Dim2 <- data$sd$Dim2 * 5
#' data$mst$Dim2 <- data$mst$Dim2 * 5
#' data$sd$cl <- kmeans(as.matrix(data$sd[, 1:2]), 8)$cluster
#' data$sd$cl <- as.character(data$sd$cl)
#' }
#' @usage data(toy_dataset)
"toy_dataset"
