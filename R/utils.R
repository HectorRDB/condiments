# Toy examples ----
.create_fork <- function(n_cells, noise) {
  rd <- matrix(0, ncol = 2, nrow = n_cells)
  common <- round(n_cells / 3)
  rd[, 1] <- c(stats::runif(common, -40, -10),
               stats::runif(n_cells - common, -10, 30))
  lineages <- stats::rbinom(n_cells, 1, .5) * 2 - 1
  for (i in seq_len(n_cells)) {
    rd[i, 2] <- tanh(rd[i, 1] / 10) * lineages[i] +
      rnorm(n = 1, mean = lineages[i], sd = noise)
  }
  colnames(rd) <- c("Dim1", "Dim2")
  return(list("rd" = rd, "lineages" = lineages))
}

#' Create Example function
#'
#' @description This creates a simulated reduced dimension dataset
#'
#' @param n_cells The number of cells in the dataset.
#' @param noise Amount of noise. Between 0 and 1.
#' @param shift How much should the top lineage shift in condition B.
#' @param unbalance_level How much should the bottom lineage be unbalanced toward
#' condition A.
#' @return A list with three components
#'  \itemize{
#'   \item \code{rd}: The reduced dimensions coordinates of every cells. An
#'   \code{n_cells} by \code{2} matrix
#'   \item \code{lineages}: A vector of length \code{n_cells}. Either 1 or 2
#'   depending of lineage assignment for each cell.
#'   \item \code{conditions}: A vector of length \code{n_cells}. Either A or B
#'   depending of condition assignment for each cell.
#' }
#' @examples
#' sd <- create_differential_topology()
#' @importFrom stats runif qnorm pnorm rbinom rnorm
#' @export
create_differential_topology <- function(n_cells = 200, noise = .15, shift = 10,
                                         unbalance_level = .9) {
  sd_1 <- .create_fork(round(n_cells / 2), noise)
  sd_2 <- .create_fork(round(n_cells / 2), noise)

  # Shift lineage 1 in the second condition
  shifted <- sd_2$rd[sd_2$lineages == 1, 1] - shift
  to_end <- which(shifted < -40)
  sd_2$rd[sd_2$lineages == 1, 1] <- shifted
  for (i in to_end) {
    new_dim_1 <- shifted[i] + 60 + shift
    new_dim_2 <- tanh(new_dim_1 / 10) + stats::rnorm(n = 1, mean = 1, sd = noise)
    sd_2$rd[sd_2$lineages == 1, 1][i] <- new_dim_1
    sd_2$rd[sd_2$lineages == 1 , 2][i] <- new_dim_2
  }
  sd <- list(rd = rbind(sd_1$rd, sd_2$rd),
             lineages = c(sd_1$lineages, sd_2$lineages),
             conditions = rep(c("A", "B"), each = round(n_cells / 2)))

  # Unbalance lineage 2 toward the first condition
  lineage_2 <- sd$lineages == -1 & sd$rd[, 1] > -10
  sd$conditions[lineage_2] <- "B"
  change <- sample(which(lineage_2),
                   size = round(sum(lineage_2) * unbalance_level))
  sd$conditions[change] <- "A"

  return(sd)
}

# Sds merge ----

#' Merge slingshots datasets
#'
#' @description If trajectory inference needs to be manually done condition per condition,
#' this allows to merge them into one. It requires manual mapping of lineages.
#'
#' @param ... Slingshot datasets
#' @param mapping a matrix, one column per dataset. Each row amounts to lineage mapping.
#' @param condition_id A vector of condition for each condition. Default to integer values
#' in order of appearance
#' @param scale If TRUE (default), lineages that are mapped are scaled to have the same
#' lenght.
#' @return A modified slingshot dataset that can be used for downstream steps.
#' @import slingshot
#' @importFrom dplyr bind_rows
#' @details The function assumes that each lineage in a dataset maps to exactly one lineage
#' in another dataset. Anything else needs to be done manually.
#' @export
merge_sds <- function(..., mapping, condition_id = seq_len(ncol(mapping)),
                      scale = TRUE) {
  sdss <- list(...)
  names(sdss) <- condition_id
  # Checking inputs ----
  classes <- unlist(lapply(sdss, class))
  if (!all(classes %in% c("SingleCellExperiment", "SlingshotDataSet"))) {
    stop("The datasets must either be SlingshotDataset or SingleCellExperiment objects")
  }
  sdss <- lapply(sdss, SlingshotDataSet)
  # Add something if one is null
  if (ncol(mapping) != length(sdss)) {
    stop("mapping should have one column per dataset")
  }
  nCurves <- lapply(sdss, function(sds) {
    length(slingCurves(sds))
  })
  if (length(unique(unlist(nCurves))) != 1) {
    stop("All datasets must have the same number of lineages")
  }
  nCurves <- nCurves[[1]]
  if (nrow(mapping) != nCurves) stop("Some lineages are not mapped")
  mapped <- apply(mapping, 2, function(maps) {
    sort(maps) == seq_len(nCurves)
  })
  if (!all(mapped)) stop("Some lineages are not mapped")

  # Merging per say ----
  sds <- sdss[[1]]
  sds@reducedDim <- do.call('rbind', lapply(sdss, reducedDim))
  sds@clusterLabels <- do.call('rbind', lapply(sdss, slingClusterLabels))
  mapping <- as.data.frame(t(mapping))
  sds@curves <- lapply(mapping, function(map) {
    curves_i <- Map(function(n_lin, n_dataset) {
      slingCurves(sdss[[n_dataset]])[[n_lin]]
    }, n_lin = map, n_dataset = seq_along(map))
    curve_i <- curves_i[[1]]
    min_ref <- min(curve_i[["lambda"]])
    length_ref <- length(curve_i[["lambda"]])
    curve_i$s <- do.call('rbind', lapply(curves_i, "[[", "s"))
    curve_i$ord <- do.call('c', lapply(curves_i, "[[", 'ord'))
    curve_i$lambda <- do.call('c',lapply(curves_i, function(crv){
     lambda <- crv[["lambda"]]
     if (scale) {
       lambda <- (lambda - min(lambda)) / (max(lambda) - min(lambda))
       lambda <- lambda * length_ref + min_ref
      }
    }))
    curve_i$dist_ind <- do.call('c', lapply(curves_i, "[[", 'dist_ind'))
    curve_i$dist <- do.call('sum', lapply(curves_i, "[[", 'dist'))
    curve_i$w <- do.call('c', lapply(curves_i, "[[", 'w'))
    return(curve_i)
  })
  return(sds)
}
