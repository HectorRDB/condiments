# Toy examples ----
.create_fork <- function(n_cells, noise, speed = 1) {
  rd <- matrix(0, ncol = 2, nrow = n_cells)
  speed <- speed^(1 / 4) * atanh(2 / 3)
  common <- round(n_cells - n_cells * tanh(speed))
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
#' @param speed How fast the cells from condition B should differentiate
#' @return A list with two components
#'  \itemize{
#'   \item \code{sd}: An \code{n_cells} by \code{4} dataframe that contains the
#'   reduced dimensions coordinates, lineage assignment (1 or 2) and condition
#'   assignment (A or B) for each cell.
#'   \item \code{mst}: a data.frame that contains the skeleton of the trajectories
#' }
#' @examples
#' sd <- create_differential_topology()
#' @importFrom stats runif qnorm pnorm rbinom rnorm
#' @importFrom dplyr bind_rows
#' @export
create_differential_topology <- function(n_cells = 200, noise = .15, shift = 10,
                                         unbalance_level = .9, speed = 1) {
  # Create shape
  mst_1 <- data.frame(Dim1 = c(-35, -15, 25, -35, -15, 25),
                      Dim2 = c(0, 0, 2, 0, 0, -2),
                      lineages = c(rep(1, 3), rep(-1, 3)))
  mst_2 <- mst_1
  mst_2[c(2, 5), 1] <- mst_2[c(2, 5), 1] - shift
  mst <- bind_rows("A" = mst_1,
                   "B" = mst_2,
                   .id = "conditions")
  # Generate cells
  sd_1 <- .create_fork(round(n_cells / 2), noise)
  sd_2 <- .create_fork(round(n_cells / 2), noise, speed = speed)

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
  sd <- data.frame(rbind(sd_1$rd, sd_2$rd),
                   "lineages" = c(sd_1$lineages, sd_2$lineages),
                   "conditions" = rep(c("A", "B"), each = round(n_cells / 2)))

  # Unbalance lineage 2 toward the first condition
  lineage_2 <- sd$lineages == -1 & sd$Dim1 > -10
  sd$conditions[lineage_2] <- "B"
  change <- sample(which(lineage_2),
                   size = round(sum(lineage_2) * unbalance_level))
  sd$conditions[change] <- "A"

  return(list("sd" = sd, "mst" = mst))
}
