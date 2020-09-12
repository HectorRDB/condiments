create_fork <- function(seed, n_cells, noise) {
  set.seed(seed)
  rd <- matrix(0, ncol = 2, nrow = n_cells)
  common <- round(n_cells / 3)
  rd[, 1] <- c(runif(common, -40, -10), runif(n_cells - common, -10, 30))
  lineages <- rbinom(n_cells, 1, .5) * 2 - 1
  for (i in seq_len(n_cells)) {
    rd[i, 2] <- tanh(rd[i, 1] / 10) * lineages[i] + 
      rnorm(n = 1, mean = lineages[i], sd = noise)
  }
  colnames(rd) <- c("Dim1", "Dim2")
  return(list("rd" = rd, "lineages" = lineages))
}

create_differential_topology <- function(seed = 5491, n_cells = 200, noise = .15,
                                shift = 10, unbalance_level = .9) {
  set.seed(seed)
  sd_1 <- create_fork(seed, round(n_cells / 2), noise)
  sd_2 <- create_fork(runif(1, -100, 100), round(n_cells / 2), noise)
  
  # Shift lineage 1 in the second condition
  shifted <- sd_2$rd[sd_2$lineages == 1, 1] - shift
  to_end <- which(shifted < -40)
  sd_2$rd[sd_2$lineages == 1, 1] <- shifted
  for (i in 1:length(to_end)) {
    new_dim_1 <- shifted[to_end[i]] + 60 + shift
    new_dim_2 <- tanh(new_dim_1 / 10) + rnorm(n = 1, mean = 1, sd = noise)
    sd_2$rd[sd_2$lineages == 1, 1][to_end[i]] <- new_dim_1
    sd_2$rd[sd_2$lineages == 1 , 2][to_end[i]] <- new_dim_2
  }
  sd <- list(rd = rbind(sd_1$rd, sd_2$rd),
             lineages = c(sd_1$lineages, sd_2$lineages),
             cl = rep(c("A", "B"), each = round(n_cells / 2)))
  
  # Unbalance lineage 2 toward the first condition
  lineage_2 <- sd$lineages == -1 & sd$rd[, 1] > -10
  sd$cl[lineage_2] <- "B"
  change <- sample(which(lineage_2),
                   size = round(sum(lineage_2) * unbalance_level))
  sd$cl[change] <- "A"

  return(sd)
}

create_differential_1D_topology <- function(seed = 5491, n_cells = 200, 
                                            noise = .15, split = 0) {
  set.seed(seed)
  rd <- matrix(0, ncol = 2, nrow = n_cells)
  common <- round(n_cells / 3)
  rd[, 1] <- c(runif(common, -40, -10), runif(n_cells - common, -10, 30))
  lineages <- rbinom(n_cells, 1, .5) * 2 - 1
  for (i in seq_len(n_cells)) {
    rd[i, 2] <- tanh(rd[i, 1] / 10) * lineages[i] * split + 
      rnorm(n = 1, mean = lineages[i] * split, sd = noise)
  }
  colnames(rd) <- c("Dim1", "Dim2")
  return(list("rd" = rd, "lineages" = lineages))
  
  return(sd)
}
