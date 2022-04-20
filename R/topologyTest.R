# Fit one slingshot object per condition ----
.condition_sling <- function(sds, conditions, verbose = TRUE) {
  sdss <- slingshot_conditions(sds, conditions, adjust_skeleton = FALSE,
                               verbose = verbose)
  psts <- lapply(sdss, slingshot::slingPseudotime, na = FALSE) %>%
    lapply(., as.data.frame) %>%
    bind_rows(., .id = "condition")
  ws <- lapply(sdss, slingshot::slingCurveWeights) %>%
    lapply(., as.data.frame) %>%
    bind_rows(., .id = "condition")
  # If a lineage is missing for some condition, we use the others to extrapolate
  uncommon <- lapply(sdss, slingLineages) %>%
    lapply(., names) %>% unlist() %>% table()
  uncommon <- names(uncommon)[uncommon != length(sdss)]
  for (lin in uncommon) {
    has_lin <- lapply(sdss, function(sds_cond){
      lin %in% names(slingLineages(sds_cond))
    }) %>% unlist()
    sdss_lin <- sdss[has_lin]
    cond_not_lin <- names(sdss)[!has_lin]
    for (cond in cond_not_lin) {
      psts[psts$condition == cond, lin] <- lapply(sdss_lin, function(sds_lin) {
        new_pst <- slingshot::predict(object = sds_lin, slingReducedDim(sdss[[cond]])) %>%
          slingPseudotime(., na = FALSE)
        new_pst <- new_pst[, lin]
      }) %>%
        Reduce(f = '+') %>%
        `/`(., sum(has_lin))
      ws[ws$condition == cond, lin] <-
        lapply(sdss_lin, function(sds_lin) {
          new_ws <- slingshot::predict(object = sds_lin, slingReducedDim(sdss[[cond]])) %>%
            slingCurveWeights()
          new_ws <- new_ws[, lin]
        }) %>%
        Reduce(f = '+') %>%
        `/`(., sum(has_lin))
    }

  }
  psts <- psts[, -which(colnames(psts) == "condition"), drop = FALSE] %>% as.matrix()
  ws <- ws[, -which(colnames(ws) == "condition"), drop = FALSE] %>% as.matrix()
  ws <- sweep(ws, 1, FUN = "/", STATS = apply(ws, 1, sum))
  return(list("psts" = psts, "ws" = ws))
}

# Run either of the methods ----
.topologyTest_ks_all <- function(permutations, og, threshs) {
  psts <- lapply(permutations, '[[', 1) %>% do.call(what = 'rbind') %>% as.vector()
  ws <- lapply(permutations, '[[', 2) %>% do.call(what = 'rbind') %>% as.vector()
  og_psts <- og$psts %>% as.vector()
  og_ws <- og$ws %>% as.vector()
  res <- lapply(threshs, function(thresh) {
    test <- ks_test(x = og_psts, w_x = og_ws,
                   y = psts, w_y = ws,
                   thresh = thresh)
    return(data.frame("statistic" = test$statistic, "p.value" = test$p.value))
  })
  names(res) <- threshs
  res <- bind_rows(res, .id = "thresh")
  return(res)
}

.topologyTest_ks_mean <- function(permutations, og, threshs, sds, rep) {
  psts <- lapply(permutations, '[[', 1) %>% Reduce(f = '+') %>% as.vector()
  psts <- psts / rep
  og_psts <- og$psts %>% as.vector()
  ws <- lapply(permutations, '[[', 2) %>%
    Reduce(f = '+') %>%
    as.vector()
  ws <- ws / rep
  og_ws <- og$ws %>% as.vector()
  res <- lapply(threshs, function(thresh) {
    test <- Ecume::ks_test(x = og_psts, w_x = og_ws,
                           y = psts, w_y = ws, thresh = thresh)
    return(data.frame("statistic" = test$statistic, "p.value" = test$p.value))
  })
  names(res) <- threshs
  res <- bind_rows(res, .id = "thresh")
  return(res)
}

.topologyTest_distinct_mean <- function(permutations, og, sds, rep, distinct_samples,
                                        conditions) {
  psts <- lapply(permutations, '[[', 1) %>% Reduce(f = '+') %>%
    as.vector()
  psts <- psts / rep
  og_psts <- og$psts %>% as.vector()
  ws <- lapply(permutations, '[[', 2) %>% Reduce(f = '+') %>% as.vector()
  ws <- ws / rep
  og_ws <- og$ws %>% as.vector()
  og_psts <- og_psts[og_ws > 0]
  psts <- psts[ws > 0]
  inputs <- .distinct_inputs(c(og_psts, psts),
                             c(distinct_samples[og_ws > 0], distinct_samples[ws > 0]),
                             c(conditions[og_ws > 0], conditions[ws > 0]))
  res <- distinct_test(x = inputs$sce, name_assays_expression = "Pseudotime",
                          name_cluster = "Cluster", name_sample = "Samples",
                          design = inputs$design,
                          column_to_test = 2, min_non_zero_cells = 0,
                          n_cores = 1)
  return(data.frame("statistic" = qnorm(res$p_val[1]), "p.value" = res$p_val[1],
                    "thresh" = "0"))
}

.topologyTest_classifier <- function(permutations, og, threshs, sds, rep,
                                     args_classifier){
  psts <- lapply(permutations, '[[', 1) %>% Reduce(f = '+')
  psts <- psts / rep
  colnames(psts) <- colnames(og$psts)
  res <- lapply(threshs, function(thresh) {
    args <- args_classifier
    args$x <- as.matrix(og$psts); args$y <- psts; args$thresh <- thresh
    test <- do.call(what = Ecume::classifier_test, args = args)
    return(data.frame("statistic" = test$statistic, "p.value" = test$p.value))
  })
  names(res) <- threshs
  res <- bind_rows(res, .id = "thresh")
  return(res)
}

.topologyTest_mmd <- function(permutations, og, sds, rep, args_mmd, n_max = 2000){
  psts <- lapply(permutations, '[[', 1) %>% Reduce(f = '+')
  psts <- psts / rep
  colnames(psts) <- colnames(og$psts)
  args <- args_mmd
  id <- sample(nrow(psts), min(nrow(psts), n_max))
  if (is.null(args_mmd$frac)) args_mmd$frac <- .5
  args$x <- as.matrix(og$psts[id, ]); args$y <- as.matrix(psts[id, ])
  test <- do.call(what = Ecume::mmd_test, args = args)
  return(data.frame("thresh" = NA,
                    "statistic" = test$statistic,
                    "p.value" = test$p.value)
         )
}

.topologyTest_wass <- function(permutations, og, sds, rep, args_wass){
  psts <- lapply(permutations, '[[', 1) %>% Reduce(f = '+')
  psts <- psts / rep
  colnames(psts) <- colnames(og$psts)
  S <- min(10^5, nrow(psts))
  args <- args_wass
  args$x <- as.matrix(og$psts); args$y <- psts; args$S <- S; args$fast <- TRUE
  test <- do.call(what = Ecume::wasserstein_permut, args = args)
  return(data.frame("thresh" = NA,
                    "statistic" = test$statistic,
                    "p.value" = test$p.value)
  )
}

# Run all selected tests ----
.topologyTest_all_selected <- function(
  permutations, og, conditions, sds,
  methods = ifelse(dplyr::n_distinct(conditions) == 2, "KS_mean", "Classifier"),
  threshs = .01, rep = 100, args_classifier = list(), args_mmd = list(),
  args_wass = list(), distinct_samples = NULL) {
  res <- list()
  if ("KS_all" %in% methods) {
    message("Running KS-all test")
    res[["KS_all"]] <- .topologyTest_ks_all(permutations, og, threshs)
  }
  if ("KS_mean" %in% methods) {
    message("Running KS-mean test")
    res[["KS_mean"]] <- .topologyTest_ks_mean(
      permutations, og, threshs, sds, rep
    )
  }
  if ("Classifier" %in% methods) {
    message("Running Classifier test")
    res[["Classifier"]] <- .topologyTest_classifier(
      permutations, og, threshs, sds, rep, args_classifier
    )
  }
  if ("distinct" %in% methods) {
    message("Running distinct test")
    res[["distinct"]] <- .topologyTest_distinct_mean(
      permutations, og, sds, rep, distinct_samples, conditions
    )
  }
  if ("mmd" %in% methods) {
    message("Running mmd test")
    res[["mmd"]] <- .topologyTest_mmd(permutations, og, sds, rep, args_mmd)
  }
  if ("wasserstein_permutation" %in% methods) {
    message("Running wassertsein permutation test")
    res[["wasserstein_permutation"]] <- .topologyTest_wass(
      permutations, og, sds, rep, args_wass
    )
  }
  return(res)
}

.generate_permutations_curves <- function(sds, conditions, rep = 100,
                                          BPPARAM = BiocParallel::bpparam(),
                                          parallel = FALSE) {
  message("Generating permuted trajectories")
  og <- .condition_sling(sds, conditions, verbose = FALSE)
  permutations <- lapply(seq_len(rep), function(i) {
    return(sample(conditions))
  })
  if (parallel) {
    permutations <- BiocParallel::bplapply(
      X = permutations, FUN = .condition_sling, sds = sds, BPPARAM = BPPARAM
    )
  } else {
    permutations <- pbapply::pblapply(
      X = permutations, FUN = .condition_sling, sds = sds
    )
  }
  order <- rownames(sds)
  if (is.null(order)) order <- rownames(reducedDim(sds))
  permutations <- lapply(permutations, function(perm){
    lapply(perm, function(dat) dat[order, ]) %>% return()
  })
  og <- lapply(og, function(dat){
    return(dat[order, ])
  })

  return(list("og" = og, "permutations" = permutations))
}
# Run the topology test ----
.topologyTest <- function(sds, conditions, rep = 100, threshs = .01,
                          methods = "KS_mean", parallel = FALSE,
                          BPPARAM = BiocParallel::bpparam(), args_mmd = list(),
                          args_classifier = list(), args_wass = list(),
                          nmax = nrow(slingshot::slingPseudotime(sds)),
                          distinct_samples = NULL) {
  curves <- .generate_permutations_curves(sds, conditions, rep, BPPARAM, parallel)
  res <- .topologyTest_all_selected(curves$permutations, curves$og, conditions,
                                    sds, methods, threshs, rep, args_classifier,
                                    args_mmd, args_wass, distinct_samples) %>%
    bind_rows(.id = "method")
  return(res)
}


#' Differential Topology Test
#'
#' @description Test whether or not slingshot should be fitted independently
#' for different conditions or not.
#'
#' @param sds A slingshot object already run on the full dataset. Can be either a
#' \code{\link{SlingshotDataSet}} or a \code{\link{SingleCellExperiment}} object.
#' @param conditions Either the vector of conditions, or a character indicating which
#' column of the metadata contains this vector.
#' @param methods The method(s) to use to test. Must be among 'KS_mean',
#' 'Classifier', "KS_all', "mmd' and 'wasserstein_permutation'. See details.
#' @param threshs the threshold(s) for the KS test or classifier test. Default to .01
#' See \code{\link{ks_test}} and \code{\link{classifier_test}}.
#' @param parallel Logical, defaults to FALSE. Set to TRUE if you want to
#' parallellize the fitting.
#' @param BPPARAM object of class \code{bpparamClass} that specifies the
#'   back-end to be used for computations. See
#'   \code{bpparam} in \code{BiocParallel} package for details.
#' @param args_classifier arguments passed to the classifier test. See \code{\link{classifier_test}}.
#' @param args_mmd arguments passed to the mmd test. See \code{\link{mmd_test}}.
#' @param args_wass arguments passed to the wasserstein permutation test. See
#' \code{\link{wasserstein_permut}}.
#' @param nmax How many samples to use to compute the mmd test. See details.
#' @param rep How many permutations to run. Default to 50.
#' @param distinct_samples The samples to which each cell belong to. Only use
#' with method `distinct`. See `\code{\link{distinct_test}}` for help.
#' @return
#' A list containing the following components:
#' \itemize{
#'   \item *method* The method used to test
#'   \item *thresh* The threshold (if relevant)
#'   \item *statistic* the value of the test statistic.
#'   \item *p.value* the p-value of the test.
#' }
#' @examples
#' data('slingshotExample', package = "slingshot")
#' rd <- slingshotExample$rd
#' cl <- slingshotExample$cl
#' condition <- factor(rep(c('A','B'), length.out = nrow(rd)))
#' condition[110:139] <- 'A'
#' sds <- slingshot::getLineages(rd, cl)
#' topologyTest(sds, condition, rep = 10)
#' @details If there is only two conditions, default to `KS_mean`. Otherwise,
#' uses a classifier.
#'
#' More than one method can be specified at once, which avoids running slingshot on
#' the permutations more than once (as it is the slowest part).
#'
#' For the `mmd_test`, if `null=unbiased`, it is recommand to set `nmax=2000` or
#' something of that order of magnitude to avoid overflowing the memory.
#' @export
#' @importFrom Ecume classifier_test ks_test wasserstein_permut
#' @importFrom slingshot SlingshotDataSet getCurves slingPseudotime slingCurveWeights
#' @importFrom dplyr n_distinct bind_rows
#' @importFrom pbapply pblapply
#' @importFrom distinct distinct_test
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom stats model.matrix
#' @rdname topologyTest
setMethod(f = "topologyTest",
          signature = c(sds = "SlingshotDataSet"),
          definition = function(sds,
                                conditions,
                                rep = 100,
                                threshs = .01,
                                methods = ifelse(dplyr::n_distinct(conditions) == 2,
                                                 "KS_mean", "Classifier"),
                                parallel = FALSE,
                                BPPARAM = BiocParallel::bpparam(),
                                args_mmd = list(),
                                args_classifier = list(),
                                args_wass = list(),
                                nmax = nrow(slingshot::slingPseudotime(sds)),
                                distinct_samples = NULL){
            if (n_distinct(conditions) > 2 && methods != "Classifier") {
              warning("Changing to methods `classifier` since more than ",
                      "two conditions are present.")
              methods <- "Classifier"
            }
            res <- .topologyTest(sds = sds,
                                 conditions = conditions,
                                 rep = rep,
                                 threshs = threshs,
                                 methods = methods,
                                 parallel = parallel,
                                 BPPARAM = BPPARAM,
                                 args_mmd = args_mmd,
                                 args_wass = args_wass,
                                 args_classifier = args_classifier,
                                 nmax = nmax,
                                 distinct_samples = distinct_samples)
            return(res)
          }
)


#' @export
#' @rdname topologyTest
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData SummarizedExperiment
setMethod(f = "topologyTest",
          signature = c(sds = "SingleCellExperiment"),
          definition = function(sds,
                                conditions,
                                rep = 100,
                                threshs = .01,
                                methods = ifelse(dplyr::n_distinct(conditions) == 2,
                                                 "KS_mean", "Classifier"),
                                parallel = FALSE,
                                BPPARAM = BiocParallel::bpparam(),
                                args_mmd = list(),
                                args_classifier = list(),
                                args_wass = list(),
                                nmax = ncol(sds),
                                distinct_samples = NULL){
            if (is.null(sds@int_metadata$slingshot) &
                is.null(colData(sds)$slingshot)) {
              stop("For now this only works downstream of slingshot")
            }
            if (length(conditions) == 1) {
              if (conditions %in% colnames(SummarizedExperiment::colData(sds))) {
                conditions <- SummarizedExperiment::colData(sds)[, conditions]
              } else {
                stop("conditions is not a column of colData(sds)")
              }
            }
            return(topologyTest(slingshot::SlingshotDataSet(sds),
                                conditions = conditions,
                                rep = rep,
                                threshs = threshs,
                                methods = methods,
                                parallel = parallel,
                                BPPARAM = BPPARAM,
                                args_mmd = args_mmd,
                                args_wass = args_wass,
                                args_classifier = args_classifier,
                                nmax = nmax,
                                distinct_samples = distinct_samples))
          }
)

#' @rdname topologyTest
#' @importClassesFrom TrajectoryUtils PseudotimeOrdering
#' @export
setMethod(f = "topologyTest",
          signature = c(sds = "PseudotimeOrdering"),
          definition = function(sds,
                                conditions,
                                rep = 100,
                                threshs = .01,
                                methods = ifelse(dplyr::n_distinct(conditions) == 2,
                                                 "KS_mean", "Classifier"),
                                parallel = FALSE,
                                BPPARAM = BiocParallel::bpparam(),
                                args_mmd = list(),
                                args_classifier = list(),
                                args_wass = list(),
                                nmax = nrow(slingshot::slingPseudotime(sds)),
                                distinct_samples = NULL){
            if (n_distinct(conditions) > 2 && methods != "Classifier") {
              warning("Changing to methods `classifier` since more than ",
                      "two conditions are present.")
              methods <- "Classifier"
            }
            res <- .topologyTest(sds = sds,
                                 conditions = conditions,
                                 rep = rep,
                                 threshs = threshs,
                                 methods = methods,
                                 parallel = parallel,
                                 BPPARAM = BPPARAM,
                                 args_mmd = args_mmd,
                                 args_wass = args_wass,
                                 args_classifier = args_classifier,
                                 nmax = nmax,
                                 distinct_samples = distinct_samples)
            return(res)
          }
)
