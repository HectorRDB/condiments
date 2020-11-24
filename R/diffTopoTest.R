.condition_sling <- function(sds, conditions) {
  psts <- NULL
  ws <- NULL
  for (cond in unique(conditions)) {
    sds_cond <- sds
    sds_cond@reducedDim <- sds_cond@reducedDim[conditions == cond, ]
    sds_cond@clusterLabels <- sds_cond@clusterLabels[conditions == cond, ]
    # TODO: handle cases where one cluster contains only one condition!!
    # missing_clusters <- which(colSums(sds_cond@clusterLabels) == 0)
    # sds_cond@clusterLabels <- sds_cond@clusterLabels[, -missing_clusters]
    # sds_cond@adjacency <- sds_cond@adjacency[-missing_clusters, -missing_clusters]
    # sds_cond@lineages <- lapply(sds_cond@lineages, function(l) {
    #   return(l[!l %in% names(missing_clusters)])
    # })
    sds_cond <- slingshot::getCurves(sds_cond, approx_points = 100)
    pst_cond <- slingshot::slingPseudotime(sds_cond)
    psts <- rbind(psts, pst_cond)
    w_conds <- slingshot::slingCurveWeights(sds_cond)
    ws <- rbind(ws, w_conds)
  }
  ws <- sweep(ws, 1, FUN = "/", STATS = apply(ws, 1, sum))
  return(list("psts" = psts, "ws" = ws))
}

.diffTopoTest_ks_all <- function(permutations, og, thresh) {
  psts <- lapply(permutations, '[[', 1) %>% do.call(what = 'rbind') %>%
    as.vector()
  ws <- lapply(permutations, '[[', 2) %>% do.call(what = 'rbind') %>%
    as.vector()
  og_psts <- og$psts %>% as.vector()
  og_ws <- og$ws %>% as.vector()
  res <- ks_test(x = og_psts, w_x = og_ws,
                 y = psts, w_y = ws,
                 thresh = thresh)
  return(res)
}

.diffTopoTest_ks_mean <- function(permutations, og, thresh, sds, rep) {
  psts <- lapply(permutations, '[[', 1) %>%
    lapply(function(df) {
      as.matrix(df[rownames(reducedDim(sds)), ])
    }) %>%
    Reduce(f = '+') %>%
    as.vector()
  psts <- psts / rep
  og_psts <- og$psts %>% as.vector()
  ws <- lapply(permutations, '[[', 2) %>%
    lapply(function(df) {
      as.matrix(df[rownames(reducedDim(sds)), ])
    }) %>%
    Reduce(f = '+') %>%
    as.vector()
  ws <- ws / rep
  og_ws <- og$ws %>% as.vector()
  res <- Ecume::ks_test(x = og_psts, w_x = og_ws,
                        y = psts, w_y = ws,
                        thresh = thresh)
  return(res)
}

.diffTopoTest_classifier <- function(permutations, og, thresh, sds, rep, ...){
  psts <- lapply(permutations, '[[', 1) %>%
    lapply(function(df) {
      as.matrix(df[rownames(reducedDim(sds)), ])
    }) %>%
    Reduce(f = '+')
  psts <- psts / rep
  colnames(psts) <- colnames(og$psts)
  res <- Ecume::classifier_test(x = og$psts, y = psts, thresh = thresh, ...)
  return(res)
}

.diffTopoTest_mmd2 <- function(permutations, og, sds, rep, ...){
  psts <- lapply(permutations, '[[', 1) %>%
    lapply(function(df) {
      as.matrix(df[rownames(reducedDim(sds)), ])
    }) %>%
    Reduce(f = '+')
  psts <- psts / rep
  colnames(psts) <- colnames(og$psts)
  frac <- 10^5 / (nrow(psts) * (nrow(psts) - 1))
  res <- Ecume::mmd_test(x = og$psts, y = psts, frac = frac, ...)
  return(res)
}

.diffTopoTest <- function(sds, conditions, rep = 200, thresh = .05,
                          method = "KS_mean", ...) {
  og <- .condition_sling(sds, conditions)
  permutations <- pbapply::pblapply(seq_len(rep), function(i) {
    condition <- sample(conditions)
    return(.condition_sling(sds, condition))
  })
  if (method == "KS_all") {
    res <- .diffTopoTest_ks_all(permutations, og, thresh)
  } else if (method == "KS_mean") {
    res <- .diffTopoTest_ks_mean(permutations, og, thresh, sds, rep)
  } else if (method == "Classifier") {
    res <- .diffTopoTest_classifier(permutations, og, thresh, sds, rep, ...)
  } else if (method == "mmd2") {
    res <- .diffTopoTest_mmd2(permutations, og, sds, rep, ...)
  }
  return(res[c("statistic", "p.value")])
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
#' @param rep How many permutations to run. Default to 50.
#' @param thresh the threshold for the KS test. See \code{\link{ks_test}}.
#' @param method The method to use to test. One of 'KS_mean', "KS_all', "mmd2'
#' and 'Classifier'. See details.
#' @param ... Other arguments passed to the test method. See details
#' @return
#' A list containing the following components:
#' \itemize{
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
#' diffTopoTest(sds, condition, rep = 20)
#' @details If there is only two conditions, default to `KS_mean`. Otherwise,
#' uses a classifier.
#' @export
#' @importFrom Ecume classifier_test ks_test
#' @importFrom slingshot SlingshotDataSet getCurves slingPseudotime slingCurveWeights
#' @importFrom dplyr n_distinct
#' @importFrom pbapply pblapply
#' @rdname diffTopoTest
setMethod(f = "diffTopoTest",
          signature = c(sds = "SlingshotDataSet"),
          definition = function(sds, conditions, rep = 200, thresh = .05,
    method = ifelse(dplyr::n_distinct(conditions) == 2, "KS_mean", "Classifier"),
    ...){
            if (n_distinct(conditions) > 2 && method != "classifier") {
              warning(paste0("Changing to method classifier since more than ",
                             "two conditions are present."))
              method <- "classifier"
            }
            res <- .diffTopoTest(sds = sds, conditions = conditions, rep = rep,
                                 thresh = thresh, method = method, ...)
            return(res)
          }
)


#' @export
#' @rdname diffTopoTest
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData
setMethod(f = "diffTopoTest",
          signature = c(sds = "SingleCellExperiment"),
          definition = function(sds, conditions, rep = 200, thresh = .05,
    method = ifelse(dplyr::n_distinct(conditions) == 2, "KS_mean", "Classifier"),
    ...){
            if (is.null(sds@int_metadata$slingshot)) {
              stop("For now this only works downstream of slingshot")
            }
            if (length(conditions) == 1) {
              if (conditions %in% colnames(SummarizedExperiment::colData(sds))) {
                conditions <- SummarizedExperiment::colData(sds)[, conditions]
              } else {
                stop("conditions is not a column of colData(sds)")
              }
            }
            return(diffTopoTest(slingshot::SlingshotDataSet(sds),
                                conditions = conditions,
                                rep = rep,
                                thresh = thresh,
                                method = method,
                                ...))
          }
)
