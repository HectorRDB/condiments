.differentiationTest <- function(ws, conditions, global = TRUE, pairwise = FALSE,
                                 method = "Classifier", thresh,
                                 classifier_method = "rf",
                                 args_mmd = list(), args_classifier = list(),
                                 args_wass = list()) {
  ws <- sweep(ws, 1, FUN = "/", STATS = apply(ws, 1, sum))
  pairs <- utils::combn(ncol(ws), 2)
  if (ncol(ws) == 2) {
    global <- FALSE
    pairwise <- TRUE
  }
  n_conditions <- dplyr::n_distinct(conditions)
  pairwise_test <- apply(pairs, 2, function(pair) {
    ws_p <- ws[, pair]
    keep <- rowSums(ws_p) > 0
    ws_p <- sweep(ws_p[keep, ], 1, FUN = "/", STATS = apply(ws_p[keep, ], 1, sum))
    nmin <- min(table(conditions[keep]))
    xs <- lapply(unique(conditions), function(cond) {
      ws_cond <- ws_p[conditions[keep] == cond, ]
      ws_cond <- ws_cond[sample(seq_len(nrow(ws_cond)), nmin), ]
      return(as.matrix(ws_cond))
    })
    if (method == "Classifier") {
      args <- args_classifier
      args$method <- classifier_method
      args$x <- xs; args$thresh <- thresh
      return(do.call(Ecume::classifier_test, args))
    }
    if (method == "mmd") {
      n <- max(table(conditions))
      frac <- 10^5 / (n * (n - 1))
      args <- args_mmd
      args$x <- xs[[1]]; args$y <- xs[[2]]; args$frac <- frac
      return(do.call(Ecume::mmd_test, args))
    }
    if (method == "wasserstein_permutation") {
      n <- max(table(conditions))
      S <- min(10^5, n)
      args <- args_wass
      args$x <- xs[[1]]; args$y <- xs[[2]]
      args$S <- S; args$fast <- TRUE
      return(do.call(Ecume::wasserstein_permut, args))
    }
  }) %>%
    dplyr::bind_rows(.id = "pair") %>%
    dplyr::mutate(pair = as.numeric(pair),
                  pair = paste0(pairs[1, pair], "vs", pairs[2, pair])) %>%
    dplyr::select(pair, statistic, p.value)
  nmin <- min(table(conditions))
  xs <- lapply(unique(conditions), function(cond) {
    ws_cond <- ws[conditions == cond, ]
    ws_cond <- ws_cond[sample(seq_len(nrow(ws_cond)), nmin), ]
  })
  if (method == "Classifier") {
    args <- args_classifier
    args$method <- classifier_method
    args$x <- xs; args$thresh <- thresh
    glob_test <- do.call(Ecume::classifier_test, args)
  }
  if (method == "mmd") {
    n <- max(table(conditions))
    frac <- 10^5 / (n * (n - 1))
    args <- args_mmd
    args$x <- xs[[1]]; args$y <- xs[[2]]; args$frac <- frac
    glob_test <- do.call(Ecume::mmd_test, args)
  }
  if (method == "wasserstein_permutation") {
    n <- max(table(conditions))
    S <- min(10^5, n)
    args <- args_wass
    args$x <- xs[[1]]; args$y <- xs[[2]]; args$S <- S; args$fast <- TRUE
    glob_test <- do.call(Ecume::wasserstein_permut, args)
  }
  glob_test <- data.frame("pair" = "All",
                          "statistic" = glob_test$statistic,
                          "p.value" = glob_test$p.value)

  if (global == TRUE & pairwise == FALSE) return(glob_test)
  if (global == FALSE & pairwise == TRUE) return(pairwise_test)
  if (global == TRUE & pairwise == TRUE) {
    return(dplyr::bind_rows(glob_test, pairwise_test))
  }
}


#' Differential Differentiation Test
#'
#' @description Test whether or not the cell repartition between lineages is
#' independent of the conditions
#'
#' @param cellWeights Can be either a \code{\link{SlingshotDataSet}}, a
#' \code{\link{SingleCellExperiment}} object or a matrix of cell weights
#' defining the probability that a cell belongs to a particular lineage.
#' Each row represents a cell and each column represents a lineage. If only a
#' single lineage, provide a matrix with one column containing all values of 1.
#' @param conditions Either the vector of conditions, or a character indicating which
#' column of the metadata contains this vector
#' @param global If TRUE, test for all pairs simultaneously.
#' @param pairwise If TRUE, test for all pairs independently.
#' @param method One of "Classifier" or "mmd".
#' @param thresh The threshold for the classifier test. See details.
#' Default to .05.
#' @param args_mmd arguments passed to the mmd test. See \code{\link{mmd_test}}.
#' @param args_wass arguments passed to the wasserstein permutation test. See
#' \code{\link{wasserstein_permut}}.
#' @param args_classifier arguments passed to the classifier test. See \code{\link{classifier_test}}.
#' @param classifier_method The method used in the classifier test. Default to
#' "rf", i.e random forest.
#' @importFrom slingshot slingshot SlingshotDataSet
#' @importFrom utils combn
#' @importFrom dplyr n_distinct
#' @importFrom Ecume classifier_test mmd_test wasserstein_permut
#' @importFrom matrixStats rowMaxs rowMins
#' @return A data frame with 3 columns:
#' \itemize{
#'   \item *pair* for individual pairs, the lineages numbers. For global,
#'   \code{"All"}.
#'   \item *p.value* the pvalue for the test at the global or pair level
#'   \item *statistic* The classifier accuracy
#' }
#' @examples
#' data('slingshotExample', package = "slingshot")
#' rd <- slingshotExample$rd
#' cl <- slingshotExample$cl
#' condition <- factor(rep(c('A','B'), length.out = nrow(rd)))
#' condition[110:139] <- 'A'
#' sds <- slingshot::slingshot(rd, cl)
#' differentiationTest(sds, condition)
#' @export
#' @rdname differentiationTest
setMethod(f = "differentiationTest",
          signature = c(cellWeights = "matrix"),
          definition = function(cellWeights, conditions, global = TRUE,
                                pairwise = FALSE,
                                method = c("Classifier", "mmd", "wasserstein_permutation"),
                                classifier_method = "rf",
                                thresh = .01, args_classifier = list(),
                                args_mmd = list(), args_wass = list()){
            method <- match.arg(method)
            if (ncol(cellWeights) == 1) {
              stop("This only works with more than one lineage.")
            }
            if (dplyr::n_distinct(conditions) > 2 && method != "Classifier") {
              method <- "Classifier"
              warning("If more than two conditions are present, ",
                      "only the Classifier method is possible.")
            }
            res <- .differentiationTest(ws = cellWeights, conditions = conditions,
                                        global = global, pairwise = pairwise,
                                        method = method, thresh = thresh,
                                        classifier_method = classifier_method,
                                        args_mmd = args_mmd, args_wass = args_wass,
                                        args_classifier = args_classifier)
            return(res)
          }
)


#' @rdname differentiationTest
#' @importFrom slingshot as.PseudotimeOrdering
setMethod(f = "differentiationTest",
          signature = c(cellWeights = "SlingshotDataSet"),
          definition = function(cellWeights, conditions, global = TRUE,
                                pairwise = FALSE,
                                method = c("Classifier", "mmd", "wasserstein_permutation"),
                                classifier_method = "rf",
                                thresh = .01, args_classifier = list(),
                                args_mmd = list(), args_wass = list()){
            method <- match.arg(method)
            if (nLineages(cellWeights) == 1) {
              stop("This only works with more than one lineage.")
            }
            if (dplyr::n_distinct(conditions) > 2 && method != "Classifier") {
              method <- "Classifier"
              warning("If more than two conditions are present, ",
                      "only the Classifier method is possible.")
            }
            res <- differentiationTest(cellWeights = as.PseudotimeOrdering(cellWeights),
                                       conditions = conditions,
                                       global = global, pairwise = pairwise,
                                       method = method, thresh = thresh,
                                       classifier_method = classifier_method,
                                       args_mmd = args_mmd, args_wass = args_wass,
                                       args_classifier = args_classifier)
            return(res)
          }
)


#' @export
#' @rdname differentiationTest
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData
setMethod(f = "differentiationTest",
          signature = c(cellWeights = "SingleCellExperiment"),
          definition = function(cellWeights, conditions, global = TRUE,
                                pairwise = FALSE,
                                method = c("Classifier", "mmd", "wasserstein_permutation"),
                                classifier_method = "rf",
                                thresh = .01, args_classifier = list(),
                                args_mmd = list(), args_wass = list()){
            if (is.null(cellWeights@int_metadata$slingshot)) {
              stop("For now this only works downstream of slingshot")
            }
            if (length(conditions) == 1) {
              if (conditions %in%
                  colnames(SummarizedExperiment::colData(cellWeights))) {
                conditions <- SummarizedExperiment::colData(cellWeights)[, conditions]
              } else {
                stop("conditions is not a column of colData(cellWeights)")
              }
            }
            return(differentiationTest(slingshot::SlingshotDataSet(cellWeights),
                                       conditions = conditions, global = global,
                                       pairwise = pairwise, method = method,
                                       classifier_method = classifier_method,
                                       thresh = thresh, args_mmd = args_mmd,
                                       args_wass = args_wass,
                                       args_classifier = args_classifier))
          }
)


#' @rdname differentiationTest
#' @importClassesFrom TrajectoryUtils PseudotimeOrdering
#' @export
setMethod(f = "differentiationTest",
          signature = c(cellWeights = "PseudotimeOrdering"),
          definition = function(cellWeights, conditions, global = TRUE,
                                pairwise = FALSE,
                                method = c("Classifier", "mmd", "wasserstein_permutation"),
                                classifier_method = "rf",
                                thresh = .01, args_classifier = list(),
                                args_mmd = list(), args_wass = list()){
            method <- match.arg(method)
            if (nLineages(cellWeights) == 1) {
              stop("This only works with more than one lineage.")
            }
            if (dplyr::n_distinct(conditions) > 2 && method != "Classifier") {
              method <- "Classifier"
              warning("If more than two conditions are present, ",
                      "only the Classifier method is possible.")
            }
            if (slingParams(cellWeights)$reweight | slingParams(cellWeights)$reassign) {
              ws <- slingshot::slingCurveWeights(cellWeights, as.probs = TRUE)
            } else {
              ws <- .sling_reassign(cellWeights)
            }

            res <- .differentiationTest(ws = ws, conditions = conditions,
                                        global = global, pairwise = pairwise,
                                        method = method, thresh = thresh,
                                        classifier_method = classifier_method,
                                        args_mmd = args_mmd, args_wass = args_wass,
                                        args_classifier = args_classifier)
            return(res)
          }
)
