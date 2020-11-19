.diffDifferentiationTest <- function(sds, conditions, global = TRUE,
                                     pairwise = FALSE, method = "Classifier",
                                     thresh, ...) {
  if (method != "Classifier") stop("Only Classifier is possible for now")
  pairs <- utils::combn(length(slingshot::slingLineages(sds)), 2)
  n_conditions <- dplyr::n_distinct(conditions)
  nmin <- min(table(conditions))
  pairwise_test <- apply(pairs, 2, function(pair) {
    if(method == "Classifier") {
      ws <- slingshot::slingCurveWeights(sds, as.probs = TRUE)[, pair]
      ws <- sweep(ws, 1, FUN = "/", STATS = apply(ws, 1, sum))
      xs <- lapply(unique(conditions), function(cond) {
        ws_cond <- ws[conditions == cond, ]
        ws_cond <- ws_cond[sample(seq_len(nrow(ws_cond)), nmin), ]
      })
      return(Ecume::classifier_test(x = xs, thresh = thresh, ...))
    }
  }) %>%
    dplyr::bind_rows(.id = "pair") %>%
    dplyr::mutate(pair = as.character(pair)) %>%
    dplyr::select(pair, statistic, p.value)
  if (method == "Classifier") {
    ws <- slingshot::slingCurveWeights(sds, as.probs = TRUE)
    xs <- lapply(unique(conditions), function(cond) {
      ws_cond <- ws[conditions == cond, ]
      ws_cond <- ws_cond[sample(seq_len(nrow(ws_cond)), nmin), ]
    })
    glob_test <- Ecume::classifier_test(xs, thresh = thresh, ...)
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
#' @param sds The final object after running slingshot. Can be either a
#' \code{\link{SlingshotDataSet}} or a \code{\link{SingleCellExperiment}} object.
#' @param conditions Either the vector of conditions, or a character indicating which
#' column of the metadata contains this vector
#' @param global If TRUE, test for all pairs simultaneously.
#' @param pairwise If TRUE, test for all pairs independently.
#' @param method For now, only "Classifier" is accepted.
#' @param thresh The threshold for the classifier test. See details.
#' Default to .05.
#' @param ... Other arguments passed to \link[Ecume]{classifier_test}.
#' @importFrom slingshot slingshot SlingshotDataSet
#' @importFrom utils combn
#' @importFrom dplyr n_distinct
#' @importFrom Ecume classifier_test
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
#' diffDifferentiationTest(sds, condition)
#' @export
#' @rdname diffDifferentiationTest
setMethod(f = "diffDifferentiationTest",
          signature = c(sds = "SlingshotDataSet"),
          definition = function(sds, conditions, global = TRUE, pairwise = FALSE,
                                method = "Classifier", thresh = .05, ...){
            res <- .diffDifferentiationTest(sds = sds,
                                            conditions = conditions,
                                            global = global,
                                            pairwise = pairwise,
                                            method = method,
                                            thresh = thresh,
                                            ...)
            return(res)
          }
)


#' @export
#' @rdname diffDifferentiationTest
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData
setMethod(f = "diffDifferentiationTest",
          signature = c(sds = "SingleCellExperiment"),
          definition = function(sds, conditions,  global = TRUE,
                                pairwise = FALSE, method = "Classifier",
                                thresh = .05, ...){
            if (is.null(sds@int_metadata$slingshot)) {
              stop("For now this only works downstream of slingshot")
            }
            if (length(conditions == 1)) {
              if (conditions %in% colnames(SummarizedExperiment::colData(sds))) {
                conditions <- SummarizedExperiment::colData(sds)[, conditions]
              } else {
                stop("conditions is not a column of colData(sds)")
              }
            }
            return(diffDifferentiationTest(slingshot::SlingshotDataSet(sds),
                                           conditions = conditions,
                                           global = global,
                                           pairwise = pairwise,
                                           method = method,
                                           thresh = thresh,
                                           ...))
          }
)
