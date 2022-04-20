.multiple_samples <- function(sds = NULL, pseudotime = NULL, cellWeights = NULL,
                              conditions, Test, Samples, ...) {
  samples <- unique(Samples)
  names(samples) <- samples
  res <- lapply(samples, function(sample) {
    if (Test == "fateSelectionTest") {
      return(fateSelectionTest(cellWeights[Samples == sample, ],
                               conditions = conditions[Samples == sample], ...))
    }
    if (Test == "progressionTest") {
      return(progressionTest(pseudotime[Samples == sample, ],
                             cellWeights[Samples == sample, ],
                             conditions = conditions[Samples == sample],
                             ...))
    }
  }) %>%
    bind_rows(.id = "Sample")
  overall_res <- res %>%
    dplyr::group_by(across(starts_with("lineage")),
                    across(starts_with("pair"))) %>%
    dplyr::summarise(statistic = Ecume::stouffer_zscore(
      pvals = p.value, weights = rep(1, length(samples)))$statistic,
                     p.value = Ecume::stouffer_zscore(
                       p.value, weights = rep(1, length(samples)))$p.value) %>%
    dplyr::mutate(Sample = "Combined")
  return(bind_rows(res, overall_res))
}

.topologyTest_multipleSamples <- function(sds, conditions, Samples, ...) {
  samples <- unique(Samples)
  names(samples) <- samples
  condition_samples <- interaction(unique(conditions), samples)
  curves <- .generate_permutations_curves(sds, condition_samples, ...)
  res <- lapply(samples, function(sample) {
    inds <- Samples == sample
    permutations <- lapply(curves$permutations, function(perm){
      lapply(perm, function(dat) dat[inds, ]) %>% return()
    })
    og <- lapply(curves$og, function(dat){
       return(dat[inds, ])
    })
    .topologyTest_all_selected(
      permutations = permutations, og = og, sds = sds,
      conditions = conditions, ...) %>%
      bind_rows(.id = "method") %>%
      return()
  }) %>%
    bind_rows(.id = "Sample")
  overall_res <- res %>%
    dplyr::group_by(across(starts_with("method")),
                    across(starts_with("thresh"))) %>%
    dplyr::summarise(statistic = Ecume::stouffer_zscore(
      pvals = p.value, weights = rep(1, length(samples)))$statistic,
      p.value = Ecume::stouffer_zscore(
        p.value, weights = rep(1, length(samples)))$p.value) %>%
    dplyr::mutate(Sample = "Combined")
  return(bind_rows(res, overall_res))
}

# fate selection test ----
#' Differential fate selection Test with multiple samples
#'
#' @description Test whether or not the cell repartition between lineages is
#' independent of the conditions, with samples not being confounded by conditions
#'
#' @param cellWeights Can be either a \code{\link{SlingshotDataSet}}, a
#' \code{\link{SingleCellExperiment}} object or a matrix of cell weights
#' defining the probability that a cell belongs to a particular lineage.
#' Each row represents a cell and each column represents a lineage. If only a
#' single lineage, provide a matrix with one column containing all values of 1.
#' @param conditions Either the vector of conditions, or a character indicating which
#' column of the metadata contains this vector.
#' @param Samples A vector assigning each cell to a sample. Samples must be shared across all conditions.
#' @param ... Other arguments passed to \code{\link{fateSelectionTest}}.
#' @return The same object has the \code{\link{fateSelectionTest}} with one more column per sample.
#' @importFrom dplyr group_by summarize mutate
#' @examples
#' data('slingshotExample', package = "slingshot")
#' rd <- slingshotExample$rd
#' cl <- slingshotExample$cl
#' condition <- factor(rep(c('A','B'), length.out = nrow(rd)))
#' condition[110:139] <- 'A'
#' sds <- slingshot::slingshot(rd, cl)
#' samples <- sample(1:2, 140, replace = TRUE)
#' fateSelectionTest_multipleSamples(cellWeights = sds, conditions = condition, Samples = samples)
#' @export
#' @rdname fateSelectionTest_multipleSamples
setMethod(f = "fateSelectionTest_multipleSamples",
          signature = c(cellWeights = "matrix"),
          definition = function(cellWeights, conditions, Samples, ...) {
            return(.multiple_samples(NULL, NULL, cellWeights, conditions,
                                     "fateSelectionTest", Samples, ...))
          })

#' @rdname fateSelectionTest_multipleSamples
#' @importFrom slingshot as.PseudotimeOrdering
setMethod(f = "fateSelectionTest_multipleSamples",
          signature = c(cellWeights = "SlingshotDataSet"),
          definition = function(cellWeights, conditions, Samples, ...){
            return(fateSelectionTest_multipleSamples(
              cellWeights = as.PseudotimeOrdering(cellWeights),
              conditions = condition, Samples = Samples, ...))
          })

#' @rdname fateSelectionTest_multipleSamples
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData
setMethod(f = "fateSelectionTest_multipleSamples",
          signature = c(cellWeights = "SingleCellExperiment"),
          definition = function(cellWeights, conditions, Samples, ...){
            if (is.null(cellWeights@int_metadata$slingshot) &
                is.null(colData(cellWeights)$slingshot)) {
              stop("For now this only works downstream of slingshot")
            }
            if (length(cellWeights) == 1) {
              if (conditions %in%
                  colnames(SummarizedExperiment::colData(cellWeights))
              ) {
                conditions <-
                  SummarizedExperiment::colData(cellWeights)[, conditions]
              } else {
                stop("conditions is not a column of colData(cellWeights)")
              }
            }
            return(fateSelectionTest_multipleSamples(
              slingshot::SlingshotDataSet(cellWeights),
              conditions = conditions, Samples = Samples, ...))
          })

#' @rdname fateSelectionTest_multipleSamples
#' @importClassesFrom TrajectoryUtils PseudotimeOrdering
#' @export
setMethod(f = "fateSelectionTest_multipleSamples",
          signature = c(cellWeights = "PseudotimeOrdering"),
          definition = function(cellWeights, conditions, Samples, ...){
            if (nLineages(cellWeights) == 1) {
              stop("This only works with more than one lineage.")
            }
            if (slingParams(cellWeights)$reweight | slingParams(cellWeights)$reassign) {
              ws <- slingshot::slingCurveWeights(cellWeights, as.probs = TRUE)
            } else {
              ws <- .sling_reassign(cellWeights)
            }
            return(.multiple_samples(NULL, NULL, ws, conditions,
                                     "fateSelectionTest", Samples, ...))
          })

# Progression test ----

#' Differential Progression Test with multiple samples
#'
#' @description Test whether or not the pseudotime distribution are identical
#' within lineages between conditions, with samples not being confounded by conditions
#'
#' @param pseudotime Can be either a \code{\link{SlingshotDataSet}} or a
#' \code{\link{SingleCellExperiment}} object or a matrix of pseudotime values,
#' each row represents a cell and each column represents a lineage.
#' @param cellWeights If `pseudotime` is a matrix of pseudotime values, this
#' represent the cell weights for each lineage. Ignored if `pseudotime` is not
#' a matrix.
#' @param conditions Either the vector of conditions, or a character indicating which
#' column of the metadata contains this vector.
#' @param Samples A vector assigning each cell to a sample. Samples must be shared across all conditions.
#' @param ... Other arguments passed to \code{\link{progressionTest}}.
#' @return The same object has the \code{\link{progressionTest}} with one more column per sample.
#' @examples
#' data('slingshotExample', package = "slingshot")
#' rd <- slingshotExample$rd
#' cl <- slingshotExample$cl
#' condition <- factor(rep(c('A','B'), length.out = nrow(rd)))
#' condition[110:139] <- 'A'
#' sds <- slingshot::slingshot(rd, cl)
#' samples <- sample(1:2, 140, replace = TRUE)
#' progressionTest_multipleSamples(pseudotime = sds, conditions = condition, Samples = samples)
#' @export
#' @rdname progressionTest_multipleSamples
setMethod(f = "progressionTest_multipleSamples",
          signature = c(pseudotime = "matrix"),
          definition = function(pseudotime, cellWeights, conditions, Samples, ...) {
            return(.multiple_samples(NULL, pseudotime, cellWeights, conditions,
                                     "progressionTest", Samples, ...))
})

#' @rdname progressionTest_multipleSamples
#' @importFrom slingshot as.PseudotimeOrdering
setMethod(f = "progressionTest_multipleSamples",
          signature = c(pseudotime = "SlingshotDataSet"),
          definition = function(pseudotime, conditions, Samples, ...){
  return(progressionTest_multipleSamples(
    pseudotime = as.PseudotimeOrdering(pseudotime),
    conditions = condition, Samples = Samples, ...))
})

#' @rdname progressionTest_multipleSamples
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData
setMethod(f = "progressionTest_multipleSamples",
          signature = c(pseudotime = "SingleCellExperiment"),
          definition = function(pseudotime, conditions, Samples, ...){
            if (is.null(pseudotime@int_metadata$slingshot) &
                is.null(colData(pseudotime)$slingshot)) {
              stop("For now this only works downstream of slingshot")
            }
            if (length(conditions) == 1) {
              if (conditions %in%
                  colnames(SummarizedExperiment::colData(pseudotime))
              ) {
                conditions <-
                  SummarizedExperiment::colData(pseudotime)[, conditions]
              } else {
                stop("conditions is not a column of colData(pseudotime)")
              }
            }
            return(progressionTest_multipleSamples(
              slingshot::SlingshotDataSet(pseudotime),
              conditions = conditions, Samples = Samples, ...))
})

#' @rdname progressionTest_multipleSamples
#' @importClassesFrom TrajectoryUtils PseudotimeOrdering
#' @export
setMethod(f = "progressionTest_multipleSamples",
          signature = c(pseudotime = "PseudotimeOrdering"),
          definition = function(pseudotime, conditions, Samples, ...){
            pst <- slingshot::slingPseudotime(pseudotime, na = FALSE)
            ws <- slingshot::slingCurveWeights(pseudotime, as.probs = TRUE)
            return(.multiple_samples(NULL, pst, ws, conditions,
                                     "progressionTest", Samples, ...))
          })


# Topology test ----

#' Differential Topology Test with multiple samples
#'
#' @description Test whether or not slingshot should be fitted independently
#' for different conditions or not, per sample, with samples not being confounded by conditions.
#' @param sds A slingshot object already run on the full dataset. Can be either a
#' \code{\link{SlingshotDataSet}} or a \code{\link{SingleCellExperiment}} object.
#' @param conditions Either the vector of conditions, or a character indicating which
#' column of the metadata contains this vector.
#' @param Samples A vector assigning each cell to a sample. Samples must be shared across all conditions.
#' @param ... Other arguments passed to \code{\link{topologyTest}}.
#' @import slingshot
#' @importFrom dplyr across starts_with
#' @return The same object has the \code{\link{topologyTest}} with one more column per sample.
#' @examples
#' data('slingshotExample', package = "slingshot")
#' rd <- slingshotExample$rd
#' cl <- slingshotExample$cl
#' condition <- factor(rep(c('A','B'), length.out = nrow(rd)))
#' condition[110:139] <- 'A'
#' sds <- slingshot::slingshot(rd, cl)
#' samples <- sample(1:2, 140, replace = TRUE)
#' topologyTest_multipleSamples(sds = sds, conditions = condition,
#'                              Samples = samples, rep = 10)
#' @export
#' @rdname topologyTest_multipleSamples
setMethod(f = "topologyTest_multipleSamples",
          signature = c(sds = "SlingshotDataSet"),
          definition = function(sds, conditions, Samples, ...){
            res <- .topologyTest_multipleSamples(sds = as.PseudotimeOrdering(sds),
                                                 conditions = conditions,
                                                 Samples, ...)
            return(res)
          }
)


#' @export
#' @rdname topologyTest_multipleSamples
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData
setMethod(f = "topologyTest_multipleSamples",
          signature = c(sds = "SingleCellExperiment"),
          definition = function(sds, conditions, Samples, ...){
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
                                conditions = conditions, Samples, ...))
          }
)

#' @rdname topologyTest_multipleSamples
#' @importClassesFrom TrajectoryUtils PseudotimeOrdering
#' @export
setMethod(f = "topologyTest_multipleSamples",
          signature = c(sds = "PseudotimeOrdering"),
          definition = function(sds, conditions, Samples, ...){
            return(
              .topologyTest_multipleSamples(sds = sds, conditions = conditions,
                                            Samples, ...))
          }
)
