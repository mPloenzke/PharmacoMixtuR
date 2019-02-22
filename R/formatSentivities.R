#' Wrapper to return pre-computed sensitivity profile
#'
#' A wrapper to return the sensitivity information for the
#' intersected PSets across cells and drugs. The sensitivity measure must be present
#' in the \code{PSet} to utilize this function (e.g. \code{'auc_recomputed'}).
#'
#' @param common List of \code{PSets}
#' @param sensitivity_measure Sensitivity measure to return
#'
#' @return Tibble containing sensitivity information
#'
#' @export
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @seealso \code{\link{intersectPSetWrapper}}
#' @keywords intersectPSet Sensitivity
#'
#' @importFrom PharmacoGx summarizeSensitivityProfiles
#' @importFrom dplyr mutate arrange bind_rows
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
formatSensitivities <- function(common, sensitivity_measure) {
  sensitivity.list <- lapply(1:length(common), function(co) {
    sensitivity.df <- PharmacoGx::summarizeSensitivityProfiles(pSet=common[[co]],
                                                               sensitivity.measure=sensitivity_measure,
                                                               summary.stat="median",
                                                               verbose=FALSE)
    drug.names <- row.names(sensitivity.df)
    sensitivity.df <- as.data.frame(sensitivity.df) %>%
      mutate(drug=drug.names) %>%
      gather(cell,value,-c(drug),na.rm=TRUE) %>%
      arrange(drug,cell) %>%
      mutate(experiment = names(common)[co],
             measure = sensitivity_measure)
  })
  rez <- do.call(bind_rows, sensitivity.list)
  return(rez)
}
