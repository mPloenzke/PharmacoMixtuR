#' Wrapper to calculate and return a sensitivity profile
#'
#' A wrapper to return the sensitivity information for the
#' intersected PSets across cells and drugs. The sensitivity measure must not be present
#' in the \code{PSet} but calculable from the dose-response information
#' (e.g. \code{'ec50', 'e_inf', 'hill_slope','raw_auc'}).
#'
#' @param common List of \code{PSets}
#' @param sensitivity_measure Sensitivity measure to return
#' @param truncateViabilities Truncate \code{viabilities>100}
#'
#' @return Tibble containing sensitivity information
#'
#' @export
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @seealso \code{\link{intersectPSetWrapper}}
#' @keywords intersectPSet Sensitivity
#'
#' @importFrom PharmacoGx computeAUC logLogisticRegression
#' @importFrom dplyr mutate group_by bind_rows do ungroup case_when
#' @importFrom magrittr %>%
computeSensitivities <- function(common, sensitivity_measure,truncateViabilities=TRUE) {
  sensitivity.list <- lapply(1:length(common), function(co) {
    tdat <- formatPSetDoseResponse(common[[co]]) %>%
      mutate(Concentration = as.numeric(Concentration),
             Viability = as.numeric(Viability)) %>%
      group_by(drug,cellid) %>%
      do(value = case_when(sensitivity_measure == 'raw_auc' ~ computeAUC(concentration=.$Concentration,
                                                                         viability=.$Viability/100,
                                                                         trunc=truncateViabilities,
                                                                         conc_as_log=FALSE,
                                                                         viability_as_pct=FALSE,
                                                                         area.type='Actual',
                                                                         verbose=FALSE),
                           sensitivity_measure == 'ec50' ~ logLogisticRegression(conc=.$Concentration,
                                                                                 viability=.$Viability/100,
                                                                                 trunc=truncateViabilities,
                                                                                 conc_as_log = FALSE,
                                                                                 viability_as_pct=FALSE,
                                                                                 family='normal')$EC50,
                           sensitivity_measure == 'e_inf' ~ logLogisticRegression(conc=.$Concentration,
                                                                                  viability=.$Viability/100,
                                                                                  trunc=truncateViabilities,
                                                                                  conc_as_log = FALSE,
                                                                                  viability_as_pct=FALSE,
                                                                                  family='normal')$E_inf,
                           sensitivity_measure == 'hill_slope' ~ logLogisticRegression(conc=.$Concentration,
                                                                                       viability=.$Viability/100,
                                                                                       trunc=truncateViabilities,
                                                                                       conc_as_log = FALSE,
                                                                                       viability_as_pct=FALSE,
                                                                                       family='normal')$HS)) %>%
      mutate(value = unlist(value)) %>%
      ungroup() %>%
      rename(cell=cellid) %>%
      mutate(measure = sensitivity_measure,
             experiment = common[[co]]@annotation$name)

  })
  rez <- do.call(bind_rows, sensitivity.list)
  return(rez)
}
