#' Wrapper to intersect PSets
#'
#' A wrapper to return the intersected PSets for cells and drugs
#' attaining a minimum number of intersections.
#'
#' @param datasets Vector of datasets to calculate intersections among
#' @param minimum_intersection Number of experiments to calculate intersection across
#' @param intersectOn Vetor of items to perform intersection across. Default is \code{c'cell.lines','drugs')}
#' @param sensitivity_measure Vector of sensitivity measures to return
#' @param remove_noisy_curves Remove noisy dose-response curves
#' @param remove_mislabelled_cells Remove incorrectly-labelled cells
#' @param mislabelled_cells Vector of cells to remove. If \code{NULL}, GDSC and CCLE mislabelled cells used.
#'
#' @return List of tibbles containing sensitivity information
#'
#' @export
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @seealso \code{\link{formatSensitivities}}
#' @keywords intersectPSet PSet
#'
#' @importFrom PharmacoGx filterNoisyCurves intersectPSet sensitivityMeasures
#' @importFrom tibble as.tibble
intersectPSetWrapper <- function(datasets,
                                 minimum_intersection = 1,
                                 intersectOn = c("cell.lines", "drugs"),
                                 sensitivity_measure = 'auc_recomputed',
                                 remove_noisy_curves = FALSE,
                                 remove_mislabelled_cells = FALSE,
                                 mislabelled_cells = NULL) {
  num_datasets <- length(datasets)
  combos <- t(combn(num_datasets,minimum_intersection))
  res <- lapply(1:nrow(combos), function(ro) {
    curr_list <- 'list('
    for (co in 1:ncol(combos)) {
      if (co>1) {curr_list <- paste(curr_list,',',sep='')}
      curr_list <- paste(curr_list,'"',datasets[combos[ro,co]],'"=',datasets[combos[ro,co]],sep='')
    }
    curr_list <- eval(parse(text=paste(curr_list,')',sep='')))
    common <- PharmacoGx::intersectPSet(curr_list,intersectOn=intersectOn, strictIntersect=TRUE)
    if (remove_mislabelled_cells) {
      str <- 'cells <- PharmacoGx::intersectList(unique(sensitivityInfo(common[[1]])$cellid)'
      if (length(common) > 1) {
        for (co in 2:length(common)) {
          str <- paste(str,', unique(PharmacoGx::sensitivityInfo(common[[',co,']])$cellid)',sep='')
        }
      }
      str <- paste(str,')',sep='')
      eval(parse(text=str))
      if (is.null(mislabelled_cells)) {
        mislabelled_cells <- c("LC-1F","HCC1937","MDA-MB-468","HuH-7","SW403","COR-L51","MOG-G-CCM","NB4")
        print('No mislabelled cells provided, using the following cells:')
        print(mislabelled_cells)
      }
      cells <- setdiff(cells, mislabelled_cells)
      common <- PharmacoGx::intersectPSet(pSets = curr_list, intersectOn = intersectOn, cells=cells)
    }
    if (remove_noisy_curves) {
      for (li in 1:length(common)) {
        ps.filter <- PharmacoGx:::filterNoisyCurves(common[[li]], nthread=detectCores())
        common[[li]]@sensitivity$info <- common[[li]]@sensitivity$info[ps.filter$ok, ]
        common[[li]]@sensitivity$raw <- common[[li]]@sensitivity$raw[ps.filter$ok, , ]
        common[[li]]@sensitivity$profiles <- common[[li]]@sensitivity$profiles[ps.filter$ok, ]
      }
    }
    meas1 <- sensitivity_measure[sensitivity_measure %in% PharmacoGx::sensitivityMeasures(common[[1]])]
    meas2 <- sensitivity_measure[!(sensitivity_measure %in% meas1)]
    if (any(!(meas2 %in% c('ec50', 'e_inf', 'hill_slope','raw_auc')))) {
      print(paste('Unknown sensitivity measure:',
                  meas2[!(meas2 %in% c('ec50', 'e_inf', 'hill_slope'))],
                  sep=' '))
    }
    tres <- lapply(sensitivity_measure, function(meas) {
      if (meas %in% meas1) {
        formatSensitivities(common, meas)
      } else if (meas %in% meas2) {
        computeSensitivities(common, meas)
      }
    })
    as.tibble(do.call(bind_rows, tres))
  })
  names(res) <- apply(combos, 1, function(ii) {
    paste(datasets[ii],collapse='_')
  })
  return(res)
}
