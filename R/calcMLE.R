#' Wrapper to calculate sensitivity profile
#'
#' A wrapper to return the sensitivity information for the
#' intersected PSets across cells and drugs.
#'
#' @param sensitivitiy.tibb
#' @param drug.names
#' @param fit
#'
#' @return List of tibbles; one containing posterior estimates and one containing estimated parameters
#'
#' @export
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @seealso \code{\link{intersectPSetWrapper}}
#' @keywords intersectPSet Sensitivity
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr summarise mutate group_by ungroup select filter
#' @importFrom Rfast beta.mle
calcMLE <- function(sensitivity.tibb, drug.names=NULL, fit='normal') {
  if (is.null(drug.names)) {drug.names <- unique(sensitivity.tibb$drug)}
  if (fit=='beta') {
    params <- sensitivity.tibb %>%
      filter(drug %in% drug.names) %>%
      group_by(experiment, measure, drug) %>%
      mutate(value = ifelse(value==0,1e-6,value)) %>%
      do(param1 = beta.mle(.$value)$param[1],
         param2 = beta.mle(.$value)$param[2]) %>%
      ungroup() %>%
      mutate(param1 = as.vector(unlist(.$param1)),
             param2 = as.vector(unlist(.$param2)))
    sens.new <- sensitivity.tibb %>%
      filter(drug %in% drug.names) %>%
      left_join(params, by=c('experiment','measure','drug')) %>%
      mutate(prob = pbeta(.$value,param1,param2),
             Sensitive=ifelse(as.logical((prob<=.025) | (prob>=.975)),1,0)) %>%
      select(-param1,-param2)
  } else if (fit=='normal') {
    params <- sensitivity.tibb %>%
      filter(drug %in% drug.names) %>%
      group_by(experiment, measure, drug) %>%
      summarise(param1=mean(value),param2=sd(value),count=n()) %>%
      ungroup()
    sens.new <- sensitivity.tibb %>%
      filter(drug %in% drug.names) %>%
      left_join(params, by=c('experiment','measure','drug')) %>%
      mutate(prob = pnorm(.$value,param1,param2),
             Sensitive=ifelse(as.logical((prob<=.025) | (prob>=.975)),1,0)) %>%
      select(-param1,-param2,-count)
  }
  return(list(sens.new,params))
}
