#' Wrapper to calculate sensitivity profile
#'
#' A wrapper to return the sensitivity information for the
#' intersected PSets across cells and drugs.
#'
#' @param sensitivitiy.tibb Tibble containing sensitivity information, drugs, and cells
#' @param drug.names List of drugs to fit model to 
#' @param prior_proportion Proportion of cells believed to be sensitive, per drug, to use as a prior
#' @param sd_na_val Avoid collapsing SD by setting poserior SD to this value
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
#' @importFrom dplyr summarise mutate group_by ungroup select filter pull arrange bind_rows
#' @importFrom mixtools normalmixEM
calcMixture <- function(sensitivity.tibb, drug.names=NULL, prior_proportion=.8, sd_na_val=.02) {
  if (is.null(drug.names)) {drug.names <- unique(sensitivity.tibb$drug)}
  sens.new <- params.new <- NULL
  for (dr in drug.names) {
    sens.short <- sensitivity.tibb %>% filter(drug==dr)
    params <- sens.short %>%
      group_by(experiment, measure) %>%
      mutate(Sensitive=ifelse(as.logical(value>=quantile(.$value,probs=prior_proportion)),'Sensitive','Resistant')) %>%
      group_by(Sensitive,experiment, measure) %>%
      summarise(mean=mean(value),sd=sd(value),count=n()) %>%
      ungroup() %>%
      mutate(Sensitive = ifelse(Sensitive=='Sensitive',1,0))
    if (any(is.na(params$sd))) {
      params <- params %>%
        mutate(sd = ifelse(is.na(sd),sd_na_val,sd))
    }
    for (stud in unique(sens.short$experiment)) {
      for(meas in unique(sens.short$measure)) {
        sens.temp <- sens.short %>%
          filter(experiment == stud,
                 measure == meas) %>%
          select(value) %>%
          as.matrix()
        param.means <- params %>%
          filter(experiment == stud,
                 measure == meas) %>%
          pull(mean)
        param.sd <- params %>%
          filter(experiment == stud,
                 measure == meas) %>%
          pull(sd)
        fit <- try(normalmixEM(sens.temp,
                           lambda = c(prior_proportion,1-prior_proportion),
                           mu = param.means,
                           sigma = param.sd),
                   silent=TRUE)
        if (class(fit) == 'try-error') {
          fit <- list(posterior = matrix(c(0,0),ncol=2),
                      mu = NA,
                      sigma = NA,
                      lambda = NA,
                      loglik = NA)
        }
        sens.new <- sens.short %>%
          filter(experiment == stud,
                 measure == meas) %>%
          mutate(posterior_Z0 = fit$posterior[,1],
                 posterior_Z1 = fit$posterior[,2],
                 posterior_Zmax = apply(fit$posterior,1,which.max)-1) %>%
          bind_rows(sens.new)
        params.new <- params %>%
          filter(experiment == stud,
                 measure == meas) %>%
          arrange(Sensitive) %>%
          mutate(posterior_mu = fit$mu,
                 posterior_sigma = fit$sigma,
                 posterior_lambda = fit$lambda,
                 loglik = fit$loglik,
                 drug = dr) %>%
          bind_rows(params.new)
      }
    }
  }
  return(list(sens.new,params.new))
}
