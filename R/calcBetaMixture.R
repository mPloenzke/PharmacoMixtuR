#' Wrapper to calculate beta mixture model
#'
#' A wrapper to return the sensitivity information for the
#' intersected PSets across cells and drugs.
#'
#' @param sensitivitiy.tibb Tibble containing sensitivity information, drugs, and cells
#' @param drug.names List of drugs to fit model to 
#' @param prior_proportion Proportion of cells believed to be sensitive, per drug, to use as a prior
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
#' @importFrom dplyr summarise mutate group_by ungroup select filter arrange bind_rows
#' @importFrom betareg betamix clusters
#' @importFrom stats plogis
#' @importFrom Rfast beta.mle
calcBetaMixture <- function(sensitivity.tibb, drug.names=NULL, prior_proportion=.8) {
  if (is.null(drug.names)) {drug.names <- unique(sensitivity.tibb$drug)}
  sens.new <- params.new <- NULL
  for (dr in drug.names) {
    sens.short <- sensitivity.tibb %>% filter(drug==dr)
    params <- sens.short %>%
      mutate(value = ifelse(value==0,1e-3,value)) %>%
      group_by(experiment, measure) %>%
      mutate(Sensitive=ifelse(as.logical(value>=quantile(.$value,probs=prior_proportion)),'Sensitive','Resistant')) %>%
      group_by(Sensitive,experiment, measure) %>%
      do(param1 = try(beta.mle(.$value)$param[1],silent=TRUE),
         param2 = try(beta.mle(.$value)$param[2],silent=TRUE)) %>%
      ungroup() %>%
      mutate(param1 = as.vector(unlist(.$param1)),
             param2 = as.vector(unlist(.$param2))) %>%
      mutate(Sensitive = ifelse(Sensitive=='Sensitive',1,0))
    for (stud in unique(sens.short$experiment)) {
      for(meas in unique(sens.short$measure)) {
        sens.temp <- sens.short %>%
          filter(experiment == stud,
                 measure == meas) %>%
          select(value) %>%
          mutate(value = ifelse(value==0,1e-3,value)) %>%
          as.data.frame()
        param1 <- params %>%
          filter(experiment == stud,
                 measure == meas) %>%
          pull(param1)
        param2 <- params %>%
          filter(experiment == stud,
                 measure == meas) %>%
          pull(param2)
        fit <- try(betamix(value~1|1, data=sens.temp, k=2),
                   silent=TRUE)
        if (class(fit) != 'try-error') {
          if (!is.null(nrow(coef(fit)))) {
            mu <- plogis(coef(fit)[,1])
            phi <- exp(coef(fit)[,2])
            alpha <- mu * phi
            beta <- (1 - mu) * phi
            sens.new <- sens.short %>%
              filter(experiment == stud,
                     measure == meas) %>%
              mutate(posterior_Z0 = fit$flexmix@posterior$scaled[,1],
                     posterior_Z1 = fit$flexmix@posterior$scaled[,2],
                     posterior_Zmax = clusters(fit)-1) %>%
              bind_rows(sens.new)
            params.new <- params %>%
              mutate(param1 = as.character(param1),
                     param2 = as.character(param2)) %>%
              filter(experiment == stud,
                     measure == meas) %>%
              arrange(Sensitive) %>%
              mutate(posterior_alpha = alpha,
                     posterior_beta = beta,
                     posterior_lambda = c(1-mean(clusters(fit)-1),mean(clusters(fit)-1)),
                     loglik = fit$flexmix@logLik,
                     drug = dr) %>%
              bind_rows(params.new)
          } else {
            mu <- plogis(coef(fit)[1])
            phi <- exp(coef(fit)[2])
            alpha <- mu * phi
            beta <- (1 - mu) * phi
            sens.new <- sens.short %>%
              filter(experiment == stud,
                     measure == meas) %>%
              mutate(posterior_Z0 = fit$flexmix@posterior$scaled[,1],
                     posterior_Z1 = 0,
                     posterior_Zmax = clusters(fit)-1) %>%
              bind_rows(sens.new)
            params.new <- params %>%
              mutate(param1 = as.character(param1),
                     param2 = as.character(param2)) %>%
              filter(experiment == stud,
                     measure == meas) %>%
              arrange(Sensitive) %>%
              mutate(posterior_alpha = alpha,
                     posterior_beta = beta,
                     posterior_lambda = 0,
                     loglik = fit$flexmix@logLik,
                     drug = dr) %>%
              bind_rows(params.new)
          }
        } else {
          mu <- NA
          phi <- NA
          alpha <- NA
          beta <- NA
          sens.new <- sens.short %>%
            filter(experiment == stud,
                   measure == meas) %>%
            mutate(posterior_Z0 = NA,
                   posterior_Z1 = NA,
                   posterior_Zmax = NA) %>%
            bind_rows(sens.new)
          params.new <- params %>%
            filter(experiment == stud,
                   measure == meas) %>%
            arrange(Sensitive) %>%
            mutate(posterior_alpha = alpha,
                   posterior_beta = beta,
                   posterior_lambda = c(NA,NA),
                   loglik = NA,
                   drug = dr) %>%
            bind_rows(params.new)
        }
      }
    }
  }
  return(list(sens.new,params.new))
}
