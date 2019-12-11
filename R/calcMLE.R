#' Wrapper to calculate sensitivity profile
#'
#' A wrapper to return the sensitivity information for the
#' intersected PSets across cells and drugs.
#'
#' @param sensitivitiy.tibb Tibble containing sensitivity information, drugs, and cells
#' @param drug.names List of drugs to fit model to 
#' @param fit Return MLE for beta or normal distribution
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
#' Wrapper to calculate sensitivity profile
#'
#' A wrapper to return the sensitivity information for the
#' intersected PSets across cells and drugs.
#'
#' @param x Number of sucessess in beta-binomial
#' @param n Number of trials
#'
#' @return List of parameter estimates
#'
#' @export
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @seealso \code{\link{iterate_em}} \code{\link{fit_beta_mle}}
#' @keywords beta binomial mle
#'
#' @importFrom VGAM dbetabinom.ab
#' @importFrom stats4 mle coef
fit_bb_mle <- function(x, n) {
  ll <- function(alpha, beta) {
    -sum(dbetabinom.ab(x, n, alpha, beta, log = TRUE))
  }
  m <- try({
    stats4::mle(ll, start = list(alpha = 3, beta = 10), method = "L-BFGS-B",
                   lower = c(0.001, .001))
  }, silent=TRUE)
  if (class(m) == "try-error") {
    ab <- c(1e-2,1e-2)
  } else {
    ab <- stats4::coef(m)
  }
  data_frame(alpha = ab[1], beta = ab[2], number = length(x))
}
#' Wrapper to calculate sensitivity profile
#'
#' A wrapper to return the sensitivity information for the
#' intersected PSets across cells and drugs.
#'
#' @param x Observed value in beta
#'
#' @return List of parameter estimates
#'
#' @export
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @seealso \code{\link{iterate_em}} \code{\link{fit_bb_mle}}
#' @keywords beta mle
#'
#' @importFrom stats dbeta
#' @importFrom stats4 mle coef
fit_beta_mle <- function(x) {
  
  ll <- function(alpha, beta) {
    -sum(dbeta(x2, alpha, beta, log = TRUE))
  }
  if (length(x) == 1) {
    x2 <- c(x, x+rnorm(1,0,sd=.05))
  } else {
    x2 <- x
  }
  m <- try(stats4::mle(ll, start = list(alpha = 3, beta = 10), method = "L-BFGS-B",
                   lower = c(0.001, .001)),silent=TRUE)
  if (class(m) == "try-error") {
    ab <- c(1e-2,1e-2)
  } else {
    ab <- stats4::coef(m)
  }
  data_frame(alpha = ab[1], beta = ab[2], number = length(x))
}
