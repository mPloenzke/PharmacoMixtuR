#' Wrapper to calculate beta mixture model
#'
#' A wrapper to return the sensitivity information for the
#' intersected PSets across cells and drugs.
#'
#' @param state
#'
#' @return List of tibbles; one containing posterior estimates and one containing estimated parameters
#'
#' @export
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @seealso \code{\link{intersectPSetWrapper}}
#' @keywords fit_bb_mle fit_beta_mle
#'
#' @importFrom magrittr %>%
#' @importFrom VGAM dbetabinom.ab
#' @importFrom dplyr group_by do mutate ungroup select top_n n_distinct left_join filter bind_rows
#' @importFrom tidyr crossing
#' @importFrom stats qbeta
iterate_em <- function(state, ...) {

  # M-step for drug-specific parameters (targetted)
  targetted_drug_fits_temp <- state$cell_assignments %>% 
    mutate(value = case_when(value==0 ~  1e-3, TRUE ~ value)) %>%
    mutate(drug_type = 'targeted') %>%
    group_by(drug, cell_type, drug_type, experiment) %>% 
    do(mutate(fit_beta_mle(.$value),number = nrow(.))) %>% 
    group_by(drug, experiment) %>%
    mutate(ndist=n_distinct(cell_type)) %>%
    ungroup()
  mdat_temp <- state$cell_assignments %>% 
    select(drug, experiment, cell_type, targeted_cell_type_prior) %>% 
    distinct() %>%
    group_by(drug, experiment) %>%
    mutate(ndist=n_distinct(cell_type)) %>%
    ungroup()
  mdat_temp <- mdat_temp %>%
    filter(ndist<2) %>%
    mutate(cell_type2 = case_when(cell_type=='resistant' ~ 'sensitive',cell_type=='sensitive' ~ 'resistant'),
           cell_type=cell_type2,
           targeted_cell_type_prior = 1-targeted_cell_type_prior) %>%
    select(-cell_type2) %>%
    bind_rows(mdat_temp) %>%
    select(-ndist)
  targetted_drug_fits <- targetted_drug_fits_temp %>%
    bind_rows(
      targetted_drug_fits_temp %>% 
        filter(ndist<2) %>% 
        mutate(cell_type2 = case_when(cell_type=='resistant' ~ 'sensitive',cell_type=='sensitive' ~ 'resistant'),
               cell_type=cell_type2,
               number=0,
               alpha=case_when(cell_type == 'resistant' ~ .5,
                               cell_type == 'sensitive' ~ 50),
               beta=case_when(cell_type == 'resistant' ~ 7.5,
                              cell_type == 'sensitive' ~ 50)) %>%
        select(-cell_type2)
    ) %>%
    select(-ndist, -number) %>%
    left_join(mdat_temp, by=c('drug','experiment','cell_type'))
  
  # M-step for drug-specific parameters (broad)
  broad_drug_fits <- state$cell_assignments %>%
    mutate(x=ifelse(value==0,1e-3,value)) %>%
    mutate(drug_type = 'broad') %>%
    group_by(drug, drug_type, experiment) %>%
    do(mutate(fit_beta_mle(.$x),number = nrow(.))) %>%
    ungroup()  %>% 
    select(-number) %>%
    left_join(state$cell_assignments %>% select(drug, experiment, broad_cell_type_prior) %>% distinct(), by=c('drug','experiment'))
  
  # E-step for cell types (targetted drugs)
  targetted_drug_list <- unique(targetted_drug_fits$drug)
  targetted_cell_assignments <- lapply(targetted_drug_list, function(dr) {
    state$cell_assignments %>% 
      filter(drug == dr) %>% 
      select(drug:measure) %>%
      mutate(drug_type='targeted') %>%
      mutate(value=ifelse(value==0,1e-3,value)) %>%
      crossing(targetted_drug_fits %>% filter(drug == dr) %>% select(-drug, -drug_type) %>% rename(experiment1=experiment)) %>% 
      filter(experiment==experiment1) %>% 
      select(-experiment1) %>% 
      mutate(median = qbeta(.5,alpha,beta),
             likelihood = dbeta(value, alpha, beta)) %>% 
      group_by(cell_type) %>%
      mutate(max.likelihood=max(likelihood)) %>%
      ungroup() %>%
      mutate(likelihood = case_when(value>median & cell_type == 'sensitive' ~ max.likelihood-likelihood + max.likelihood,
                                    value<median & cell_type == 'resistant' ~ max.likelihood-likelihood + max.likelihood, 
                                    TRUE~likelihood)) %>%
      mutate(likelihood = targeted_cell_type_prior * likelihood) %>%
      group_by(drug, cell, cell_type) %>%
      mutate(likelihood = sum(likelihood)) %>%
      ungroup() %>%
      group_by(drug, cell, experiment) %>%
      mutate(posterior=likelihood/sum(likelihood)) %>%
      top_n(1, likelihood) %>% 
      ungroup()
  })
  targetted_cell_assignments <- do.call(rbind, targetted_cell_assignments)
  
  # E-step for cell types (broad drugs)
  broad_drug_list <- unique(broad_drug_fits$drug)
  broad_cell_assignments <- state$cell_assignments %>%
    filter(drug %in% broad_drug_list) %>%
    select(drug:measure) %>%
    mutate(drug_type='broad') %>%  
    mutate(value=ifelse(value==0,1e-3,value)) %>%
    left_join(broad_drug_fits, by=c('drug','drug_type','experiment')) %>%
    mutate(likelihood = 1 * dbeta(value, alpha, beta),
           cell_type = ifelse((broad_cell_type_prior*pbeta(value, alpha, beta)) >
                                            ((1-broad_cell_type_prior)*(1-pbeta(value, alpha, beta))),
                                          'sensitive','resistant'),
           posterior = (broad_cell_type_prior*pbeta(value, alpha, beta)) +
             ((1-broad_cell_type_prior)*(1-pbeta(value, alpha, beta))), 
           posterior = 1-posterior)
  
  # Format all cell assignments
  cell_assignments <- bind_rows(targetted_cell_assignments,broad_cell_assignments) %>%
    group_by(drug, cell_type, experiment) %>%
    fill(targeted_cell_type_prior, .direction='down') %>%
    ungroup() %>%
    group_by(drug, experiment) %>% 
    fill(broad_cell_type_prior, .direction='up') %>%
    ungroup() %>% 
    right_join(state$drug_assignments %>% select(drug, drug_type), by=c('drug','drug_type')) %>%
    select(drug, cell, value, experiment, measure, drug_type, cell_type, broad_cell_type_prior, targeted_cell_type_prior, posterior, likelihood)
  
  # M-step for drug-type parameters
  drug_types <- cell_assignments %>%
    group_by(drug, measure, drug_type) %>%
    summarize(sensitive = sum(cell_type=='sensitive'),
              cell_count = n()) %>%
    group_by(drug_type) %>%
    ungroup() %>% 
    group_by(drug_type) %>%
    do(mutate(fit_bb_mle(.$sensitive, .$cell_count), number = nrow(.))) %>%
    ungroup() %>%
    left_join(state$drug_assignments %>% distinct(drug_type, drug_type_prior), by='drug_type')
  
  # E-step for drug types
  targetteds <- state$drug_assignments %>%
    select(drug, measure) %>%
    right_join(targetted_cell_assignments %>% select(drug, cell, cell_type), by='drug') %>% 
    group_by(drug, measure) %>%
    summarise(sensitive = sum(cell_type=='sensitive'),
              cell_count = n()) %>%
    ungroup() %>% 
    mutate(type='targeted')
  broads <- state$drug_assignments %>%
    select(drug, measure) %>%
    right_join(broad_cell_assignments %>% select(drug, cell, cell_type), by='drug') %>% 
    group_by(drug, measure) %>%
    summarise(sensitive = sum(cell_type=='sensitive'),
              cell_count = n()) %>%
    ungroup() %>%
    mutate(type='broad')
  drugs.tibble <- bind_rows(targetteds, broads)
  drug_assignments <- drugs.tibble %>% 
    rename(drug_type = type) %>%
    right_join(drug_types %>% select(-number), by='drug_type') %>%
    mutate(likelihood = drug_type_prior*VGAM::dbetabinom.ab(sensitive, cell_count, alpha, beta)) %>%
    group_by(drug) %>% 
    mutate(posterior=likelihood/sum(likelihood)) %>%
    top_n(1, likelihood) %>%
    ungroup()
  #if (drug_assignments %>% filter(drug=='Crizotinib') %>% distinct(drug_type) %>% pull() == 'broad') {browser()}
  
  drug_fits <- targetted_drug_fits %>%
    select(drug, drug_type, cell_type, experiment, alpha, beta) %>%
    bind_rows(broad_drug_fits %>% select(drug, drug_type, experiment, alpha, beta) %>% mutate(cell_type = 'resistant')) %>%
    right_join(drug_assignments %>% select(drug, drug_type), by=c('drug','drug_type'))
  
  list(drug_assignments = drug_assignments,
       drug_types = drug_types,
       cell_assignments = cell_assignments,
       drug_fits = drug_fits)
}
