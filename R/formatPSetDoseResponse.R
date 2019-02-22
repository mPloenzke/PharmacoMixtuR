#' Format dose-response measurements
#'
#' Returns formatted raw viability and concetration measurements.
#'
#' @param pset Individual \code{PSet}
#'
#' @return Tibble containing raw dose-response measurements
#'
#' @export
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @seealso \code{\link{intersectPSetWrapper}}
#' @keywords intersectPSet Sensitivity
#'
#' @importFrom dplyr mutate select distinct left_join rename ends_with mutate_if
#' @importFrom tibble as.tibble rownames_to_column
#' @importFrom tidyr gather separate
#' @importFrom magrittr %>%
formatPSetDoseResponse <- function(pset) {
  if (pset@annotation$name %in% c('CCLE','GDSC')) {
    drugs <- rownames_to_column(as.data.frame(pset@sensitivity$info)) %>%
      as.tibble() %>%
      mutate(drug.name = if(exists('drugid',where=.)) drugid else drug.name) %>%
      select(drugid, drug.name) %>%
      distinct() %>%
      mutate(drug.name = tolower(drug.name))
    sens.info <- rownames_to_column(as.data.frame(pset@sensitivity$raw)) %>%
      as.tibble() %>%
      separate(rowname,into=c('cellid','drugid'),sep='_')
    concentrations <- sens.info %>%
      select(ends_with('Dose'),cellid,drugid) %>%
      mutate_if(is.factor, as.character) %>%
      gather(key='temp',value='Concentration',-cellid,-drugid) %>%
      mutate(study = pset@annotation$name) %>%
      separate(temp,into=c('dose','temp'),sep='.D') %>%
      select(-temp)
    sens.curves <-  sens.info %>%
      select(ends_with('Viability'),cellid,drugid) %>%
      mutate_if(is.factor, as.character) %>%
      gather(key='temp',value='Viability',-cellid,-drugid) %>%
      mutate(study = pset@annotation$name) %>%
      separate(temp,into=c('dose','temp'),sep='.V') %>%
      select(-temp) %>%
      left_join(concentrations,by=c('dose','cellid','drugid','study')) %>%
      left_join(drugs, by=c('drugid')) %>%
      select(-drugid, -dose) %>%
      rename(drug=drug.name)
    return(sens.curves)
  } else if (pset@annotation$name == 'gCSI') {
    drugs <- as.data.frame(pset@sensitivity$info) %>%
      as.tibble() %>%
      mutate(drug.name = if(exists('drugid',where=.)) drugid else drug.name) %>%
      select(drugid, drug.name) %>%
      distinct() %>%
      mutate(drug.name = tolower(drug.name)) %>%
      mutate_if(is.factor, as.character)
    sens.info <- rownames_to_column(as.data.frame(stud@sensitivity$raw)) %>%
      as.tibble() %>%
      separate(rowname,into=c('cellid','drugid'),sep='_')
    concentrations <- sens.info %>%
      select(ends_with('Dose'),cellid,drugid) %>%
      mutate_if(is.factor, as.character) %>%
      gather(key='temp',value='Concentration',-cellid,-drugid) %>%
      mutate(study = pset@annotation$name) %>%
      separate(temp,into=c('dose','temp'),sep='.D') %>%
      select(-temp)
    sens.curves <-  sens.info %>%
      select(ends_with('Viability'),cellid,drugid) %>%
      mutate_if(is.factor, as.character) %>%
      gather(key='temp',value='Viability',-cellid,-drugid) %>%
      mutate(study = pset@annotation$name) %>%
      separate(temp,into=c('dose','temp'),sep='.V') %>%
      select(-temp) %>%
      left_join(concentrations,by=c('dose','cellid','drugid','study')) %>%
      left_join(drugs, by=c('drugid')) %>%
      select(-drugid, -dose) %>%
      rename(drug=drug.name) %>%
      mutate_if(is.numeric, as.character)
    return(sens.curves)
  } else if (pset@annotation$name == 'FIMM') {
    sens.info <- rownames_to_column(as.data.frame(pset@sensitivity$raw)) %>%
      as.tibble() %>%
      separate(rowname,into=c('drugid','cellid'),sep='_')
    concentrations <- sens.info %>%
      select(ends_with('Dose'),cellid,drugid) %>%
      mutate_if(is.factor, as.character) %>%
      gather(key='temp',value='Concentration',-cellid,-drugid) %>%
      mutate(study = pset@annotation$name) %>%
      separate(temp,into=c('dose','temp'),sep='.D') %>%
      select(-temp)
    sens.curves <-  sens.info %>%
      select(ends_with('Viability'),cellid,drugid) %>%
      mutate_if(is.factor, as.character) %>%
      gather(key='temp',value='Viability',-cellid,-drugid) %>%
      mutate(study = pset@annotation$name) %>%
      separate(temp,into=c('dose','temp'),sep='.V') %>%
      select(-temp) %>%
      left_join(concentrations,by=c('dose','cellid','drugid','study')) %>%
      select(-dose) %>%
      rename(drug=drugid) %>%
      mutate(drug=tolower(drug)) %>%
      mutate_if(is.numeric, as.character)
    return(sens.curves)
  } else if (pset@annotation$name == 'GDSC1000') {
    drugs <- as.data.frame(pset@sensitivity$info) %>%
      as.tibble() %>%
      mutate(drug.name = if(exists('drugid',where=.)) drugid else drug.name) %>%
      select(drug.name,DRUG_ID) %>%
      distinct() %>%
      mutate(drug.name = tolower(drug.name),
             drugid = DRUG_ID) %>%
      select(-DRUG_ID) %>%
      mutate_if(is.factor, as.character) %>%
      mutate_if(is.numeric, as.character)
    sens.info <- rownames_to_column(as.data.frame(pset@sensitivity$raw)) %>%
      as.tibble() %>%
      separate(rowname,into=c('cellid','drugid'),sep='_')
    concentrations <- sens.info %>%
      select(ends_with('Dose'),cellid,drugid) %>%
      mutate_if(is.factor, as.character) %>%
      gather(key='temp',value='Concentration',-cellid,-drugid) %>%
      mutate(study = pset@annotation$name) %>%
      separate(temp,into=c('dose','temp'),sep='.D') %>%
      select(-temp)
    sens.curves <-  sens.info %>%
      select(ends_with('Viability'),cellid,drugid) %>%
      mutate_if(is.factor, as.character) %>%
      gather(key='temp',value='Viability',-cellid,-drugid) %>%
      mutate(study = pset@annotation$name) %>%
      separate(temp,into=c('dose','temp'),sep='.V') %>%
      select(-temp) %>%
      left_join(concentrations,by=c('dose','cellid','drugid','study')) %>%
      left_join(drugs, by=c('drugid')) %>%
      select(-dose) %>%
      rename(drug=drugid) %>%
      mutate(drug=tolower(drug)) %>%
      mutate_if(is.numeric, as.character)
    return(sens.curves)
  } else {
    print('Unknown PSet')
  }
}
