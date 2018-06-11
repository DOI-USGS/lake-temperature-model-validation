# wants and needs helper functions

wqp_calc_wants <- function(wqp_sites, wqp_variables) {
  
  vars <- names(wqp_variables$characteristicName)
  df <- expand.grid(wqp_sites, vars, stringsAsFactors = FALSE)
  names(df) <- c('site', 'variable')
  return(df)
}

# wqp_wants - wqp_haves = wqp_needs
wqp_calc_needs <- function(wqp_wants, wqp_haves) {
  
  haves <- readRDS(wqp_haves)
  
  # leaves only the rows in wqp_wants that don't exist in wqp_haves
  # these are site/variable combinations where we have no data
  diffed_cells <- dplyr::anti_join(wqp_wants, haves,  by = c('site', 'variable')) 
  
  if (any(is.na(diffed_cells))){
    # shouldn't be any NAs, but if there are, throw an error
    stop('found NA(s) in NLDAS cell diff. Check ', wqp_haves, ' and `wq_wants` data')
  }
  
  # return the site/variable combinations what we need to pull
  return(diffed_cells)
}
