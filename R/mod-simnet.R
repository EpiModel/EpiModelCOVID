#' @rdname moduleset-vaxDecisions
#' @export
resim_nets_covid_vax_decisions <- function(dat, at) {
 
  dat$num.nw <- dat$num.nw - 1 # disregard household nw for this module
  
  ## Edges correction
  dat <- edges_correct(dat, at)

  ## Network Re-simulation
  dat.updates <- NVL(get_control(dat, "dat.updates"), function(dat, ...) dat)
  dat <- dat.updates(dat = dat, at = at, network = 0L)
  for (network in seq_len(dat$num.nw)) {
    dat <- simulate_dat(dat = dat, at = at, network = network)
    dat <- dat.updates(dat = dat, at = at, network = network)
  }
  
  dat <- summary_nets(dat, at)
  
  dat$num.nw <- dat$num.nw + 1

  return(dat)
}


