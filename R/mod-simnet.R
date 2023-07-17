#' @rdname moduleset-vaxDecisions
#' @export
resim_nets_covid_vax_decisions <- function(dat, at) {
 
  ## Edges correction
  dat <- edges_correct(dat, at)

  ## network resimulation
  dat.updates <- NVL(get_control(dat, "dat.updates"), function(dat, ...) dat)
  dat <- dat.updates(dat = dat, at = at, network = nets[1L] - 1L)
  for (network in nets) {
    dat <- simulate_dat(dat = dat, at = at, network = network)
    dat <- dat.updates(dat = dat, at = at, network = network)
  }

  return(dat)
}


