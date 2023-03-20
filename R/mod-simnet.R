
#' @rdname moduleset-ship
#' @export
resim_nets_covid_ship <- function(dat, at) {
  ## Edges correction
  dat <- edges_correct(dat, at)

  if (at < dat$param$network.lockdown.time) {
    nets <- 1:3
  } else {
    nets <- 4:6
  }

  ## network resimulation
  dat.updates <- NVL(get_control(dat, "dat.updates"), function(dat, ...) dat)
  dat <- dat.updates(dat = dat, at = at, network = nets[1L] - 1L)
  for (network in nets) {
    dat <- simulate_dat(dat = dat, at = at, network = network)
    dat <- dat.updates(dat = dat, at = at, network = network)
  }

  return(dat)
}


#' @rdname moduleset-corporate
#' @export
resim_nets_covid_corporate <- resim_nets
