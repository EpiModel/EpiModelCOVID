
#' @rdname moduleset-common
#' @export
vax_covid <- function(dat, at) {
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  vax <- get_attr(dat, "vax")
  vax1Time <- get_attr(dat, "vax1Time")
  vax2Time <- get_attr(dat, "vax2Time")

  vax.start <- get_param(dat, "vax.start")
  vax1.rate <- get_param(dat, "vax1.rate")
  vax2.interval <- get_param(dat, "vax2.interval")
  vax1.immune <- get_param(dat, "vax1.immune")
  vax2.immune <- get_param(dat, "vax2.immune")

  ## First vax
  nVax <- 0
  if (at >= vax.start) {
    idsElig.vax1 <- which(active == 1 & status == "s" & vax == 0)
    nElig.vax1 <- length(idsElig.vax1)
    if (nElig.vax1 > 0) {
      vecVax <- which(rbinom(nElig.vax1, 1, vax1.rate) == 1)
      idsVax <- idsElig.vax1[vecVax]
      nVax <- length(idsVax)
      if (nVax > 0) {
        vax[idsVax] <- 1
        vax1Time[idsVax] <- at
      }
    }
  }

  # Partial Immunity after first shot
  idsvaximmunePartial <- which(active == 1 & vax == 1 & at - vax1Time >= vax1.immune)
  nvaximmunePartial <- length(idsvaximmunePartial)
  if (nvaximmunePartial > 0){
    vax[idsvaximmunePartial] <- 2
  }

  ## Second vax
  idsvaxFull <- which(active == 1 & vax == 2 & (at - vax1Time >= vax2.interval))
  nvaxFull <- length(idsvaxFull)
  if (nvaxFull > 0) {
    vax[idsvaxFull] <- 3
    vax2Time[idsvaxFull] <- at
  }

  # Partial Immunity after first shot
  idsvaximmuneFull <- which(active == 1 & vax == 3 & at - vax2Time >= vax2.immune)
  nvaximmuneFull <- length(idsvaximmuneFull)
  if (nvaximmuneFull > 0){
    vax[idsvaximmuneFull] <- 4
  }

  ## Replace attr
  dat <- set_attr(dat, "vax", vax)
  dat <- set_attr(dat, "vax1Time", vax1Time)
  dat <- set_attr(dat, "vax2Time", vax2Time)

  ## Summary statistics ##
  dat <- set_epi(dat, "nVax1", at, nVax)
  dat <- set_epi(dat, "nVax2", at, nvaxFull)
  dat <- set_epi(dat, "nVaxImmunePart", at, nvaximmunePartial)

  return(dat)
}
