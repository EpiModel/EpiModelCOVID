
#' @rdname moduleset-common
#' @export
vax_covid <- function(dat, at) {
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  age <- get_attr(dat, "age")
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
      age.group <- pmin((floor(age[idsElig.vax1] / 10)) + 1, 8)
      vax1.rate.vec <- vax1.rate[age.group]
      vecVax <- which(rbinom(nElig.vax1, 1, vax1.rate.vec) == 1)
      idsVax <- idsElig.vax1[vecVax]
      nVax <- length(idsVax)
      if (nVax > 0) {
        vax[idsVax] <- 1
        vax1Time[idsVax] <- at
      }
    }
  }

  idsVax1.gt65 <- which(active == 1 & vax == 1 & vax1Time == at & age >= 65)
  nidsVax1.gt65 <- length(idsVax1.gt65)

  idsVax1.15to65 <- which(active == 1 & vax == 1 & vax1Time == at & age < 65 &
                        age >= 15)
  nidsVax1.15to65 <- length(idsVax1.15to65)

  idsVax1.lt15 <- which(active == 1 & vax == 1 & vax1Time == at & age < 15)
  nidsVax1.lt15 <- length(idsVax1.lt15)

  # Partial Immunity after first shot
  idsvaximmunePartial <- which(active == 1 & vax == 1 & at - vax1Time >= vax1.immune)
  nvaximmunePartial <- length(idsvaximmunePartial)
  if (nvaximmunePartial > 0) {
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
  if (nvaximmuneFull > 0) {
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
  dat <- set_epi(dat, "nVax1gt65", at, nidsVax1.gt65)
  dat <- set_epi(dat, "nVax115to65", at, nidsVax1.15to65)
  dat <- set_epi(dat, "nVax1lt15", at, nidsVax1.lt15)

  return(dat)
}
