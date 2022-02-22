
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
  dat <- set_epi(dat, "nVax1gt65", at, nidsVax1.gt65)
  dat <- set_epi(dat, "nVax115to65", at, nidsVax1.15to65)
  dat <- set_epi(dat, "nVax1lt15", at, nidsVax1.lt15)

  return(dat)
}

#' @rdname moduleset-boost
#' @export
#'
vax_covid_boost <- function(dat, at) {
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  age <- get_attr(dat, "age")
  vax <- get_attr(dat, "vax")
  vax1Time <- get_attr(dat, "vax1Time")
  vax2Time <- get_attr(dat, "vax2Time")
  vax3Time <- get_attr(dat, "vax3Time")

  vax.start <- get_param(dat, "vax.start")
  vax1.rate <- get_param(dat, "vax1.rate")
  vax2.rate <- get_param(dat, "vax2.rate")
  vax3.rate <- get_param(dat, "vax3.rate")
  vax2.interval <- get_param(dat, "vax2.interval")
  vax3.interval <- get_param(dat, "vax3.interval")
  vax1.immune <- get_param(dat, "vax1.immune")
  vax2.immune <- get_param(dat, "vax2.immune")
  vax3.immune <- get_param(dat, "vax3.immune")

  ## Vax1
  nVax <- 0
  if (at >= vax.start) {
    idsElig.vax1 <- which(active == 1 & status == "s" & vax == 0)
    nElig.vax1 <- length(idsElig.vax1)
    if (nElig.vax1 > 0) {
      age.group <- ifelse(floor(age[idsElig.vax1]) <= 5, 1,
                          ifelse(floor(age[idsElig.vax1]) <= 11, 2,
                                 ifelse(floor(age[idsElig.vax1]) <= 17, 3,
                                        ifelse(floor(age[idsElig.vax1]) <= 64, 4,
                                               5))))
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

  # Partial Immunity after Vax1
  idsimmune.vax1 <- which(active == 1 & vax == 1 & at - vax1Time >= vax1.immune)
  nimmune.vax1 <- length(idsimmune.vax1)
  if (nimmune.vax1 > 0){
    vax[idsimmune.vax1] <- 2
  }

  ## Vax2
  nVax2 <- 0
  idsElig.vax2 <- which(active == 1 & vax == 2 & (at - vax1Time >= vax2.interval))
  nElig.vax2 <- length(idsElig.vax2)
  if (nElig.vax2 > 0){
    age.group <- ifelse(floor(age[idsElig.vax2]) <= 5, 1,
                        ifelse(floor(age[idsElig.vax2]) <= 11, 2,
                               ifelse(floor(age[idsElig.vax2]) <= 17, 3,
                                      ifelse(floor(age[idsElig.vax2]) <= 64, 4,
                                             5))))
    vax2.rate.vec <- vax2.rate[age.group]
    vecVax2 <- which(rbinom(nElig.vax2, 1, vax2.rate.vec) == 1)
    idsVax2 <- idsElig.vax2[vecVax2]
    nVax2 <- length(idsVax2)
    if (nVax2 > 0){
      vax[idsVax2] <- 3
      vax2Time[idsVax2] <- at
    }
  }
  # Immunity after Vax2
  idsimmune.vax2 <- which(active == 1 & vax == 3 & at - vax2Time >= vax2.immune)
  nimmune.vax2 <- length(idsimmune.vax2)
  if (nimmune.vax2 > 0){
    vax[idsimmune.vax2] <- 4
  }

  ## Vax3 - Boost
  nVax3 <- 0
  nVax3.65 <- 0
  nVax3.18 <- 0
  nVax3.50 <- 0
  idsElig.vax3 <- which(active == 1 & vax == 4 & (at - vax2Time >= vax3.interval))
  nElig.vax3 <- length(idsElig.vax3)
  if (nElig.vax3 > 0){
    age.group <- ifelse(floor(age[idsElig.vax3]) <= 18, 1,
                        ifelse(floor(age[idsElig.vax2]) <= 50, 2,
                                             3))
    vax3.rate.vec <- vax3.rate[age.group]
    vecVax3 <- which(rbinom(nElig.vax3, 1, vax3.rate) == 1)
    idsVax3 <- idsElig.vax3[vecVax3]
    nVax3 <- length(idsVax3)
    if (nVax3 > 0){
      vax[idsVax3] <- 5
      vax3Time[idsVax3] <- at
    }
  }


  # Immunity after Vax3
  idsimmune.vax3 <- which(active == 1 & vax == 5 & at - vax3Time >= vax3.immune)
  nimmune.vax3 <- length(idsimmune.vax3)
  if (nimmune.vax3 > 0){
    vax[idsimmune.vax3] <- 6
  }

  ## Replace attr
  dat <- set_attr(dat, "vax", vax)
  dat <- set_attr(dat, "vax1Time", vax1Time)
  dat <- set_attr(dat, "vax2Time", vax2Time)
  dat <- set_attr(dat, "vax3Time", vax3Time)

  ## Summary statistics ##
  dat <- set_epi(dat, "nVax1", at, nVax)
  dat <- set_epi(dat, "nVax2", at, nVax2)
  dat <- set_epi(dat, "nVax3", at, nVax3)
  dat <- set_epi(dat, "nVaxImmunePart", at, nimmune.vax2)

  return(dat)
}
