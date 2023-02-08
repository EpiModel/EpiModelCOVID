
#' @rdname moduleset-common
#' @export

# USING J&J VACCINE AS A REFERENCE

vax_covid <- function(dat, at) {
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  age <- get_attr(dat, "age")
  vax <- get_attr(dat, "vax")
  vaxTime <- get_attr(dat, "vaxTime")

  vax.start <- get_param(dat, "vax.start")
  vax.rate <- get_param(dat, "vax.rate")
  vax.partial.immune <- get_param(dat, "vax.partial.immune")
  vax.full.immune <- get_param(dat, "vax.full.immune")

  ## First vax
  nVax <- 0
  if (at >= vax.start) {
    idsElig.vax <- which(active == 1 & status == "s" & vax == 0)
    nElig.vax <- length(idsElig.vax)
    if (nElig.vax > 0) {
      age.breaks <- c(17,20,30,40,50,60,85)
      age.group <- as.numeric(cut(age[idsElig.vax], age.breaks, labels = FALSE, right = FALSE))
      vax.rate.vec <- vax.rate[age.group]
      vecVax <- which(rbinom(nElig.vax, 1, vax.rate.vec) == 1)
      idsVax <- idsElig.vax[vecVax]
      nVax <- length(idsVax)
      if (nVax > 0) {
        vax[idsVax] <- 1
        vaxTime[idsVax] <- at
      }
    }
  }

  idsVax.gt65 <- which(active == 1 & vax == 1 & vaxTime == at & age >= 65)
  nidsVax.gt65 <- length(idsVax.gt65)

  idsVax.17to65 <- which(active == 1 & vax == 1 & vaxTime == at & age < 65 &
                        age >= 17)
  nidsVax.17to65 <- length(idsVax.17to65)

  # Partial Immunity after Single Dose
  idsvaximmunePartial <- which(active == 1 & vax == 1 & at - vaxTime >= vax.partial.immune)
  nvaximmunePartial <- length(idsvaximmunePartial)
  if (nvaximmunePartial > 0) {
    vax[idsvaximmunePartial] <- 2
  }

  # Full Immunity after Single Dose
  idsvaximmuneFull <- which(active == 1 & vax == 2 & at - vaxTime >= vax.full.immune)
  nvaximmuneFull <- length(idsvaximmuneFull)
  if (nvaximmuneFull > 0) {
    vax[idsvaximmuneFull] <- 3
  }
    
  ## Replace attr
  dat <- set_attr(dat, "vax", vax)
  dat <- set_attr(dat, "vaxTime", vaxTime)

  ## Summary statistics ##
  dat <- set_epi(dat, "nVax", at, nVax)
  dat <- set_epi(dat, "nVaxImmunePart", at, nvaximmunePartial)
  dat <- set_epi(dat, "nVaxImmuneFull", at, nvaximmuneFull)
  dat <- set_epi(dat, "nVaxgt65", at, nidsVax.gt65)
  dat <- set_epi(dat, "nVax17to65", at, nidsVax.17to65)

  return(dat)
}
