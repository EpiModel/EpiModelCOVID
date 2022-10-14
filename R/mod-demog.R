
#' @rdname moduleset-common
#' @export
aging_covid <- function(dat, at) {

  age <- get_attr(dat, "age")
  age <- age + 1 / 365
  dat <- set_attr(dat, "age", age)

  return(dat)
}


#' @rdname moduleset-ship
#' @export
deaths_covid_ship <- function(dat, at) {

  ## Attributes ##
  active <- dat$attr$active
  age <- dat$attr$age
  race <- dat$attr$race
  status <- dat$attr$status

  ## Parameters ##
  mort.rates <- dat$param$mort.rates
  mort.dis.mult <- dat$param$mort.dis.mult

  idsElig <- which(active == 1)
  nElig <- length(idsElig)
  nDeaths <- nDeathsIC <- 0

  if (nElig > 0) {

    whole_ages_of_elig <- pmin(ceiling(age[idsElig]), 86)
    death_rates_of_elig <- mort.rates[whole_ages_of_elig, race]

    idsElig.inf <- which(status[idsElig] == "ic")
    death_rates_of_elig[idsElig.inf] <- death_rates_of_elig[idsElig.inf] * mort.dis.mult

    vecDeaths <- which(rbinom(nElig, 1, death_rates_of_elig) == 1)
    idsDeaths <- idsElig[vecDeaths]
    nDeaths <- length(idsDeaths)
    nDeathsIC <- length(intersect(idsDeaths, idsElig.inf))

    if (nDeaths > 0) {
      dat$attr$active[idsDeaths] <- 0
      inactive <- which(dat$attr$active == 0)
      dat$attr <- deleteAttr(dat$attr, inactive)
      for (i in seq_along(dat$el)) {
        dat$el[[i]] <- delete_vertices(dat$el[[i]], inactive)
      }
    }
  }

  ## Summary statistics ##
  dat$epi$d.flow[at] <- nDeaths
  dat$epi$d.ic.flow[at] <- nDeathsIC

  return(dat)
}


#' @rdname moduleset-corporate
#' @export
deaths_covid_corporate <- function(dat, at) {

  ## Attributes ##
  active <- get_attr(dat, "active")
  age <- get_attr(dat, "age")
  race <- get_attr(dat, "race")
  status <- get_attr(dat, "status")

  ## Parameters ##
  mort.rates <- get_param(dat, "mort.rates")
  mort.dis.mult <- get_param(dat, "mort.dis.mult")

  idsElig <- which(active == 1)
  nElig <- length(idsElig)
  nDeaths <- nDeathsH <- 0

  if (nElig > 0) {

    whole_ages_of_elig <- pmin(ceiling(age[idsElig]), 86)
    death_rates_of_elig <- mort.rates[whole_ages_of_elig]

    idsElig.inf <- which(status[idsElig] == "h")
    death_rates_of_elig[idsElig.inf] <- death_rates_of_elig[idsElig.inf] *
                                        mort.dis.mult

    vecDeaths <- which(rbinom(nElig, 1, death_rates_of_elig) == 1)
    idsDeaths <- idsElig[vecDeaths]
    nDeaths <- length(idsDeaths)
    nDeathsH <- length(intersect(idsDeaths, idsElig.inf))

    if (nDeaths > 0) {
      dat$attr$active[idsDeaths] <- 0
      inactive <- which(dat$attr$active == 0)
      dat$attr <- deleteAttr(dat$attr, inactive)
      for (i in seq_along(dat$el)) {
        dat$el[[i]] <- delete_vertices(dat$el[[i]], inactive)
      }
    }
  }

  ## Summary statistics ##
  dat <- set_epi(dat, "d.flow", at, nDeaths)
  dat <- set_epi(dat, "d.h.flow", at, nDeathsH)

  return(dat)
}


#' @rdname moduleset-ship
#' @export
offload_covid_ship <- function(dat, at) {

  ## Attributes ##
  active <- dat$attr$active
  status <- dat$attr$status
  dxStatus <- dat$attr$dxStatus
  type <- dat$attr$type

  ## Parameters ##
  exit.rate.pass <- dat$param$exit.rate.pass
  exit.rate.crew <- dat$param$exit.rate.crew

  exit.elig.status <- dat$param$exit.elig.status
  require.dx <- dat$param$exit.require.dx

  idsElig <- which(active == 1 & status %in% exit.elig.status)
  if (require.dx == TRUE) {
    idsElig <- intersect(idsElig, which(dxStatus == 2))
  }
  nElig <- length(idsElig)
  nExits <- 0

  if (nElig > 0) {
    exit.rates <- ifelse(type[idsElig] == "p", exit.rate.pass, exit.rate.crew)
    vecExits <- which(rbinom(nElig, 1, exit.rates) == 1)
    idsExits <- idsElig[vecExits]
    nExits <- length(idsExits)

    if (nExits > 0) {
      active[idsExits] <- 0
      inactive <- which(active == 0)
      dat$attr <- deleteAttr(dat$attr, inactive)
      for (i in seq_along(dat$el)) {
        dat$el[[i]] <- delete_vertices(dat$el[[i]], inactive)
      }
    }
  }

  ## Summary statistics ##
  dat$epi$exit.flow[at] <- nExits

  return(dat)
}

#' @rdname moduleset-corporate
#' @export
arrival_covid_corporate <- function(dat, at) {

  # Parameters
  a.rate   <- get_param(dat, "a.rate")

  ## Process
  num <- dat$epi$num[1]
  nNew <- rpois(1, a.rate * num)

  ## Update Attr
  if (nNew > 0) {
    dat <- setNewAttr_covid_corporate(dat, at, nNew)
  }

  # Update Networks
  if (nNew > 0) {
    for (i in seq_along(dat$el)) {
      dat$el[[i]] <- add_vertices(dat$el[[i]], nNew)
    }
  }

  ## Output
  dat <- set_epi(dat, "nNew", at, nNew)

  return(dat)
}


setNewAttr_covid_corporate <- function(dat, at, nNew) {

  dat <- append_core_attr(dat, at, nNew)

  arrival.age <- get_param(dat, "arrival.age")
  newAges <- rep(arrival.age, nNew)
  dat <- append_attr(dat, "age", newAges, nNew)

  age.breaks <- seq(0, 200, 10)
  attr_age.grp <- cut(newAges, age.breaks, labels = FALSE, right = FALSE)
  dat <- append_attr(dat, "age.grp", attr_age.grp, nNew)

  # Disease status and related
  dat <- append_attr(dat, "status", "s", nNew)
  dat <- append_attr(dat, "infTime", NA, nNew)

  dat <- append_attr(dat, "statusTime", 0, nNew)
  dat <- append_attr(dat, "clinical", NA, nNew)
  dat <- append_attr(dat, "hospit", NA, nNew)
  dat <- append_attr(dat, "dxStatus", NA, nNew)
  dat <- append_attr(dat, "vax", 0, nNew)
  dat <- append_attr(dat, "vax1Time", NA, nNew)
  dat <- append_attr(dat, "vax2Time", NA, nNew)

  return(dat)
}
