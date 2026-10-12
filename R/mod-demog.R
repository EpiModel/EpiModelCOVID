
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
  active <- get_attr(dat, "active")
  age <- get_attr(dat, "age")
  status <- get_attr(dat, "status")

  ## Parameters ##
  mort.rates <- get_param(dat, "mort.rates")
  mort.dis.mult <- get_param(dat, "mort.dis.mult")

  idsElig <- which(active == 1)
  nElig <- length(idsElig)
  nDeaths <- nDeathsIC <- 0

  if (nElig > 0) {

    whole_ages_of_elig <- pmin(ceiling(age[idsElig]), 86)
    death_rates_of_elig <- mort.rates[whole_ages_of_elig]

    idsElig.inf <- which(status[idsElig] == "ic")
    death_rates_of_elig[idsElig.inf] <- death_rates_of_elig[idsElig.inf] * mort.dis.mult

    vecDeaths <- which(rbinom(nElig, 1, death_rates_of_elig) == 1)
    idsDeaths <- idsElig[vecDeaths]
    nDeaths <- length(idsDeaths)
    nDeathsIC <- length(intersect(idsDeaths, idsElig.inf))

    if (nDeaths > 0) {
      active[idsDeaths] <- 0
      dat <- set_attr(dat, "active", active)
      inactive <- which(active == 0)
      dat <- depart_nodes(dat, inactive)
    }
  }

  ## Summary statistics ##
  dat <- set_epi(dat, "d.flow", nDeaths, at = at)
  dat <- set_epi(dat, "d.ic.flow", nDeathsIC, at = at)

  return(dat)
}


#' @rdname moduleset-corporate
#' @export
deaths_covid_corporate <- function(dat, at) {

  ## Attributes ##
  active <- get_attr(dat, "active")
  age <- get_attr(dat, "age")
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
      active[idsDeaths] <- 0
      dat <- set_attr(dat, "active", active)
      inactive <- which(active == 0)
      dat <- depart_nodes(dat, inactive)
    }
  }

  ## Summary statistics ##
  dat <- set_epi(dat, "d.flow", nDeaths, at = at)
  dat <- set_epi(dat, "d.h.flow", nDeathsH, at = at)

  return(dat)
}


#' @rdname moduleset-ship
#' @export
offload_covid_ship <- function(dat, at) {

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  dxStatus <- get_attr(dat, "dxStatus")
  type <- get_attr(dat, "type")

  ## Parameters ##
  exit.rate.pass <- get_param(dat, "exit.rate.pass")
  exit.rate.crew <- get_param(dat, "exit.rate.crew")

  exit.elig.status <- get_param(dat, "exit.elig.status")
  require.dx <- get_param(dat, "exit.require.dx")

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
      dat <- set_attr(dat, "active", active)
      inactive <- which(active == 0)
      dat <- depart_nodes(dat, inactive)
    }
  }

  ## Summary statistics ##
  dat <- set_epi(dat, "exit.flow", nExits, at = at)

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
  dat <- arrive_nodes(dat, nNew)

  ## Output
  dat <- set_epi(dat, "nNew", nNew, at = at)

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
