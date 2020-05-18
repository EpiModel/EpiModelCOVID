
#' @rdname moduleset-ship
#' @export
aging_covid_ship <- function(dat, at) {

  dat$attr$age <- dat$attr$age + 1/365

  return(dat)
}


#' @rdname moduleset-ship
#' @export
deaths_covid_ship <- function(dat, at) {

  ## Attributes ##
  active <- dat$attr$active
  age <- dat$attr$age
  status <- dat$attr$status

  ## Parameters ##
  mort.rates <- dat$param$mort.rates
  mort.dis.mult <- dat$param$mort.dis.mult

  idsElig <- which(active == 1)
  nElig <- length(idsElig)
  nDeaths <- 0

  if (nElig > 0) {

    whole_ages_of_elig <- pmin(ceiling(age[idsElig]), 86)
    death_rates_of_elig <- mort.rates[whole_ages_of_elig]

    idsElig.inf <- which(status[idsElig] == "ic")
    death_rates_of_elig[idsElig.inf] <- death_rates_of_elig[idsElig.inf] * mort.dis.mult

    vecDeaths <- which(rbinom(nElig, 1, death_rates_of_elig) == 1)
    idsDeaths <- idsElig[vecDeaths]
    nDeaths <- length(idsDeaths)

    if (nDeaths > 0) {
      dat$attr$active[idsDeaths] <- 0
      inactive <- which(dat$attr$active == 0)
      dat$attr <- deleteAttr(dat$attr, inactive)
      for (i in 1:length(dat$el)) {
        dat$el[[i]] <- delete_vertices(dat$el[[i]], inactive)
      }
    }
  }

  ## Summary statistics ##
  dat$epi$d.flow[at] <- nDeaths

  return(dat)
}


#' @rdname moduleset-ship
#' @export
offload_covid_ship <- function(dat, at) {

  ## Attributes ##
  active <- dat$attr$active
  age <- dat$attr$age
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
      for (i in 1:length(dat$el)) {
        dat$el[[i]] <- delete_vertices(dat$el[[i]], inactive)
      }
    }
  }

  ## Summary statistics ##
  dat$epi$exit.flow[at] <- nExits

  return(dat)
}
