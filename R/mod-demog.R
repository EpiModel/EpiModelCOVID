
#' @rdname moduleset-common
#' @export
aging_covid <- function(dat, at) {

  age <- get_attr(dat, "age")
  vaxType <- get_attr(dat, "vaxType")
  
  age <- age + 1 / 365
  
  idsNewAdults <- which(age >= 18 & is.na(vaxType))
  vax.willing.prob <- get_param(dat, "vax.willing.prob")
  vaxType.new <- rbinom(length(idsNewAdults), 1, vax.willing.prob[1])
  vaxType[idsNewAdults] <- vaxType.new
  
  dat <- set_attr(dat, "vaxType", vaxType)
  dat <- set_attr(dat, "age", age)

  return(dat)
}

#' @rdname moduleset-vaxDecisions
#' @export
deaths_covid_vax_decisions <- function(dat, at) {

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

#' @rdname moduleset-vaxDecisions
#' @export
arrival_covid_vax_decisions <- function(dat, at) {

  # Parameters
  a.rate   <- get_param(dat, "a.rate")

  ## Process
  num <- dat$epi$num[1]
  nNew <- rpois(1, a.rate * num)

  ## Update Attr
  if (nNew > 0) {
    dat <- setNewAttr_covid_vax_decisions(dat, at, nNew)
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


setNewAttr_covid_vax_decisions <- function(dat, at, nNew) {
  
  dat <- append_core_attr(dat, at, nNew)

  arrival.age <- get_param(dat, "arrival.age")
  newAges <- rep(arrival.age, nNew)
  dat <- append_attr(dat, "age", newAges, nNew)

  age.breaks <- seq(0, 200, 10)
  attr_age.grp <- cut(newAges, age.breaks, labels = FALSE, right = FALSE)
  dat <- append_attr(dat, "age.grp", attr_age.grp, nNew)
  
  vax.age.group <- rep(NA, length(newAges))
  vax.age.group[newAges < 5] <- 1
  vax.age.group[newAges >= 5 & newAges < 18] <- 2
  vax.age.group[newAges >= 18 & newAges < 50] <- 3
  vax.age.group[newAges >= 50 & newAges < 65] <- 4
  vax.age.group[newAges >= 65] <- 5
  
  dat <- append_attr(dat, "vax.age.group", vax.age.group, nNew)

  # Disease status and related
  dat <- append_attr(dat, "status", "s", nNew)
  dat <- append_attr(dat, "infTime", NA, nNew)

  dat <- append_attr(dat, "statusTime", 0, nNew)
  dat <- append_attr(dat, "symptStartTime", NA, nNew)
  dat <- append_attr(dat, "clinical", NA, nNew)
  dat <- append_attr(dat, "hospit", NA, nNew)
  dat <- append_attr(dat, "dxStatus", 0, nNew)
  dat <- append_attr(dat, "dxTime", NA, nNew)
  dat <- append_attr(dat, "vax", 0, nNew)
  dat <- append_attr(dat, "vax1Time", NA, nNew)
  dat <- append_attr(dat, "vax2Time", NA, nNew)
  dat <- append_attr(dat, "vax3Time", NA, nNew)
  dat <- append_attr(dat, "vax4Time", NA, nNew)
  dat <- append_attr(dat, "vaxSE", NA, nNew)
  
  ## Vaccine willingness vs. resistance
  vax.willing.prob <- get_param(dat, "vax.willing.prob")
  vaxType.new <- rep(NA, nNew)
  idsnewAdults <- which(vax.age.group >= 3)
  
  vaxType.newAdults <- rbinom(length(idsnewAdults), 1, 
                              vax.willing.prob[vax.age.group[idsnewAdults] - 2])
  vaxType.new[idsnewAdults] <- vaxType.newAdults
  dat <- append_attr(dat, "vaxType", vaxType.new, nNew)

  return(dat)
}
