
#' @rdname moduleset-common
#' @export
aging_covid <- function(dat, at) {

  age <- get_attr(dat, "age")
  age <- age + 1 / 365
  dat <- set_attr(dat, "age", age)

  return(dat)
}

#' @rdname moduleset-netjail
#' @export
deaths_covid_netjail <- function(dat, at) {

  ## Attributes ##
  active <- get_attr(dat, "active")
  age <- get_attr(dat, "age")
  race <- get_attr(dat, "race")
  status <- get_attr(dat, "status")

  ## Parameters ##
  mort.rates.w <- get_param(dat, "mort.rates.w")
  mort.rates.b <- get_param(dat, "mort.rates.b")
  mort.rates.o <- get_param(dat, "mort.rates.o")
  mort.dis.mult <- get_param(dat, "mort.dis.mult")

  idsElig_w <- which(active == 1 & race == "white")
  idsElig_b <- which(active == 1 & race == "black")
  idsElig_o <- which(active == 1 & race == "other")
  idsElig <- c(idsElig_w, idsElig_b, idsElig_o)
  nElig_w <- length(idsElig_w)
  nElig_b <- length(idsElig_b)
  nElig_o <- length(idsElig_o)
  nElig <- length(idsElig)
  nDeaths <- nDeathsH <- 0
  nExits <- 0

  if (nElig > 0) {

    whole_ages_of_elig_w <- pmin(ceiling(age[idsElig_w]), 100)
    whole_ages_of_elig_b <- pmin(ceiling(age[idsElig_b]), 100)
    whole_ages_of_elig_o <- pmin(ceiling(age[idsElig_o]), 100)
    death_rates_of_elig_w <- mort.rates.w[whole_ages_of_elig_w]
    death_rates_of_elig_b <- mort.rates.b[whole_ages_of_elig_b]
    death_rates_of_elig_o <- mort.rates.o[whole_ages_of_elig_o]
    
    idsElig.inf.w <- which(status[idsElig_w] == "h")
    idsElig.inf.b <- which(status[idsElig_b] == "h")
    idsElig.inf.o <- which(status[idsElig_o] == "h")
    idsElig.inf <- c(idsElig.inf.w, idsElig.inf.b, idsElig.inf.o)
    
    death_rates_of_elig_w[idsElig.inf.w] <- death_rates_of_elig_w[idsElig.inf.w] * mort.dis.mult
    death_rates_of_elig_b[idsElig.inf.b] <- death_rates_of_elig_b[idsElig.inf.b] * mort.dis.mult
    death_rates_of_elig_o[idsElig.inf.o] <- death_rates_of_elig_o[idsElig.inf.o] * mort.dis.mult
    death_rates_of_elig <- c(death_rates_of_elig_w, death_rates_of_elig_b, death_rates_of_elig_o)
    
    vecDeaths.w <- which(rbinom(nElig_w, 1, death_rates_of_elig_w) == 1)
    vecDeaths.b <- which(rbinom(nElig_b, 1, death_rates_of_elig_b) == 1)
    vecDeaths.o <- which(rbinom(nElig_o, 1, death_rates_of_elig_o) == 1)
    vecDeaths <- c(vecDeaths.w, vecDeaths.b, vecDeaths.o)
    
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


#' @rdname moduleset-netjail
#' @export
covid_release_netjail <- function(dat, at) {

  ## Attributes ##
  active <- dat$attr$active
  status <- dat$attr$status
  dxStatus <- dat$attr$dxStatus
  type <- dat$attr$type

  ## Parameters ##
  jail.exit.rate <- get_param(dat, "jail.exit.rate")

  idsElig_exit <- which(active == 1)
  nElig_exit <- length(idsElig_exit)
  nExits <- 0

  if (nElig_exit > 0) {
    vecExits <- which(rbinom(nElig_exit, 1, jail.exit.rate) == 1)
    idsExits <- idsElig_exit[vecExits]
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

#' @rdname moduleset-netjail
#' @export
arrival_covid_netjail <- function(dat, at) {

  # Parameters
  a.rate   <- get_param(dat, "a.rate")

  ## Process
  num <- dat$epi$num[1]
  nNew <- rpois(1, a.rate * num)

  ## Update Attr
  if (nNew > 0) {
    dat <- setNewAttr_covid_netjail(dat, at, nNew)
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


setNewAttr_covid_netjail <- function(dat, at, nNew) {

  dat <- append_core_attr(dat, at, nNew)
  
  newIds <- which(dat$attr$entrTime == at)
  
  ages.list <- c(dat$attr$age)
  prob.age <- proportions(dat$attr$age)
  arrival.ages <- sample(ages.list, nNew, prob = prob.age, replace = TRUE)
  dat <- append_attr(dat, "age", arrival.ages, nNew)

  age.breaks <- c(17,20,30,40,50,60,85)
  attr_age.grp <- cut(arrival.ages, age.breaks, labels = FALSE, right = FALSE)
  dat <- append_attr(dat, "age.grp", attr_age.grp, nNew)
  
  race_fact <- c("white", "black", "other")
  attr_race <- sample(race_fact, nNew, prob = c(0.093, 0.897, 0.01), replace = TRUE)
  dat <- append_attr(dat, "race", attr_race, nNew)
  
  status.list <- c("s", "e")
  prob.status <- c(0.95, 0.05)
  arrival.status <- sample(status.list, nNew, prob = prob.status, replace = TRUE)
  dat <- append_attr(dat, "status", arrival.status, nNew)
  
  newIdsInf <- which(arrival.status == "e")
  nNewInf <- length(newIdsInf)
  
  vax.list <- c(0, 3)
  prob.vax <- c(0.5, 0.5)
  vax.status <- sample(vax.list, nNew, prob = prob.vax, replace = TRUE)
  dat <- append_attr(dat, "vax", vax.status, nNew)
  
  newIdsVax <- which(vax.status == 3)
  nNewVax <- length(newIdsVax)
  
  # Disease status and related
  dat <- append_attr(dat, "infTime", NA, nNew)
  dat <- append_attr(dat, "statusTime", 0, nNew)
  dat <- append_attr(dat, "statusTime.Ic", NA, nNew)
  dat <- append_attr(dat, "clinical", NA, nNew)
  dat <- append_attr(dat, "hospit", NA, nNew)
  dat <- append_attr(dat, "dxStatus", NA, nNew)
  dat <- append_attr(dat, "dxTime", NA, nNew)
  dat <- append_attr(dat, "eligible.case", NA, nNew)
  dat <- append_attr(dat, "eligible.cc", NA, nNew)
  dat <- append_attr(dat, "symendTime", NA, nNew)
  dat <- append_attr(dat, "traced.cc", NA, nNew)
  dat <- append_attr(dat, "quar", NA, nNew)
  dat <- append_attr(dat, "tracedTime", NA, nNew)
  dat <- append_attr(dat, "quarEnd", NA, nNew)
  dat <- append_attr(dat, "iso.end", NA, nNew)
  dat <- append_attr(dat, "vaxTime", NA, nNew)

  # set attributes for infected individuals brought into the jail
  dat <- set_attr(dat, "infTime", at, newIdsInf)
  dat <- set_attr(dat, "statusTime", at, newIdsInf)
  
  # Summary characteristic for infected individuals brought into the jail
  dat <- set_epi(dat, "se.flow.incarcerated", at, nNewInf)
  dat <- set_epi(dat, "new.Vax", at, nNewVax)
  
  return(dat)
}
