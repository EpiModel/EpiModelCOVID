
#' @rdname moduleset-common
#' @export
aging_covid <- function(dat, at) {

  age <- get_attr(dat, "age")
  vax.age.group <- get_attr(dat, "vax.age.group")

  #Update age
  age <- age + 1 / 365

  #Update age.grp
  age.breaks <- seq(0, 200, 10)
  age.grp <- cut(age, age.breaks, labels = FALSE, right = FALSE)

  #Update vax.age.group
  vax.age.group[which(age >= 5 & age < 18)] <- 2
  vax.age.group[which(age >= 18 & age < 50)] <- 3
  vax.age.group[which(age >= 50 & age < 65)] <- 4
  vax.age.group[which(age >= 65)] <- 5

  dat <- set_attr(dat, "vax.age.group", vax.age.group)
  dat <- set_attr(dat, "age", age)
  dat <- set_attr(dat, "age.grp", age.grp)

  return(dat)
}


#' @rdname moduleset-corporate
#' @export
deaths_covid_corporate <- function(dat, at) {

  ## Attributes ##
  active <- dat$attr$active
  age <- dat$attr$age
  status <- dat$attr$status

  ## Parameters ##
  mort.rates <- dat$param$mort.rates
  mort.dis.mult <- dat$param$mort.dis.mult

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
      dat <- depart_nodes(dat, inactive)
    }
  }

  ## Summary statistics ##
  dat <- set_epi(dat, "d.flow", at, nDeaths)
  dat <- set_epi(dat, "d.h.flow", at, nDeathsH)

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
  dat <- set_epi(dat, "nNew", at, nNew)

  return(dat)
}


setNewAttr_covid_corporate <- function(dat, at, nNew) {

  dat <- append_core_attr(dat, at, nNew)

  # Age related
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

  # Assign new nodes to households
  age.grp <- get_attr(dat, "age.grp")
  household <- get_attr(dat, "household")

  newHH <- rep(NA, length(newAges))
  for (i in seq_along(unique(attr_age.grp))) {
    newHH[attr_age.grp == unique(attr_age.grp)[i]] <-
      sample(household[which(age.grp == unique(attr_age.grp)[i])],
             sum(attr_age.grp == unique(attr_age.grp)[i]), replace = TRUE)
  }

  # Update household edgelist
  heads <- cbind((length(household) + 1):(length(household) + nNew), newHH)
  tails <- cbind(which(household %in% newHH), household[which(household %in% newHH)])
  new.edges <- merge(heads, tails, by.x = 2, by.y = 2)[, 2:3]
  new.edgelist <- as.matrix(rbind(dat$el[[dat$num.nw]], setNames(new.edges, c(".head", ".tail"))))
  attr(new.edgelist, 'n') <- attr(dat$el[[dat$num.nw]], 'n')
  dat$el[[dat$num.nw]] <- new.edgelist

  dat <- append_attr(dat, "household", newHH, nNew)

  # Disease status and related
  dat <- append_attr(dat, "status", "s", nNew)
  dat <- append_attr(dat, "infTime", NA, nNew)

  dat <- append_attr(dat, "statusTime", 0, nNew)
  dat <- append_attr(dat, "clinical", NA, nNew)
  dat <- append_attr(dat, "hospit", NA, nNew)
  dat <- append_attr(dat, "dxStatus", NA, nNew)
  dat <- append_attr(dat, "dxTime", NA, nNew)
  dat <- append_attr(dat, "vax", 0, nNew)
  dat <- append_attr(dat, "vax1Time", NA, nNew)
  dat <- append_attr(dat, "vax2Time", NA, nNew)
  dat <- append_attr(dat, "vax3Time", NA, nNew)
  dat <- append_attr(dat, "vax4Time", NA, nNew)
  dat <- append_attr(dat, "isolate", NA, nNew)
  dat <- append_attr(dat, "isoTime", NA, nNew)
  dat <- append_attr(dat, "non.office", 1, nNew)

  return(dat)
}
