#' @rdname moduleset-gmc19
#' @export
deaths_covid_gmc19 <- function(dat, at) {
  
  if (at>2) browser()

  ## Input
  # Attributes
  active <- get_attr(dat, "active")
  age <- get_attr(dat, "age")
  status <- get_attr(dat, "status")

  # age <- floor(age) # from HIV

  ## Parameters ##
  mort.rates <-  get_param(dat, "mort.rates")
  mort.dis.mult <- get_param(dat, "mort.dis.mult")

  ## General deaths
  idsElig <-which(as.logical(active))

  nElig <- length(idsElig)
  nDeaths <-
    nDeathsH <- 0

  if (nElig > 0) {

    whole_ages_of_elig <- pmin(ceiling(age[idsElig]), 86) # take upper bound of numeric age as integer age
    death_rates_of_elig <- mort.rates[whole_ages_of_elig] # age-dependent death rate

    idsElig.inf <- which(status[idsElig] == "h") # if hospitalized, scale up the mortality rate
    death_rates_of_elig[idsElig.inf] <- death_rates_of_elig[idsElig.inf] *mort.dis.mult

    vecDeaths <- which(rbinom(nElig, 1, death_rates_of_elig) == 1) # among the number of active nodes, use death rate of each individual to determine who dies
    idsDeaths <- idsElig[vecDeaths] # id of individual who dies
    nDeaths <- length(idsDeaths) # number of people who dies
    nDeathsH <- length(intersect(idsDeaths, idsElig.inf)) # number of people who dies and is hospitalized

    if (nDeaths > 0) { # record the ids of those who dies at this time step
      # dat$attr$active[idsDeaths] <- 0, this doesn't work
      #inactive <- which(dat$attr$active == 0) this doesn't work
      dat <- depart_nodes(dat, departures = idsDeaths) # only the things in "run" are updated
      attr.length <- unique(vapply(get_attr_list(dat), length, numeric(1)))
      if (attr.length != attributes(dat$run$el[[1]])[["n"]]) {
        stop("mismatch between el and attr length in departures mod")
      }



    }
  }

  ## Summary Output
  dat <- set_epi(dat, "d.flow", at, nDeaths)
  dat <- set_epi(dat, "d.h.flow", at, nDeathsH)

  return(dat)
}
