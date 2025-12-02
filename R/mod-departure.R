#' @rdname moduleset-gmc19
#' @export
deaths_covid_gmc19 <- function(dat, at) {
  
  #if (at>=4) browser()

  ## Input
  # Attributes
  active <- get_attr(dat, "active")
  age <- get_attr(dat, "age")
  status <- get_attr(dat, "status")

  #age <- floor(age) 

  ## Parameters ##
  mort.rates <-  get_param(dat, "mort.rates")
  #mort.dis.mult <- get_param(dat, "mort.dis.mult")

  ## General deaths
  idsElig <- which(as.logical(active))
  
  if (length(idsElig) > 0) {
  age_idx_elig <- pmin(ceiling(age[idsElig]), 86) # take upper bound of numeric age as integer age (0.5 -> 1); this define a age index look up mortality rate

  death_rates_of_elig <- mort.rates[age_idx_elig] # mortality rate of each node based on their age index

  idsDep <- idsElig[ runif(length(idsElig)) < # vector of random nunmber btw 0-1 for each node
                       death_rates_of_elig] # when the number < death rate of that node, death occur
  
    if (length(idsDep) > 0) { # record the ids of those who dies at this time step
      dat <- set_attr(dat, "active", 0, posit_ids = idsDep)
      dat <- depart_nodes(dat, departures = idsDep) 
      attr.length <- unique(vapply(get_attr_list(dat), length, numeric(1)))
      if (attr.length != attributes(dat$run$el[[1]])[["n"]]) {
        stop("mismatch between el and attr length in departures mod")
      }



  }
}
  ## Summary Output
  #dat <- set_epi(dat, "d.flow", at, nDeaths)
  #dat <- set_epi(dat, "d.h.flow", at, nDeathsH)

  return(dat)
}
