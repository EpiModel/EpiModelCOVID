#' @rdname moduleset-vaxDecisions
#' @export
init_covid_vax_decisions <- function(x, param, init, control, s) {

  # Main Data List Setup
  dat <- create_dat_object(param, init, control)

  dat <- init_nets(dat, x)

  # simulate first time step
  dat <- sim_nets_t1(dat)
  dat <- summary_nets(dat, at = 1L)
  
  # Add household network edgelist
  dat$num.nw <- dat$num.nw + 1
  dat$el[[dat$num.nw]] <- as.matrix(dat$param$hh.pairs)
  attr(dat$el[[dat$num.nw]], 'n') <- get_epi(dat, "sim.num", at = 1)
  dat$net_attr[[dat$num.nw]] <- list(n = get_epi(dat, "sim.num", at = 1))
  dat$control[["tergmLite.track.duration"]][[dat$num.nw]] <- FALSE

  ## Infection Status and Time Modules
  dat <- init_status_covid_vax_decisions(dat)

  ## Get initial prevalence
  dat <- prevalence_covid_vax_decisions(dat, at = 1)

  return(dat)
}


init_status_covid_vax_decisions <- function(dat) {

  e.num <- get_init(dat, "e.num")

  active <- get_attr(dat, "active")
  num <- sum(active)

  ## Disease status
  status <- rep("s", num)
  if (e.num > 0) {
    status[sample(which(active == 1), size = e.num)] <- "e"
  }

  dat <- set_attr(dat, "status", status)
  
  # Age group for vaccination processes
  age <- get_attr(dat, "age")
  
  vax.age.group <- rep(NA, length(age))
  vax.age.group[age < 5] <- 1
  vax.age.group[age >= 5 & age < 18] <- 2
  vax.age.group[age >= 18 & age < 50] <- 3
  vax.age.group[age >= 50 & age < 65] <- 4
  vax.age.group[age >= 65] <- 5
  
  dat <- set_attr(dat, "vax.age.group", vax.age.group)
  
  ## Vaccine willingness vs. resistance
  vax.willing.prob <- get_param(dat, "vax.willing.prob")
  vaxType <- rep(NA, num)
  idsAdults <- which(vax.age.group >= 3)
  vaxType.adults <- rbinom(length(idsAdults), 1, 
                           vax.willing.prob[vax.age.group[idsAdults] - 2])
  vaxType[idsAdults] <- vaxType.adults
  dat <- set_attr(dat, "vaxType", vaxType)

  # Infection Time and related attributes
  idsInf <- which(status == "e")
  infTime <- rep(NA, num)
  symptStartTime <- rep(NA, num)
  clinical <- rep(NA, num)
  hospit <- rep(NA, num)
  statusTime <- rep(NA, num)
  statusTime[idsInf] <- 1
  dxStatus <- rep(0, num)
  dxTime <- rep(NA, num)
  vax <- rep(0, num)
  vax1Time <- rep(NA, num)
  vax2Time <- rep(NA, num)
  vax3Time <- rep(NA, num)
  vax4Time <- rep(NA, num)
  vaxSE <- rep(NA, num)

  dat <- set_attr(dat, "statusTime", statusTime)
  dat <- set_attr(dat, "infTime", infTime)
  dat <- set_attr(dat, "symptStartTime", symptStartTime)
  dat <- set_attr(dat, "clinical", clinical)
  dat <- set_attr(dat, "hospit", hospit)
  dat <- set_attr(dat, "dxStatus", dxStatus)
  dat <- set_attr(dat, "dxTime", dxTime)
  dat <- set_attr(dat, "vax", vax)
  dat <- set_attr(dat, "vax1Time", vax1Time)
  dat <- set_attr(dat, "vax2Time", vax2Time)
  dat <- set_attr(dat, "vax3Time", vax3Time)
  dat <- set_attr(dat, "vax4Time", vax4Time)
  dat <- set_attr(dat, "vaxSE", vaxSE)

  return(dat)
}
