
#' @rdname moduleset-common
#' @export
progress_covid <- function(dat, at) {

  ## Attributes
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  statusTime <- get_attr(dat, "statusTime") # time step of transitions
  clinical <- get_attr(dat, "clinical")
  hospit <- get_attr(dat, "hospit")
  age <- get_attr(dat, "age")
  vax <- get_attr(dat, "vax")

  ## Parameters
  prop.clinical <- get_param(dat, "prop.clinical")
    # vector of probabilities corresponding to age - asymptomatic disease
      # decreases over age
  vax.rr.clinical <- get_param(dat, "vax.rr.clinical")
  prop.hospit <- get_param(dat, "prop.hospit")
  ea.rate <- get_param(dat, "ea.rate")
  ar.rate <- get_param(dat, "ar.rate")
  at.rate <- get_param(dat, "at.rate")
  eip.rate <- get_param(dat, "eip.rate")
  ipic.rate <- get_param(dat, "ipic.rate")
  ich.rate <- get_param(dat, "ich.rate")
  icr.rate <- get_param(dat, "icr.rate")
  icicu.rate <- get_param(dat, "icicu.rate")
  ict.rate <- get_param(dat, "ict.rate")
  hr.rate <- get_param(dat, "hr.rate")
  hicu.rate <- get_param(dat, "hicu.rate")
  icur.rate <- get_param(dat, "icur.rate")
  
  # do I need to include rates into the death compartment here?
  # add to demog script

  

  ## Determine Subclinical (E to A) or Clinical (E to Ip to Ic) pathway
  ids.newInf <- which(active == 1 & status == "e" & statusTime <= at & is.na(clinical))
  # why is statusTime here allowed to be = to at? 
  num.newInf <- length(ids.newInf) # counting new infections
  if (num.newInf > 0) {
    age.group <- pmin((round(age[ids.newInf], -1)/10) + 1, 8)
    prop.clin.vec <- prop.clinical[age.group]
    prop.clin.vec[vax[ids.newInf] == 2] <- prop.clin.vec[vax[ids.newInf] == 2] *
                                            vax.rr.clinical
    if (any(is.na(prop.clin.vec))) stop("error in prop.clin.vec")
    vec.new.clinical <- rbinom(num.newInf, 1, prop.clin.vec)
    clinical[ids.newInf] <- vec.new.clinical
  }
  # if (any(status == "e" & is.na(clinical))) browser()

  ## Subclinical Pathway
  # E to A: latent move to asymptomatic infectious
  num.new.EtoA <- 0
  ids.Es <- which(active == 1 & status == "e" & statusTime < at & clinical == 0)
  # prevention of immediate transitions across multiple states within single time step
  num.Es <- length(ids.Es)
  if (num.Es > 0) {
    vec.new.A <- which(rbinom(num.Es, 1, ea.rate) == 1)
    if (length(vec.new.A) > 0) {
      ids.new.A <- ids.Es[vec.new.A]
      num.new.EtoA <- length(ids.new.A)
      status[ids.new.A] <- "a"
      statusTime[ids.new.A] <- at # tracks the time step at which transition occurred
    }
  }

  # A to T: asymptomatic infectious move to infected with positive test
  num.new.AtoT <- 0
  ids.A 
  
  # A to R: asymptomatic infectious move to recovered
  num.new.AtoR <- 0
  ids.A <- which(active == 1 & status == "a" & statusTime < at & clinical == 0)
  num.A <- length(ids.A)
  if (num.A > 0) {
    vec.new.R <- which(rbinom(num.A, 1, ar.rate) == 1)
    if (length(vec.new.R) > 0) {
      ids.new.R <- ids.A[vec.new.R]
      num.new.AtoR <- length(ids.new.R)
      status[ids.new.R] <- "r"
      statusTime[ids.new.R] <- at
    }
  }

  ## Clinical Pathway
  # E to Ip: latent move to preclinical infectious
  num.new.EtoIp <- 0
  ids.Ec <- which(active == 1 & status == "e" & statusTime < at & clinical == 1)
  num.Ec <- length(ids.Ec)
  if (num.Ec > 0) {
    vec.new.Ip <- which(rbinom(num.Ec, 1, eip.rate) == 1)
    if (length(vec.new.Ip) > 0) {
      ids.new.Ip <- ids.Ec[vec.new.Ip]
      num.new.EtoIp <- length(ids.new.Ip)
      status[ids.new.Ip] <- "ip"
      statusTime[ids.new.Ip] <- at
    }
  }

  # Ip to Ic: preclinical infectious move to clinical infectious
  num.new.IptoIc <- 0
  ids.Ip <- which(active == 1 & status == "ip" & statusTime < at & clinical == 1)
  num.Ip <- length(ids.Ip)
  if (num.Ip > 0) {
    vec.new.Ic <- which(rbinom(num.Ip, 1, ipic.rate) == 1)
    if (length(vec.new.Ic) > 0) {
      ids.new.Ic <- ids.Ip[vec.new.Ic]
      num.new.IptoIc <- length(ids.new.Ic)
      status[ids.new.Ic] <- "ic"
      statusTime[ids.new.Ic] <- at
    }
  }

  ## Determine Hospitalized (Ic to H) or Recovered (Ic to R) pathway
  ids.newIc <- which(active == 1 & status == "ic" & statusTime <= at & is.na(hospit))
  num.newIc <- length(ids.newIc)
  if (num.newIc > 0) {
    age.group <- pmin((round(age[ids.newInf], -1)/10) + 1, 8)
    prop.hosp.vec <- prop.hospit[age.group]
    if (any(is.na(prop.hosp.vec))) stop("error in prop.clin.vec")
    vec.new.hospit <- rbinom(num.newIc, 1, prop.hospit)
    hospit[ids.newIc] <- vec.new.hospit
  }

  # Ic to H: clinical infectious move to hospitalized
  num.new.IctoH <- 0
  ids.Ich <- which(active == 1 & status == "ic" & statusTime < at & hospit == 1)
  num.Ich <- length(ids.Ich)
  if (num.Ich > 0) {
    vec.new.H <- which(rbinom(num.Ich, 1, ich.rate) == 1)
    if (length(vec.new.H) > 0) {
      ids.new.H <- ids.Ich[vec.new.H]
      num.new.IctoH <- length(ids.new.H)
      status[ids.new.H] <- "h"
      statusTime[ids.new.H] <- at
    }
  }

  # Ic to H: clinical infectious move to hospitalized
  num.new.IctoH <- 0
  ids.Ich <- which(active == 1 & status == "ic" & statusTime < at & hospit == 1)
  num.Ich <- length(ids.Ich)
  if (num.Ich > 0) {
    vec.new.H <- which(rbinom(num.Ich, 1, ich.rate) == 1)
    if (length(vec.new.H) > 0) {
      ids.new.H <- ids.Ich[vec.new.H]
      num.new.IctoH <- length(ids.new.H)
      status[ids.new.H] <- "h"
      statusTime[ids.new.H] <- at
    }
  }

  # H to R: hospitalized move to recovered
  num.new.HtoR <- 0
  ids.H <- which(active == 1 & status == "h" & statusTime < at & hospit == 1)
  num.H <- length(ids.H)
  if (num.H > 0) {
    vec.new.R <- which(rbinom(num.H, 1, hr.rate) == 1)
    if (length(vec.new.R) > 0) {
      ids.new.R <- ids.H[vec.new.R]
      num.new.HtoR <- length(ids.new.R)
      status[ids.new.R] <- "r"
      statusTime[ids.new.R] <- at
    }
  }

  # Ic to R: clinical infectious move to recovered
  num.new.IctoR <- 0
  ids.Icr <- which(active == 1 & status == "ic" & statusTime < at & hospit == 0)
  num.Icr <- length(ids.Icr)
  if (num.Icr > 0) {
    vec.new.R <- which(rbinom(num.Icr, 1, icr.rate) == 1)
    if (length(vec.new.R) > 0) {
      ids.new.R <- ids.Icr[vec.new.R]
      num.new.IctoR <- length(ids.new.R)
      status[ids.new.R] <- "r"
      statusTime[ids.new.R] <- at
    }
  }

  ## Save updated attributes
  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "statusTime", statusTime)
  dat <- set_attr(dat, "clinical", clinical)
  dat <- set_attr(dat, "hospit", hospit)

  ## Save summary statistics
  dat <- set_epi(dat, "ea.flow", at, num.new.EtoA)
  dat <- set_epi(dat, "ar.flow", at, num.new.AtoR)
  dat <- set_epi(dat, "eip.flow", at, num.new.EtoIp)
  dat <- set_epi(dat, "ipic.flow", at, num.new.IptoIc)
  dat <- set_epi(dat, "icr.flow", at, num.new.IctoR)
  dat <- set_epi(dat, "hr.flow", at, num.new.HtoR)

  return(dat)
}


progress_covid_contacttrace <- function(dat, at) {
  
  ## Attributes
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  statusTime <- get_attr(dat, "statusTime") # time step of transitions
  clinical <- get_attr(dat, "clinical")
  hospit <- get_attr(dat, "hospit")
  age <- get_attr(dat, "age")
  vax <- get_attr(dat, "vax")
  
  ## Parameters
  prop.clinical <- get_param(dat, "prop.clinical")
  # vector of probabilities corresponding to age - asymptomatic disease
  # decreases over age
  # vax.rr.clinical <- get_param(dat, "vax.rr.clinical")
  prop.hospit <- get_param(dat, "prop.hospit")
  ea.rate <- get_param(dat, "ea.rate")
  ar.rate <- get_param(dat, "ar.rate")
  at.rate <- get_param(dat, "at.rate")
  eip.rate <- get_param(dat, "eip.rate")
  ipic.rate <- get_param(dat, "ipic.rate")
  ich.rate <- get_param(dat, "ich.rate")
  icr.rate <- get_param(dat, "icr.rate")
  icicu.rate <- get_param(dat, "icicu.rate")
  ict.rate <- get_param(dat, "ict.rate")
  hr.rate <- get_param(dat, "hr.rate")
  hicu.rate <- get_param(dat, "hicu.rate")
  icur.rate <- get_param(dat, "icur.rate")
  
  # do I need to include rates into the death compartment here?
  # there is no rate from S -> E compartment? 
  
  
  ## Determine Subclinical (E to A) or Clinical (E to Ip to Ic) pathway
  ids.newInf <- which(active == 1 & status == "e" & statusTime <= at & is.na(clinical))
  # why is statusTime here allowed to be = to at? 
  num.newInf <- length(ids.newInf) # counting new infections
  if (num.newInf > 0) {
    age.group <- pmin((round(age[ids.newInf], -1)/10) + 1, 8)
    prop.clin.vec <- prop.clinical[age.group]
    # prop.clin.vec[vax[ids.newInf] == 2] <- prop.clin.vec[vax[ids.newInf] == 2] *
    #                                        vax.rr.clinical
    if (any(is.na(prop.clin.vec))) stop("error in prop.clin.vec")
    vec.new.clinical <- rbinom(num.newInf, 1, prop.clinical)
    # how does this know which prop.clinical probability to use in the draw?
    clinical[ids.newInf] <- vec.new.clinical
  }
  # if (any(status == "e" & is.na(clinical))) browser()
  
  ## Subclinical Pathway
  # E to A: latent move to asymptomatic infectious
  num.new.EtoA <- 0
  ids.Es <- which(active == 1 & status == "e" & statusTime < at & clinical == 0)
  # prevention of immediate transitions across multiple states within single time step
  num.Es <- length(ids.Es)
  if (num.Es > 0) {
    vec.new.A <- which(rbinom(num.Es, 1, ea.rate) == 1)
    if (length(vec.new.A) > 0) {
      ids.new.A <- ids.Es[vec.new.A]
      num.new.EtoA <- length(ids.new.A)
      status[ids.new.A] <- "a"
      statusTime[ids.new.A] <- at # tracks the time step at which transition occurred
    }
  }
  
  # A to T: asymptomatic infectious move to infected with positive test
  num.new.AtoT <- 0
  ids.A 
  
  # A to R: asymptomatic infectious move to recovered
  num.new.AtoR <- 0
  ids.A <- which(active == 1 & status == "a" & statusTime < at & clinical == 0)
  num.A <- length(ids.A)
  if (num.A > 0) {
    vec.new.R <- which(rbinom(num.A, 1, ar.rate) == 1)
    if (length(vec.new.R) > 0) {
      ids.new.R <- ids.A[vec.new.R]
      num.new.AtoR <- length(ids.new.R)
      status[ids.new.R] <- "r"
      statusTime[ids.new.R] <- at
    }
  }
  
  ## Clinical Pathway
  # E to Ip: latent move to preclinical infectious
  num.new.EtoIp <- 0
  ids.Ec <- which(active == 1 & status == "e" & statusTime < at & clinical == 1)
  num.Ec <- length(ids.Ec)
  if (num.Ec > 0) {
    vec.new.Ip <- which(rbinom(num.Ec, 1, eip.rate) == 1)
    if (length(vec.new.Ip) > 0) {
      ids.new.Ip <- ids.Ec[vec.new.Ip]
      num.new.EtoIp <- length(ids.new.Ip)
      status[ids.new.Ip] <- "ip"
      statusTime[ids.new.Ip] <- at
    }
  }
  
  # Ip to Ic: preclinical infectious move to clinical infectious
  num.new.IptoIc <- 0
  ids.Ip <- which(active == 1 & status == "ip" & statusTime < at & clinical == 1)
  num.Ip <- length(ids.Ip)
  if (num.Ip > 0) {
    vec.new.Ic <- which(rbinom(num.Ip, 1, ipic.rate) == 1)
    if (length(vec.new.Ic) > 0) {
      ids.new.Ic <- ids.Ip[vec.new.Ic]
      num.new.IptoIc <- length(ids.new.Ic)
      status[ids.new.Ic] <- "ic"
      statusTime[ids.new.Ic] <- at
    }
  }
  
  ## Determine Hospitalized (Ic to H) or Recovered (Ic to R) pathway
  ids.newIc <- which(active == 1 & status == "ic" & statusTime <= at & is.na(hospit))
  num.newIc <- length(ids.newIc)
  if (num.newIc > 0) {
    age.group <- pmin((round(age[ids.newInf], -1)/10) + 1, 8)
    prop.hosp.vec <- prop.hospit[age.group]
    if (any(is.na(prop.hosp.vec))) stop("error in prop.clin.vec")
    vec.new.hospit <- rbinom(num.newIc, 1, prop.hospit)
    hospit[ids.newIc] <- vec.new.hospit
  }
  
  # Ic to H: clinical infectious move to hospitalized
  num.new.IctoH <- 0
  ids.Ich <- which(active == 1 & status == "ic" & statusTime < at & hospit == 1)
  num.Ich <- length(ids.Ich)
  if (num.Ich > 0) {
    vec.new.H <- which(rbinom(num.Ich, 1, ich.rate) == 1)
    if (length(vec.new.H) > 0) {
      ids.new.H <- ids.Ich[vec.new.H]
      num.new.IctoH <- length(ids.new.H)
      status[ids.new.H] <- "h"
      statusTime[ids.new.H] <- at
    }
  }
  
  # Ic to H: clinical infectious move to hospitalized
  num.new.IctoH <- 0
  ids.Ich <- which(active == 1 & status == "ic" & statusTime < at & hospit == 1)
  num.Ich <- length(ids.Ich)
  if (num.Ich > 0) {
    vec.new.H <- which(rbinom(num.Ich, 1, ich.rate) == 1)
    if (length(vec.new.H) > 0) {
      ids.new.H <- ids.Ich[vec.new.H]
      num.new.IctoH <- length(ids.new.H)
      status[ids.new.H] <- "h"
      statusTime[ids.new.H] <- at
    }
  }
  
  # H to R: hospitalized move to recovered
  num.new.HtoR <- 0
  ids.H <- which(active == 1 & status == "h" & statusTime < at & hospit == 1)
  num.H <- length(ids.H)
  if (num.H > 0) {
    vec.new.R <- which(rbinom(num.H, 1, hr.rate) == 1)
    if (length(vec.new.R) > 0) {
      ids.new.R <- ids.H[vec.new.R]
      num.new.HtoR <- length(ids.new.R)
      status[ids.new.R] <- "r"
      statusTime[ids.new.R] <- at
    }
  }
  
  # Ic to R: clinical infectious move to recovered
  num.new.IctoR <- 0
  ids.Icr <- which(active == 1 & status == "ic" & statusTime < at & hospit == 0)
  num.Icr <- length(ids.Icr)
  if (num.Icr > 0) {
    vec.new.R <- which(rbinom(num.Icr, 1, icr.rate) == 1)
    if (length(vec.new.R) > 0) {
      ids.new.R <- ids.Icr[vec.new.R]
      num.new.IctoR <- length(ids.new.R)
      status[ids.new.R] <- "r"
      statusTime[ids.new.R] <- at
    }
  }
  
  ## Save updated attributes
  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "statusTime", statusTime)
  dat <- set_attr(dat, "clinical", clinical)
  dat <- set_attr(dat, "hospit", hospit)
  
  ## Save summary statistics
  dat <- set_epi(dat, "ea.flow", at, num.new.EtoA)
  dat <- set_epi(dat, "ar.flow", at, num.new.AtoR)
  dat <- set_epi(dat, "eip.flow", at, num.new.EtoIp)
  dat <- set_epi(dat, "ipic.flow", at, num.new.IptoIc)
  dat <- set_epi(dat, "icr.flow", at, num.new.IctoR)
  dat <- set_epi(dat, "hr.flow", at, num.new.HtoR)
  
  return(dat)
}
