
#' @rdname moduleset-vaxDecisions
#' @export
progress_covid_vax_decisions <- function(dat, at) {

  ## Attributes
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  statusTime <- get_attr(dat, "statusTime")
  clinical <- get_attr(dat, "clinical")
  hospit <- get_attr(dat, "hospit")
  age <- get_attr(dat, "age")
  vax <- get_attr(dat, "vax")
  symptStartTime <- get_attr(dat, "symptStartTime")
  dxStatus <- get_attr(dat, "dxStatus")
  vax1Time <- get_attr(dat, "vax1Time")
  vax2Time <- get_attr(dat, "vax2Time")
  vax3Time <- get_attr(dat, "vax3Time")
  vax4Time <- get_attr(dat, "vax4Time")
  

  ## Parameters
  prop.clinical <- get_param(dat, "prop.clinical")
  vax1.rr.clinical <- get_param(dat, "vax1.rr.clinical")
  vax2.rr.clinical <- get_param(dat, "vax2.rr.clinical")
  vax3.rr.clinical <- get_param(dat, "vax3.rr.clinical")
  vax4.rr.clinical <- get_param(dat, "vax4.rr.clinical")
  prop.hospit <- get_param(dat, "prop.hospit")
  vax1.rr.hosp <- get_param(dat, "vax1.rr.hosp")
  vax2.rr.hosp <- get_param(dat, "vax2.rr.hosp")
  vax3.rr.hosp <- get_param(dat, "vax3.rr.hosp")
  vax4.rr.hosp <- get_param(dat, "vax4.rr.hosp")
  ea.rate <- get_param(dat, "ea.rate")
  ar.rate <- get_param(dat, "ar.rate")
  eip.rate <- get_param(dat, "eip.rate")
  ipic.rate <- get_param(dat, "ipic.rate")
  ich.rate <- get_param(dat, "ich.rate")
  icr.rate <- get_param(dat, "icr.rate")
  hr.rate <- get_param(dat, "hr.rate")
  rs.rate <- get_param(dat, "rs.rate")
  half.life <- get_param(dat, "half.life")

  ## Determine Subclinical (E to A) or Clinical (E to Ip to Ic) pathway
  ids.newInf <- which(active == 1 & status == "e" & statusTime <= at & is.na(clinical))
  num.newInf <- length(ids.newInf)
  if (num.newInf > 0) {
    age.group <- pmin((floor(age[ids.newInf] / 10)) + 1, 9)
    prop.clin.vec <- prop.clinical[age.group]
    
    #vaccination reduces risk of clinical disease
    prop.clin.vec[vax[ids.newInf] == 1] <- 
      prop.clin.vec[vax[ids.newInf] == 1] * vax1.rr.clinical
    prop.clin.vec[vax[ids.newInf] == 2] <- 
      prop.clin.vec[vax[ids.newInf] == 2] * vax2.rr.clinical
    prop.clin.vec[vax[ids.newInf] == 3] <- 
      prop.clin.vec[vax[ids.newInf] == 3] * vax3.rr.clinical
    prop.clin.vec[vax[ids.newInf] == 4] <- 
      prop.clin.vec[vax[ids.newInf] == 4] * vax4.rr.clinical
    
    #waning immunity from vaccination
    sinceVax1 <- at - vax1Time[ids.newInf]
    sinceVax2 <- at - vax2Time[ids.newInf]
    sinceVax3 <- at - vax3Time[ids.newInf]
    sinceVax4 <- at - vax4Time[ids.newInf]
    
    latest.vax <- pmin(sinceVax1, sinceVax2, sinceVax3, sinceVax4, na.rm = TRUE)
    latest.vax[is.na(latest.vax)] <- 0
    
    prop.clin.vec <- pmin(prop.clin.vec * (2 ^ (latest.vax / half.life)), prop.clinical[age.group])
    
    #assign pathway
    if (any(is.na(prop.clin.vec))) stop("error in prop.clin.vec")
    vec.new.clinical <- rbinom(num.newInf, 1, prop.clin.vec)
    clinical[ids.newInf] <- vec.new.clinical
  }

  ## Subclinical Pathway
  # E to A: latent move to asymptomatic infectious
  num.new.EtoA <- 0
  ids.Es <- which(active == 1 & status == "e" & statusTime < at & clinical == 0)
  num.Es <- length(ids.Es)
  if (num.Es > 0) {
    vec.new.A <- which(rbinom(num.Es, 1, ea.rate) == 1)
    if (length(vec.new.A) > 0) {
      ids.new.A <- ids.Es[vec.new.A]
      num.new.EtoA <- length(ids.new.A)
      status[ids.new.A] <- "a"
      statusTime[ids.new.A] <- at
    }
  }

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
      symptStartTime[ids.new.Ic] <- at
    }
  }

  ## Determine Hospitalized (Ic to H) or Recovered (Ic to R) pathway
  ids.newIc <- which(active == 1 & status == "ic" & statusTime <= at & is.na(hospit))
  num.newIc <- length(ids.newIc)
  if (num.newIc > 0) {
    age.group <- pmin((floor(age[ids.newIc] / 10)) + 1, 9)
    prop.hosp.vec <- prop.hospit[age.group]
    
    #Vaccination reduces risk of hospitalization
    prop.hosp.vec[vax[ids.newIc] == 1] <- prop.hosp.vec[vax[ids.newIc] == 1] * vax1.rr.hosp
    prop.hosp.vec[vax[ids.newIc] == 2] <- prop.hosp.vec[vax[ids.newIc] == 2] * vax2.rr.hosp
    prop.hosp.vec[vax[ids.newIc] == 3] <- prop.hosp.vec[vax[ids.newIc] == 3] * vax3.rr.hosp
    prop.hosp.vec[vax[ids.newIc] == 4] <- prop.hosp.vec[vax[ids.newIc] == 4] * vax4.rr.hosp
    
    #Waning immunity from vaccination
    sinceVax1 <- at - vax1Time[ids.newIc]
    sinceVax2 <- at - vax2Time[ids.newIc]
    sinceVax3 <- at - vax3Time[ids.newIc]
    sinceVax4 <- at - vax4Time[ids.newIc]
    
    latest.vax <- pmin(sinceVax1, sinceVax2, sinceVax3, sinceVax4, na.rm = TRUE)
    latest.vax[is.na(latest.vax)] <- 0
    
    prop.hosp.vec <- pmin(prop.hosp.vec * (2 ^ (latest.vax / half.life)), prop.hospit[age.group])
    
    #Set pathway
    if (any(is.na(prop.hosp.vec))) stop("error in prop.hosp.vec")
    vec.new.hospit <- rbinom(num.newIc, 1, prop.hosp.vec)
    hospit[ids.newIc] <- vec.new.hospit
  }

  # Ic to H: clinical infectious move to hospitalized
  ids.Ich <- which(active == 1 & status == "ic" & statusTime < at & hospit == 1)
  num.Ich <- length(ids.Ich)
  if (num.Ich > 0) {
    vec.new.H <- which(rbinom(num.Ich, 1, ich.rate) == 1)
    if (length(vec.new.H) > 0) {
      ids.new.H <- ids.Ich[vec.new.H]
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
  
  #R to S: become susceptible again after infection
  num.new.RtoS <- 0
  ids.RS <- which(active == 1 & status == "r" & statusTime < at)
  num.RS <- length(ids.RS)
  if (num.RS > 0) {
    vec.new.RS <- which(rbinom(num.RS, 1, rs.rate) == 1)
    if (length(vec.new.RS) > 0) {
      ids.new.RS <- ids.RS[vec.new.RS]
      num.new.RtoS <- length(ids.new.RS)
      status[ids.new.RS] <- "s"
      statusTime[ids.new.RS] <- at
      dxStatus[ids.new.RS] <- 0
    }
  }

  ## Save updated attributes
  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "statusTime", statusTime)
  dat <- set_attr(dat, "clinical", clinical)
  dat <- set_attr(dat, "hospit", hospit)
  dat <- set_attr(dat, "symptStartTime", symptStartTime)
  dat <- set_attr(dat, "dxStatus", dxStatus)

  ## Save summary statistics
  dat <- set_epi(dat, "ea.flow", at, num.new.EtoA)
  dat <- set_epi(dat, "ar.flow", at, num.new.AtoR)
  dat <- set_epi(dat, "eip.flow", at, num.new.EtoIp)
  dat <- set_epi(dat, "ipic.flow", at, num.new.IptoIc)
  dat <- set_epi(dat, "icr.flow", at, num.new.IctoR)
  dat <- set_epi(dat, "hr.flow", at, num.new.HtoR)
  dat <- set_epi(dat, "rs.flow", at, num.new.RtoS)

  return(dat)
}


