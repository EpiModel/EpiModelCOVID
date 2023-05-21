
#' @rdname moduleset-common
#' @export
progress_covid <- function(dat, at) {

  ## Attributes
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  statusTime <- get_attr(dat, "statusTime")
  clinical <- get_attr(dat, "clinical")
  hospit <- get_attr(dat, "hospit")
  age <- get_attr(dat, "age")
  vax <- get_attr(dat, "vax")
  isolate <- get_attr(dat, "isolate")
  isoTime <- get_attr(dat, "isoTime")
  vax1Time <- get_attr(dat, "vax1Time")
  vax2Time <- get_attr(dat, "vax2Time")
  vax3Time <- get_attr(dat, "vax3Time")
  dxStatus <- get_attr(dat, "dxStatus")
  dxTime <- get_attr(dat, "dxTime")

  ## Parameters
  prop.clinical <- get_param(dat, "prop.clinical")
  vax.rr.clinical <- get_param(dat, "vax.rr.clinical")
  vax1.immune <- get_param(dat, "vax1.immune")
  vax2.immune <- get_param(dat, "vax2.immune")
  vax3.immune <- get_param(dat, "vax3.immune")
  half.life <- get_param(dat, "half.life")
  prop.hospit <- get_param(dat, "prop.hospit")
  vax1.rr.hosp <- get_param(dat, "vax1.rr.hosp")
  vax2.rr.hosp <- get_param(dat, "vax2.rr.hosp")
  vax3.rr.hosp <- get_param(dat, "vax3.rr.hosp")
  ea.rate <- get_param(dat, "ea.rate")
  ar.rate <- get_param(dat, "ar.rate")
  eip.rate <- get_param(dat, "eip.rate")
  ipic.rate <- get_param(dat, "ipic.rate")
  ich.rate <- get_param(dat, "ich.rate")
  icr.rate <- get_param(dat, "icr.rate")
  hr.rate <- get_param(dat, "hr.rate")
  rs.rate <- get_param(dat, "rs.rate")
  iso.prob <- get_param(dat, "iso.prob")
  notif.prob <- get_param(dat, "notif.prob")

  # Determine exposures to infected contacts for masking guidelines
  num.elig.exp <- 0
  num.new.notif.iso <- 0
  # pull ids for non-symtomatic, not diagnosed, not isolating
  ids.elig.exp <- which(active == 1 & status %in% c("s","a","e","ip","r") &
                       is.na(isolate) & dxStatus %in% 0:1)
  num.elig.exp <- length(ids.elig.exp)
  if (num.elig.exp > 0) {
    # pull their partners and only keep those who have been diagnosed
    # in the past 3 time steps
    ids.exp <- get_partners(dat,ids.elig.exp, only.active.nodes = TRUE)
    ids.exp$index_posit_ids <- get_posit_ids(dat,ids.exp$index)
    ids.exp$partner_posit_ids <- get_posit_ids(dat,ids.exp$partner)
    ids.dx <- which(active == 1 & dxStatus == 2 & dxTime >= at-2)
    ids.dx.time <- data.frame(partner_posit_ids = ids.dx,
                              dxTime = dxTime[ids.dx])
    ids.exp.dx <- merge(ids.dx.time,ids.exp,by = "partner_posit_ids")
    # keep only contacts that occurred before the diagnosis, by network
    ids.exp.dx2 <- subset(ids.exp.dx,
                          (network == 1 & ids.exp.dx$dxTime == at-1) |
                            (network == 2 & ids.exp.dx$dxTime >= at-3) |
                            (network == 3 &
                               ids.exp.dx$dxTime >= ids.exp.dx$stop &
                               ids.exp.dx$dxTime <= ids.exp.dx$stop + 3))
    ids.index <- ids.exp.dx2$index_posit_ids
    ids.index <- unique(ids.index)
    if (length(ids.index) > 0) {
      vec.new.notif <- which(rbinom(length(ids.index),1,notif.prob) == 1)
      if (length(vec.new.notif) > 0) {
        ids.new.notif <- ids.index[vec.new.notif]
        vec.new.notif.iso <- which(rbinom(length(ids.new.notif), 1, iso.prob) == 1)
        if (length(vec.new.notif.iso) > 0) {
          ids.new.notif.iso <- ids.new.notif[vec.new.notif.iso]
          num.new.notif.iso <- length(ids.new.notif.iso)
          partner.dx.time <- subset(ids.exp.dx2, ids.exp.dx2$index_posit_ids %in% ids.new.notif.iso)
          partner.dx.time <- partner.dx.time[ave(partner.dx.time$dxTime,
                                                 partner.dx.time$index_posit_ids,
                                                 FUN = function(x) x == min(x)) == 1, ]
          # need to keep one dxtime per index
          # set exposure time based on network
          isolate[ids.new.notif.iso] <- 4 # masking due to exposure notification
          isoTime[ids.new.notif.iso] <- partner.dx.time$dxTime # change to time of exposure
        }
      }
    }
  }


  ## Determine Subclinical (E to A) or Clinical (E to Ip to Ic) pathway
  ids.newInf <- which(active == 1 & status == "e" & statusTime <= at & is.na(clinical))
  num.newInf <- length(ids.newInf)
  if (num.newInf > 0) {
    age.group <- pmin((floor(age[ids.newInf] / 10)) + 1, 8)
    prop.clin.vec <- prop.clinical[age.group]

    # Vaccination reducing risk of developing clinical disease
    prop.clin.vec[vax[ids.newInf] == 4] <- prop.clin.vec[vax[ids.newInf] == 4] *
                                           vax.rr.clinical
    if (any(is.na(prop.clin.vec))) stop("error in prop.clin.vec")

    # Waning vaccine provided protection
    sinceVax1 <- ifelse((at - vax1Time[ids.newInf] - vax1.immune) >= 0,
                        at - vax1Time[ids.newInf] - vax1.immune, 0)
    sinceVax2 <- ifelse((at - vax2Time[ids.newInf] - vax2.immune) >= 0,
                        at - vax2Time[ids.newInf] - vax2.immune, 0)
    sinceVax3 <- ifelse((at - vax3Time[ids.newInf] - vax3.immune) >= 0,
                        at - vax3Time[ids.newInf] - vax3.immune, 0)
    latest.vax <- pmin(sinceVax1, sinceVax2, sinceVax3, na.rm = TRUE)
    latest.vax[is.na(latest.vax)] <- 0

    latest.vax.newInf <- latest.vax
    prop.clin.vec <- prop.clin.vec *
      (2 ^ (latest.vax.newInf / half.life))

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
  num.new.iso1 <- 0
  ids.Ip <- which(active == 1 & status == "ip" & statusTime < at & clinical == 1)
  num.Ip <- length(ids.Ip)
  if (num.Ip > 0) {
    vec.new.Ic <- which(rbinom(num.Ip, 1, ipic.rate) == 1)
    if (length(vec.new.Ic) > 0) {
      ids.new.Ic <- ids.Ip[vec.new.Ic]
      num.new.IptoIc <- length(ids.new.Ic)
      status[ids.new.Ic] <- "ic"
      statusTime[ids.new.Ic] <- at
      vec.new.iso <- which(rbinom(length(vec.new.Ic), 1, iso.prob) == 1)
      if (length(vec.new.iso) > 0) {
        ids.new.iso1 <- ids.new.Ic[vec.new.iso]
        num.new.iso1 <- length(ids.new.iso1)
        isolate[ids.new.iso1] <- 1 # isolation for mild infection
        isoTime[ids.new.iso1] <- at
      }
    }
  }

  ## Determine Hospitalized (Ic to H) or Recovered (Ic to R) pathway
  ids.newIc <- which(active == 1 & status == "ic" & statusTime <= at & is.na(hospit))
  num.newIc <- length(ids.newIc)
  if (num.newIc > 0) {
    age.group <- pmin((floor(age[ids.newIc] / 10)) + 1, 8)
    prop.hosp.vec <- prop.hospit[age.group]
    if (any(is.na(prop.hosp.vec))) stop("error in prop.hosp.vec")

    # Vaccination reducing risk of hospitalization
    prop.hosp.vec[vax[ids.newIc] %in% 2:3] <- prop.hosp.vec[vax[ids.newIc] %in% 2:3] *
      vax1.rr.hosp
    prop.hosp.vec[vax[ids.newIc] %in% 4:5] <- prop.hosp.vec[vax[ids.newIc] %in% 4:5] *
      vax2.rr.hosp
    prop.hosp.vec[vax[ids.newIc] == 6] <- prop.hosp.vec[vax[ids.newIc] == 6] *
      vax3.rr.hosp

    # Waning vaccine provided protection
    sinceVax1 <- ifelse((at - vax1Time[ids.newIc] - vax1.immune) >= 0,
                        at - vax1Time[ids.newIc] - vax1.immune, 0)
    sinceVax2 <- ifelse((at - vax2Time[ids.newIc] - vax2.immune) >= 0,
                        at - vax2Time[ids.newIc] - vax2.immune, 0)
    sinceVax3 <- ifelse((at - vax3Time[ids.newIc] - vax3.immune) >= 0,
                        at - vax3Time[ids.newIc] - vax3.immune, 0)
    latest.vax <- pmin(sinceVax1, sinceVax2, sinceVax3, na.rm = TRUE)
    latest.vax[is.na(latest.vax)] <- 0

    latest.vax.newIc <- latest.vax
    prop.hosp.vec <- prop.hosp.vec *
      (2 ^ (latest.vax.newIc / half.life))

    vec.new.hospit <- rbinom(num.newIc, 1, prop.hosp.vec)
    hospit[ids.newIc] <- vec.new.hospit
  }

  # Ic to H: clinical infectious move to hospitalized
  num.new.iso2 <- 0
  ids.Ich <- which(active == 1 & status == "ic" & statusTime < at & hospit == 1)
  num.Ich <- length(ids.Ich)
  if (num.Ich > 0) {
    vec.new.H <- which(rbinom(num.Ich, 1, ich.rate) == 1)
    if (length(vec.new.H) > 0) {
      ids.new.H <- ids.Ich[vec.new.H]
      status[ids.new.H] <- "h"
      statusTime[ids.new.H] <- at
      ids.new.iso2 <- which(status == "h" & statusTime == at & isolate == 1)
      num.new.iso2 <- length(ids.new.iso2)
      if (num.new.iso2 > 0) {
        isolate[ids.new.iso2] <- 2 # isolation for severe infection
      }
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
      dxStatus[ids.new.RS] <- NA
    }
  }

  # Move mild isolation to masking among R
  num.new.iso.mask <- 0
  ids.new.iso.mask <- which(active == 1 & status == "r" & isolate == 1 & (at - isoTime) > 5)
  num.new.iso.mask <- length(ids.new.iso.mask)
  if (num.new.iso.mask > 0) {
    isolate[ids.new.iso.mask] <- 3 # post-isolation masking
  }

  # End isolation pathway
  num.new.iso.end <- 0
  ids.new.iso.end <- which(active == 1 &
                             ((status == "r" & isolate %in% c(2,3)) | isolate == 4) &
                             (at - isoTime) > 10)
  num.new.iso.end <- length(ids.new.iso.end)
  if (num.new.iso.end > 0) {
    isolate[ids.new.iso.end] <- NA # end isolation pathway
    isoTime[ids.new.iso.end] <- NA
  }


  ## Save updated attributes
  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "statusTime", statusTime)
  dat <- set_attr(dat, "clinical", clinical)
  dat <- set_attr(dat, "hospit", hospit)
  dat <- set_attr(dat, "isolate", isolate)
  dat <- set_attr(dat, "isoTime", isoTime)

  ## Save summary statistics
  dat <- set_epi(dat, "ea.flow", at, num.new.EtoA)
  dat <- set_epi(dat, "ar.flow", at, num.new.AtoR)
  dat <- set_epi(dat, "eip.flow", at, num.new.EtoIp)
  dat <- set_epi(dat, "ipic.flow", at, num.new.IptoIc)
  dat <- set_epi(dat, "icr.flow", at, num.new.IctoR)
  dat <- set_epi(dat, "hr.flow", at, num.new.HtoR)
  dat <- set_epi(dat, "rs.flow", at, num.new.RtoS)
  dat <- set_epi(dat, "iso1.flow", at, num.new.iso1)
  dat <- set_epi(dat, "iso2.flow", at, num.new.iso2)
  dat <- set_epi(dat, "isomask.flow", at, num.new.iso.mask)
  dat <- set_epi(dat, "isoend.flow", at, num.new.iso.end)

  return(dat)
}


#' @rdname moduleset-ship
#' @export
progress_covid_ship <- function(dat, at) {

  ## Attributes
  active <- dat$attr$active
  status <- dat$attr$status
  statusTime <- dat$attr$statusTime
  clinical <- dat$attr$clinical
  age <- dat$attr$age

  ## Parameters
  prop.clinical <- dat$param$prop.clinical
  ea.rate <- dat$param$ea.rate
  ar.rate <- dat$param$ar.rate
  eip.rate <- dat$param$eip.rate
  ipic.rate <- dat$param$ipic.rate
  icr.rate <- dat$param$icr.rate

  ## Determine Subclinical (E to A) or Clinical (E to Ip to Ic) pathway
  ids.newInf <- which(active == 1 & status == "e" & statusTime <= at & is.na(clinical))
  num.newInf <- length(ids.newInf)
  if (num.newInf > 0) {
    age.group <- pmin((floor(age[ids.newInf] / 10)) + 1, 8)
    prop.clin.vec <- prop.clinical[age.group]
    if (any(is.na(prop.clin.vec))) stop("error in prop.clin.vec")
    vec.new.clinical <- rbinom(num.newInf, 1, prop.clinical)
    clinical[ids.newInf] <- vec.new.clinical
  }
  if (any(status == "e" & is.na(clinical))) browser()

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
    }
  }

  # Ic to R: clinical infectious move to recovered (if not mortality first)
  num.new.IctoR <- 0
  ids.Ic <- which(active == 1 & status == "ic" & statusTime < at & clinical == 1)
  num.Ic <- length(ids.Ic)
  if (num.Ic > 0) {
    vec.new.R <- which(rbinom(num.Ic, 1, icr.rate) == 1)
    if (length(vec.new.R) > 0) {
      ids.new.R <- ids.Ic[vec.new.R]
      num.new.IctoR <- length(ids.new.R)
      status[ids.new.R] <- "r"
      statusTime[ids.new.R] <- at
    }
  }

  ## Save updated status attribute
  dat$attr$status <- status
  dat$attr$statusTime <- statusTime
  dat$attr$clinical <- clinical

  ## Save summary statistics
  dat$epi$ea.flow[at] <- num.new.EtoA
  dat$epi$ar.flow[at] <- num.new.AtoR

  dat$epi$eip.flow[at] <- num.new.EtoIp
  dat$epi$ipic.flow[at] <- num.new.IptoIc
  dat$epi$icr.flow[at] <- num.new.IctoR

  return(dat)
}
