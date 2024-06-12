
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
  vax4Time <- get_attr(dat, "vax4Time")
  dxStatus <- get_attr(dat, "dxStatus")
  dxTime <- get_attr(dat, "dxTime")
  household <- get_attr(dat, "household")
  non.office <- get_attr(dat, "non.office")

  ## Parameters
  prop.clinical <- get_param(dat, "prop.clinical")
  vax1.rr.clinical <- get_param(dat, "vax1.rr.clinical")
  vax2.rr.clinical <- get_param(dat, "vax2.rr.clinical")
  vax3.rr.clinical <- get_param(dat, "vax3.rr.clinical")
  vax4.rr.clinical <- get_param(dat, "vax4.rr.clinical")
  half.life <- get_param(dat, "half.life")
  prop.hospit <- get_param(dat, "prop.hospit")
  hosp.boost.mult <- get_param(dat, "hosp.boost.mult")
  hosp.boost.start <- get_param(dat, "hosp.boost.start")
  hosp.boost.stop <- get_param(dat, "hosp.boost.stop")
  hosp.supp.mult <- get_param(dat, "hosp.supp.mult")
  hosp.supp.start <- get_param(dat, "hosp.supp.start")
  hosp.supp.stop <- get_param(dat, "hosp.supp.stop")
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
  iso.prob.office <- get_param(dat, "iso.prob.office")
  iso.prob.other <- get_param(dat, "iso.prob.other")
  notif.prob <- get_param(dat, "notif.prob")

  # Number of networks, excluding households
  nLayers <- length(dat$el) - 1

  # Determine exposures to infected contacts for masking guidelines
  num.elig.exp <- 0
  num.new.iso4.other <- 0
  num.new.iso4.office <- 0

  # pull ids for non-symptomatic, not diagnosed, not isolating
  ids.elig.exp <- which(active == 1 & status %in% c("s","a","e","ip","r") &
                       is.na(isolate) & dxStatus %in% 0:1)
  num.elig.exp <- length(ids.elig.exp)
  if (num.elig.exp > 0) {

    # pull partners of ids identified above (only looks at work/community layers)
    ids.exp <- get_partners(dat,ids.elig.exp, only.active.nodes = TRUE)
    ids.exp$index_posit_ids <- get_posit_ids(dat,ids.exp$index)
    ids.exp$partner_posit_ids <- get_posit_ids(dat,ids.exp$partner)
    ids.exp$hh.index <- household[ids.exp$index_posit_ids]

    # identify people who have been diagnosed in the past 3 time steps
    ids.dx <- which(active == 1 & dxStatus == 2 & dxTime >= at-2)
    ids.dx.time <- data.frame(partner_posit_ids = ids.dx,
                              dxTime = dxTime[ids.dx],
                              hh.partner = household[ids.dx],
                              non.office = non.office[ids.dx])
    ids.exp.dx <- merge(ids.dx.time,ids.exp,by = "partner_posit_ids")

    # pull household members of those diagnosed in past 3 time steps
    ids.exp.dx.hh <- data.frame(
      partner_posit_ids = numeric(),
      dxTime = numeric(),
      hh.partner = numeric(),
      non.office = numeric(),
      index = numeric(),
      partner = numeric(),
      start = numeric(),
      stop = numeric(),
      network = numeric(),
      index_posit_ids = numeric(),
      hh.index = numeric()
    )
    for (i in ids.dx) {
      l <- length(which(household==household[i]))
      new_row <- list(
        partner_posit_ids = rep(i,l),
        dxTime = rep(dxTime[i],l),
        hh.partner = rep(household[i],l),
        non.office = rep(non.office[i],l),
        index = get_unique_ids(dat,which(household==household[i])),
        partner = rep(get_unique_ids(dat,i),l),
        start = rep(1,l),
        stop = rep(NA,l),
        network = rep(3,l),
        index_posit_ids = which(household==household[i]),
        hh.index = household[which(household==household[i])]
      )
      ids.exp.dx.hh <- rbind(ids.exp.dx.hh,new_row)
    }
    ids.exp.dx.hh <- ids.exp.dx.hh[ids.exp.dx.hh$index_posit_ids %in% ids.exp$index_posit_ids,]

    # add household edges
    ids.exp.dx <- rbind(ids.exp.dx,ids.exp.dx.hh)

    # keep only contacts that occurred before the diagnosis, by network
    ids.exp.dx2 <- subset(ids.exp.dx,
                          (network == 3 & ids.exp.dx$dxTime == at-1) |
                            (network == 1 & ids.exp.dx$dxTime >= at-3) |
                            (network == 2 &
                               ids.exp.dx$dxTime >= ids.exp.dx$stop &
                               ids.exp.dx$dxTime <= ids.exp.dx$stop + 3))

    if (length(unique(ids.exp.dx2$index_posit_ids)) > 0) {
      # loop through networks
      for (j in 1:3) {
        ids.index <- unique(ids.exp.dx2$index_posit_ids[ids.exp.dx2$network == j])
        
        # identify those who were notified of the exposure by their contact
        vec.new.notif <- which(rbinom(length(ids.index),1,notif.prob[j]) == 1)

        if (length(vec.new.notif) > 0) {
          ids.new.notif <- ids.index[vec.new.notif]
          
          # identify those who decide to follow isolation guidelines
          ### office workers
          ids.new.notif.office <- intersect(ids.new.notif, which(is.na(isolate) & non.office==0))
          vec.new.iso4.office <- which(rbinom(length(ids.new.notif.office), 1, iso.prob.office) == 1)
          if (length(vec.new.iso4.office) > 0) {
            ids.new.iso4.office <- ids.new.notif.office[vec.new.iso4.office]
            num.new.iso4.office <- num.new.iso4.office + length(ids.new.iso4.office)
            isolate[ids.new.iso4.office] <- 4 # masking due to exposure notification

            partner.dx.time <- subset(ids.exp.dx2, ids.exp.dx2$index_posit_ids %in% ids.new.iso4.office)
            if (j == 2) {
              # for community contacts, start isolation at time of contact
              partner.dx.time <- partner.dx.time[,c("start",
                                                    "index_posit_ids")]
              partner.dx.time <- partner.dx.time[ave(partner.dx.time$start,
                                                     partner.dx.time$index_posit_ids,
                                                     FUN = function(x) x == min(x)) == 1, ]
              partner.dx.time <- unique(partner.dx.time)
              # if (length(ids.new.iso4.office) != length(partner.dx.time$start)) browser()
              isoTime[ids.new.iso4.office] <- partner.dx.time$start
            } else {
              # for HH and office contacts, start isolation at time of diagnosis
              # since contacts are 'permanent'
              partner.dx.time <- partner.dx.time[,c("dxTime",
                                                    "index_posit_ids")]
              partner.dx.time <- partner.dx.time[ave(partner.dx.time$dxTime,
                                                     partner.dx.time$index_posit_ids,
                                                     FUN = function(x) x == min(x)) == 1, ]
              partner.dx.time <- unique(partner.dx.time)
              # if (length(ids.new.iso4.office) != length(partner.dx.time$dxTime)) browser()
              isoTime[ids.new.iso4.office] <- partner.dx.time$dxTime
            }

          }
          ### non-office workers
          ids.new.notif.other <- intersect(ids.new.notif, which(is.na(isolate) & non.office==1))
          vec.new.iso4.other <- which(rbinom(length(ids.new.notif.other), 1, iso.prob.other) == 1)
          if (length(vec.new.iso4.other) > 0) {
            ids.new.iso4.other <- ids.new.notif.other[vec.new.iso4.other]
            num.new.iso4.other <- num.new.iso4.other + length(ids.new.iso4.other)
            isolate[ids.new.iso4.other] <- 4 # masking due to exposure notification
            
            partner.dx.time <- subset(ids.exp.dx2, ids.exp.dx2$index_posit_ids %in% ids.new.iso4.other)
            if (j == 2) {
              # for community contacts, start isolation at time of contact
              partner.dx.time <- partner.dx.time[,c("start",
                                                    "index_posit_ids")]
              partner.dx.time <- partner.dx.time[ave(partner.dx.time$start,
                                                     partner.dx.time$index_posit_ids,
                                                     FUN = function(x) x == min(x)) == 1, ]
              partner.dx.time <- unique(partner.dx.time)
              # if (length(ids.new.iso4.other) != length(partner.dx.time$start)) browser()
              isoTime[ids.new.iso4.other] <- partner.dx.time$start
            } else {
              # for HH and office contacts, start isolation at time of diagnosis
              # since contacts are 'permanent'
              partner.dx.time <- partner.dx.time[,c("dxTime",
                                                    "index_posit_ids")]
              partner.dx.time <- partner.dx.time[ave(partner.dx.time$dxTime,
                                                     partner.dx.time$index_posit_ids,
                                                     FUN = function(x) x == min(x)) == 1, ]
              partner.dx.time <- unique(partner.dx.time)
              # if (length(ids.new.iso4.other) != length(partner.dx.time$dxTime)) browser()
              isoTime[ids.new.iso4.other] <- partner.dx.time$dxTime
            }
            
          }
        }
      }

    }
  }


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
    if (any(is.na(prop.clin.vec))) stop("error in prop.clin.vec")

    # Waning vaccine provided protection
    sinceVax1 <- at - vax1Time[ids.newInf]
    sinceVax2 <- at - vax2Time[ids.newInf]
    sinceVax3 <- at - vax3Time[ids.newInf]
    sinceVax4 <- at - vax4Time[ids.newInf]

    latest.vax <- pmin(sinceVax1, sinceVax2, sinceVax3, sinceVax4, na.rm = TRUE)
    latest.vax[is.na(latest.vax)] <- 0

    prop.clin.vec <- pmin(prop.clin.vec * (2 ^ (latest.vax / half.life)), prop.clinical[age.group])

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
  num.new.IptoIc.w <- 0
  num.new.iso1.other <- 0
  num.new.iso1.office <- 0
  ids.Ip <- which(active == 1 & status == "ip" & statusTime < at & clinical == 1)
  num.Ip <- length(ids.Ip)
  if (num.Ip > 0) {
    vec.new.Ic <- which(rbinom(num.Ip, 1, ipic.rate) == 1)
    if (length(vec.new.Ic) > 0) {
      ids.new.Ic <- ids.Ip[vec.new.Ic]
      num.new.IptoIc <- length(ids.new.Ic)
      status[ids.new.Ic] <- "ic"
      statusTime[ids.new.Ic] <- at
      ids.new.Ic.other <- intersect(ids.new.Ic,which(non.office == 1))
      ids.new.Ic.office <- intersect(ids.new.Ic,which(non.office == 0))
      num.new.IptoIc.w <- length(ids.new.Ic.office)
      ### office workers
      vec.new.iso.office <- which(rbinom(length(ids.new.Ic.office), 1, iso.prob.office) == 1)
      if (length(vec.new.iso.office) > 0) {
        ids.new.iso1.office <- ids.new.Ic.office[vec.new.iso.office]
        num.new.iso1.office <- length(ids.new.iso1.office)
        isolate[ids.new.iso1.office] <- 1 # isolation for mild infection
        isoTime[ids.new.iso1.office] <- at
      }
      ### non-office workers
      vec.new.iso.other <- which(rbinom(length(ids.new.Ic.other), 1, iso.prob.other) == 1)
      if (length(vec.new.iso.other) > 0) {
        ids.new.iso1.other <- ids.new.Ic.other[vec.new.iso.other]
        num.new.iso1.other <- length(ids.new.iso1.other)
        isolate[ids.new.iso1.other] <- 1 # isolation for mild infection
        isoTime[ids.new.iso1.other] <- at
      }
    }
  }

  ## Determine Hospitalized (Ic to H) or Recovered (Ic to R) pathway
  ids.newIc <- which(active == 1 & status == "ic" & statusTime <= at & is.na(hospit))
  num.newIc <- length(ids.newIc)
  if (num.newIc > 0) {
    age.group <- pmin((floor(age[ids.newIc] / 10)) + 1, 9)
    prop.hosp.vec <- prop.hospit[age.group]
    if (any(is.na(prop.hosp.vec))) stop("error in prop.hosp.vec")

    # Vaccination reducing risk of hospitalization
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

    if (at >= hosp.boost.start & at <= hosp.boost.stop) {
      prop.hosp.vec <- prop.hosp.vec * hosp.boost.mult
    }

    if (at >= hosp.supp.start & at <= hosp.supp.stop) {
      prop.hosp.vec <- prop.hosp.vec * hosp.supp.mult
    }

    #Set pathway
    if (any(is.na(prop.hosp.vec))) stop("error in prop.hosp.vec")
    vec.new.hospit <- rbinom(num.newIc, 1, prop.hosp.vec)
    hospit[ids.newIc] <- vec.new.hospit
  }

  # Ic to H: clinical infectious move to hospitalized
  num.new.IctoH <- 0
  num.new.IctoH.w <- 0
  num.new.iso2 <- 0
  num.new.iso2.w <- 0
  ids.Ich <- which(active == 1 & status == "ic" & statusTime < at & hospit == 1)
  num.Ich <- length(ids.Ich)
  if (num.Ich > 0) {
    vec.new.H <- which(rbinom(num.Ich, 1, ich.rate) == 1)
    if (length(vec.new.H) > 0) {
      ids.new.H <- ids.Ich[vec.new.H]
      num.new.IctoH <- length(ids.new.H)
      status[ids.new.H] <- "h"
      statusTime[ids.new.H] <- at
      ids.new.H.w <- intersect(ids.new.H, which(non.office==0))
      num.new.IctoH.w <- length(ids.new.H.w)
      ids.new.iso2 <- which(status == "h" & statusTime == at & isolate == 1)
      num.new.iso2 <- length(ids.new.iso2)
      if (num.new.iso2 > 0) {
        isolate[ids.new.iso2] <- 2 # isolation for severe infection
        ids.new.iso2.w <- intersect(ids.new.iso2,which(non.office==0))
        num.new.iso2.w <- length(ids.new.iso2.w)
      }
    }
  }

  # H to R: hospitalized move to recovered
  num.new.HtoR <- 0
  num.new.HtoR.w <- 0
  ids.H <- which(active == 1 & status == "h" & statusTime < at & hospit == 1)
  num.H <- length(ids.H)
  if (num.H > 0) {
    vec.new.R <- which(rbinom(num.H, 1, hr.rate) == 1)
    if (length(vec.new.R) > 0) {
      ids.new.R <- ids.H[vec.new.R]
      num.new.HtoR <- length(ids.new.R)
      status[ids.new.R] <- "r"
      statusTime[ids.new.R] <- at
      ids.new.R.w <- intersect(ids.new.R,which(non.office == 0))
      num.new.HtoR.w <- length(ids.new.R.w)
    }
  }

  # Ic to R: clinical infectious move to recovered
  num.new.IctoR <- 0
  num.new.IctoR.w <- 0
  ids.Icr <- which(active == 1 & status == "ic" & statusTime < at & hospit == 0)
  num.Icr <- length(ids.Icr)
  if (num.Icr > 0) {
    vec.new.R <- which(rbinom(num.Icr, 1, icr.rate) == 1)
    if (length(vec.new.R) > 0) {
      ids.new.R <- ids.Icr[vec.new.R]
      num.new.IctoR <- length(ids.new.R)
      status[ids.new.R] <- "r"
      statusTime[ids.new.R] <- at
      ids.new.R.w <- intersect(ids.new.R, which(non.office == 0))
      num.new.IctoR.w <- length(ids.new.R.w)
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
  num.new.iso3 <- 0
  num.new.iso3.w <- 0
  ids.new.iso3 <- which(active == 1 & status == "r" & isolate == 1 & (at - isoTime) > 5)
  num.new.iso3 <- length(ids.new.iso3)
  if (num.new.iso3 > 0) {
    isolate[ids.new.iso3] <- 3 # post-isolation masking
    ids.new.iso3.w <- intersect(ids.new.iso3,which(non.office==0))
    num.new.iso3.w <- length(ids.new.iso3.w)
  }

  # End isolation pathway
  num.new.iso.end <- 0
  num.new.iso.end.w <- 0
  ids.new.iso.end <- which(active == 1 &
                             ((status == "r" & isolate %in% c(2,3)) | isolate == 4) &
                             (at - isoTime) > 10)
  num.new.iso.end <- length(ids.new.iso.end)
  if (num.new.iso.end > 0) {
    isolate[ids.new.iso.end] <- NA # end isolation pathway
    isoTime[ids.new.iso.end] <- NA
    ids.new.iso.end.w <- intersect(ids.new.iso.end,which(non.office==0))
    num.new.iso.end.w <- length(ids.new.iso.end.w)
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
  dat <- set_epi(dat, "ich.flow", at, num.new.IctoH)
  dat <- set_epi(dat, "hr.flow", at, num.new.HtoR)
  dat <- set_epi(dat, "rs.flow", at, num.new.RtoS)
  
  dat <- set_epi(dat, "ipic.flow.w", at, num.new.IptoIc.w)
  dat <- set_epi(dat, "icr.flow.w", at, num.new.IctoR.w)
  dat <- set_epi(dat, "ich.flow.w", at, num.new.IctoH.w)
  dat <- set_epi(dat, "hr.flow.w", at, num.new.HtoR.w)
  
  dat <- set_epi(dat, "iso1.flow", at, num.new.iso1.other + num.new.iso1.office)
  dat <- set_epi(dat, "iso2.flow", at, num.new.iso2)
  dat <- set_epi(dat, "iso3.flow", at, num.new.iso3)
  dat <- set_epi(dat, "iso4.flow", at, num.new.iso4.other + num.new.iso4.office)
  dat <- set_epi(dat, "isoend.flow", at, num.new.iso.end)
  
  dat <- set_epi(dat, "iso1.flow.w", at, num.new.iso1.office)
  dat <- set_epi(dat, "iso2.flow.w", at, num.new.iso2.w)
  dat <- set_epi(dat, "iso3.flow.w", at, num.new.iso3.w)
  dat <- set_epi(dat, "iso4.flow.w", at, num.new.iso4.office)
  dat <- set_epi(dat, "isoend.flow.w", at, num.new.iso.end.w)

  return(dat)
}
