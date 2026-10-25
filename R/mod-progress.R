#' @rdname moduleset-common
#' @export
progress_general <- function(dat, at) {

  ## Attributes
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  statusTime <- get_attr(dat, "statusTime")
  clinical <- get_attr(dat, "clinical")
  hospit <- get_attr(dat, "hospit")
  fatal <- get_attr(dat, "fatal")
  recTime <- get_attr(dat, "recTime")
  inf.count <- get_attr(dat, "inf.count")
  age <- get_attr(dat, "age")
  vax.age.group <- vax_age_group_for(dat)
  vax <- get_attr(dat, "vax")
  last.dose.time <- get_attr(dat, "last.dose.time")
  deg_work <- get_attr(dat, "deg_work")
  non.office <- ifelse(deg_work > 0, 0, 1)  # 0=office, 1=non-office
  is_infant <- get_attr(dat, "is_infant")

  ## Parameters
  # Severity per infection by severity band (GlobalMix India #71; see
  # mod-immunity.R): the symptomatic share, and the admission and death
  # probabilities per infection, which a symptomatic case takes divided by the
  # symptomatic share
  prop.clinical <- get_param(dat, "prop.clinical")
  ihr <- get_param(dat, "ihr")
  ifr <- get_param(dat, "ifr")
  death.hosp.rr <- param_or(dat, "death.hosp.rr", 1)
  ea.rate <- get_param(dat, "ea.rate")
  ar.rate <- get_param(dat, "ar.rate")
  eip.rate <- get_param(dat, "eip.rate")
  ipic.rate <- get_param(dat, "ipic.rate")
  ich.rate <- get_param(dat, "ich.rate")
  icr.rate <- get_param(dat, "icr.rate")
  hr.rate <- get_param(dat, "hr.rate")
  # Days of full protection after recovery (#69); Inf, no return to s
  immune.full <- param_or(dat, "immune.full", Inf)
  age.breaks <- get_param(dat, "age.breaks")
  vax.schedule <- build_vax_schedule(dat, get_pathogen(dat)) # get vaccination schedule
  iprog <- infant_prog_state(dat)  # the infant product (#76)
  sevdiag <- sev_diagnostics(dat)
  vaxdiag <- vax_diagnostics(dat)
  nonfatal <- !(fatal %in% 1)

  ## Determine Subclinical (E to A) or Clinical (E to Ip to Ic) pathway, once per
  ## episode (#69): start_episode() clears clinical at every infection
  ids.newInf <- which(active == 1 & status == "e" & statusTime <= at & is.na(clinical))
  num.newInf <- length(ids.newInf)
  if (num.newInf > 0) {
    prop.clin.vec <- prop.clinical[sev_band(dat, age[ids.newInf])]

    # Vaccine protection against symptomatic disease (compute_ve(), mod-vax.R;
    # dose allocation is in mod-vax.R) and protection from an earlier infection
    # in the run (#69), each a relative risk given infection
    vax_eff_clinical <- compute_ve(
      at = at,
      ids = ids.newInf,
      vax = vax,
      last.dose.time = last.dose.time,
      vax.schedule = vax.schedule,
      outcome = "clinical",
      vax.age.group = vax.age.group,
      iprog = iprog
    )
    prop.clin.vec <- prop.clin.vec * vax_eff_clinical$rr * immune_rr(dat, ids.newInf, "clinical")

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
      recTime[ids.new.R] <- at
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
  num.new.IptoIc.first <- 0

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
      ## Issue #46. Every surveillance number we could be compared against is a
      ## count of symptomatic care-seeking PEOPLE, while `prop.clinical` was an
      ## input that never surfaced as an output. Flag the first symptomatic
      ## episode per person and emit onset by age band, so the clinical fraction
      ## becomes observable and the age distribution of symptomatic disease is
      ## readable rather than inferred from total incidence. ipic.flow.first
      ## counts a person's first symptomatic episode (#69, #77).
      ever.sympt <- get_attr(dat, "ever.sympt")
      if (!is.null(ever.sympt)) {
        num.new.IptoIc.first <- sum(ever.sympt[ids.new.Ic] == 0, na.rm = TRUE)
        ever.sympt[ids.new.Ic] <- 1
        dat <- set_attr(dat, "ever.sympt", ever.sympt)
      }
      ag_ic <- cut(age[ids.new.Ic], age.breaks, labels = FALSE, right = FALSE)
      for (g in seq_len(length(age.breaks) - 1)) {
        dat <- set_epi(dat, paste0("ipic.flow.age", g), sum(ag_ic == g, na.rm = TRUE), at = at)
      }
      if (!is.null(is_infant)) {
        i_ic <- is_infant[ids.new.Ic]; a_ic <- age[ids.new.Ic]
        dat <- set_epi(dat, "ipic.flow.infant", sum(i_ic, na.rm = TRUE), at = at)
        dat <- set_epi(dat, "ipic.flow.infant.young", sum(i_ic & a_ic <  0.5, na.rm = TRUE), at = at)
        dat <- set_epi(dat, "ipic.flow.infant.old",   sum(i_ic & a_ic >= 0.5, na.rm = TRUE), at = at)
      }
      if (sevdiag) dat <- sev_flows(dat, at, "sympt", ids.new.Ic, inf.count[ids.new.Ic] == 1)
      if (vaxdiag) dat <- vd_flows(dat, at, "sympt", ids.new.Ic)
    }
  }

  ## Admission and death, drawn once per episode at onset (#71). A symptomatic
  ## case is admitted with probability ihr / prop.clinical and dies with
  ## probability ifr / prop.clinical, each scaled by the vaccine's and by prior
  ## infection's protection given symptomatic disease. A fatal case dies when it
  ## would otherwise have recovered, at the end of its symptomatic stay or of
  ## its hospital stay (deaths_covid_gmc19(), mod-departure.R). The daily
  ## disease-death hazard in ic and h this replaces made the IFR depend on how
  ## long cases stayed symptomatic.
  ids.newIc <- which(active == 1 & status == "ic" & statusTime <= at & is.na(hospit))
  num.newIc <- length(ids.newIc)
  if (num.newIc > 0) {
    b <- sev_band(dat, age[ids.newIc])
    pc <- prop.clinical[b]
    prop.hosp.vec <- ifelse(pc > 0, ihr[b] / pc, 0)
    prop.death.vec <- ifelse(pc > 0, ifr[b] / pc, 0)
    if (any(is.na(prop.hosp.vec)) || any(is.na(prop.death.vec))) stop("error in prop.hosp.vec")

    # Vaccine protection against admission and death given symptomatic disease
    # (compute_ve(), mod-vax.R), and protection from an earlier infection (#69)
    vax_eff_hosp <- compute_ve(
      at = at,
      ids = ids.newIc,
      vax = vax,
      last.dose.time = last.dose.time,
      vax.schedule = vax.schedule,
      outcome = "hosp",
      vax.age.group = vax.age.group,
      iprog = iprog
    )
    vax_eff_death <- compute_ve(
      at = at, ids = ids.newIc, vax = vax, last.dose.time = last.dose.time,
      vax.schedule = vax.schedule, outcome = "death", vax.age.group = vax.age.group, iprog = iprog)
    prop.hosp.vec <- pmin(1, prop.hosp.vec * vax_eff_hosp$rr * immune_rr(dat, ids.newIc, "hosp"))
    prop.death.vec <- prop.death.vec * vax_eff_death$rr * immune_rr(dat, ids.newIc, "death")

    vec.new.hospit <- rbinom(num.newIc, 1, prop.hosp.vec)
    hospit[ids.newIc] <- vec.new.hospit
    fatal[ids.newIc] <- rbinom(num.newIc, 1,
                               death_given_admission(prop.death.vec, prop.hosp.vec, vec.new.hospit,
                                                     death.hosp.rr))
    nonfatal <- !(fatal %in% 1)
  }

  # Ic to H: clinical infectious move to hospitalized
  num.new.IctoH <- 0
  num.new.IctoH.w <- 0

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
      ## Issue #46: hospitalisation by age band and infant sub-band. `ich.flow`
      ## was a single pooled count, so `cum_hosp` could not be expressed as the
      ## per-100,000 per-season age-specific rate that every hospital-based
      ## surveillance series reports. The denominators come from the n.age*
      ## series emitted in mod-prevalence.R.
      ag_h <- cut(age[ids.new.H], age.breaks, labels = FALSE, right = FALSE)
      for (g in seq_len(length(age.breaks) - 1)) {
        dat <- set_epi(dat, paste0("ich.flow.age", g), sum(ag_h == g, na.rm = TRUE), at = at)
      }
      if (!is.null(is_infant)) {
        i_h <- is_infant[ids.new.H]; a_h <- age[ids.new.H]
        dat <- set_epi(dat, "ich.flow.infant", sum(i_h, na.rm = TRUE), at = at)
        dat <- set_epi(dat, "ich.flow.infant.young", sum(i_h & a_h <  0.5, na.rm = TRUE), at = at)
        dat <- set_epi(dat, "ich.flow.infant.old",   sum(i_h & a_h >= 0.5, na.rm = TRUE), at = at)
      }
      if (sevdiag) dat <- sev_flows(dat, at, "hosp", ids.new.H, inf.count[ids.new.H] == 1)
      if (vaxdiag) dat <- vd_flows(dat, at, "hosp", ids.new.H)
    }
  }

  # H to R: hospitalized move to recovered; a fatal case dies instead
  # (deaths_covid_gmc19())
  num.new.HtoR <- 0
  num.new.HtoR.w <- 0
  ids.H <- which(active == 1 & status == "h" & statusTime < at & hospit == 1 & nonfatal)
  num.H <- length(ids.H)
  if (num.H > 0) {
    vec.new.R <- which(rbinom(num.H, 1, hr.rate) == 1)
    if (length(vec.new.R) > 0) {
      ids.new.R <- ids.H[vec.new.R]
      num.new.HtoR <- length(ids.new.R)
      status[ids.new.R] <- "r"
      statusTime[ids.new.R] <- at
      recTime[ids.new.R] <- at
      ids.new.R.w <- intersect(ids.new.R,which(non.office == 0))
      num.new.HtoR.w <- length(ids.new.R.w)
    }
  }

  # Ic to R: clinical infectious move to recovered; a fatal case dies instead
  num.new.IctoR <- 0
  num.new.IctoR.w <- 0
  ids.Icr <- which(active == 1 & status == "ic" & statusTime < at & hospit == 0 & nonfatal)
  num.Icr <- length(ids.Icr)
  if (num.Icr > 0) {
    vec.new.R <- which(rbinom(num.Icr, 1, icr.rate) == 1)
    if (length(vec.new.R) > 0) {
      ids.new.R <- ids.Icr[vec.new.R]
      num.new.IctoR <- length(ids.new.R)
      status[ids.new.R] <- "r"
      statusTime[ids.new.R] <- at
      recTime[ids.new.R] <- at
      ids.new.R.w <- intersect(ids.new.R, which(non.office == 0))
      num.new.IctoR.w <- length(ids.new.R.w)
    }
  }

  # R to S (#69): back to susceptible immune.full days after recovery, with the
  # partial protection of mod-immunity.R. It replaces a daily return
  # probability (rs.rate), under which a quarter to almost half of the
  # recovered were fully susceptible again within 90 days. The prior-immune
  # (#70) have no recovery in the run and stay in r.
  num.new.RtoS <- 0
  ids.new.RS <- which(active == 1 & status == "r" & !is.na(recTime) & at - recTime >= immune.full)
  if (length(ids.new.RS) > 0) {
    num.new.RtoS <- length(ids.new.RS)
    status[ids.new.RS] <- "s"
    statusTime[ids.new.RS] <- at
  }

  ## Save updated attributes
  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "statusTime", statusTime)
  dat <- set_attr(dat, "clinical", clinical)
  dat <- set_attr(dat, "hospit", hospit)
  dat <- set_attr(dat, "fatal", fatal)
  dat <- set_attr(dat, "recTime", recTime)

  ## Save summary statistics
  dat <- set_epi(dat, "ea.flow", num.new.EtoA, at = at)
  dat <- set_epi(dat, "ar.flow", num.new.AtoR, at = at)
  dat <- set_epi(dat, "eip.flow", num.new.EtoIp, at = at)
  dat <- set_epi(dat, "ipic.flow", num.new.IptoIc, at = at)
  dat <- set_epi(dat, "ipic.flow.first", num.new.IptoIc.first, at = at)
  dat <- set_epi(dat, "icr.flow", num.new.IctoR, at = at)
  dat <- set_epi(dat, "ich.flow", num.new.IctoH, at = at)
  dat <- set_epi(dat, "hr.flow", num.new.HtoR, at = at)
  dat <- set_epi(dat, "rs.flow", num.new.RtoS, at = at)

  dat <- set_epi(dat, "ipic.flow.w", num.new.IptoIc.w, at = at)
  dat <- set_epi(dat, "icr.flow.w", num.new.IctoR.w, at = at)
  dat <- set_epi(dat, "ich.flow.w", num.new.IctoH.w, at = at)
  dat <- set_epi(dat, "hr.flow.w", num.new.HtoR.w, at = at)


  return(dat)
}


progress_covid_corporate <- function(dat, at) {
  #if (at>=10) browser()
  
  ## Attributes
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  statusTime <- get_attr(dat, "statusTime")
  clinical <- get_attr(dat, "clinical")
  hospit <- get_attr(dat, "hospit")
  age <- get_attr(dat, "age")
  vax <- get_attr(dat, "vax")
  vax1Time <- get_attr(dat, "vax1Time")
  vax2Time <- get_attr(dat, "vax2Time")
  vax3Time <- get_attr(dat, "vax3Time")
  vax4Time <- get_attr(dat, "vax4Time")
  dxStatus <- get_attr(dat, "dxStatus")
  deg_work <- get_attr(dat, "deg_work")
  non.office <- ifelse(deg_work > 0, 0, 1)  # 0=office, 1=non-office
  
  
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
  age.breaks<- get_param(dat, "age.breaks")
  
  ## Determine Subclinical (E to A) or Clinical (E to Ip to Ic) pathway
  ids.newInf <- which(active == 1 & status == "e" & statusTime <= at & is.na(clinical))
  num.newInf <- length(ids.newInf)
  if (num.newInf > 0) {
    age.group <- 
      cut(
        age[ids.newInf],
        breaks = age.breaks,
        right = F,
        labels = 1:6
       
      ) |> as.character() |> as.integer()
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
    }
  }
  
  ## Determine Hospitalized (Ic to H) or Recovered (Ic to R) pathway
  ids.newIc <- which(active == 1 & status == "ic" & statusTime <= at & is.na(hospit))
  num.newIc <- length(ids.newIc)
  if (num.newIc > 0) {
    age.group <- 
      cut(
        age[ids.newIc],
        breaks = age.breaks,
        right = F,
        labels = 1:6
        
      ) |> as.character() |> as.integer()
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
  
  
  
  ## Save updated attributes
  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "statusTime", statusTime)
  dat <- set_attr(dat, "clinical", clinical)
  dat <- set_attr(dat, "hospit", hospit)
  dat <- set_attr(dat, "dxStatus", dxStatus)
  
  ## Save summary statistics
  dat <- set_epi(dat, "ea.flow", num.new.EtoA, at = at)
  dat <- set_epi(dat, "ar.flow", num.new.AtoR, at = at)
  dat <- set_epi(dat, "eip.flow", num.new.EtoIp, at = at)
  dat <- set_epi(dat, "ipic.flow", num.new.IptoIc, at = at)
  dat <- set_epi(dat, "icr.flow", num.new.IctoR, at = at)
  dat <- set_epi(dat, "ich.flow", num.new.IctoH, at = at)
  dat <- set_epi(dat, "hr.flow", num.new.HtoR, at = at)
  dat <- set_epi(dat, "rs.flow", num.new.RtoS, at = at)
  
  dat <- set_epi(dat, "ipic.flow.w", num.new.IptoIc.w, at = at)
  dat <- set_epi(dat, "icr.flow.w", num.new.IctoR.w, at = at)
  dat <- set_epi(dat, "ich.flow.w", num.new.IctoH.w, at = at)
  dat <- set_epi(dat, "hr.flow.w", num.new.HtoR.w, at = at)
  
  
  return(dat)
}


