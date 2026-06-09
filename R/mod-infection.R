#' @rdname moduleset-gmc19
#' @export
infect_general <- function(dat, at) {

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  dxStatus <- get_attr(dat, "dxStatus")
  isolate <- get_attr(dat, "isolate")

  vax <- get_attr(dat, "vax")
  last.dose.time <- get_attr(dat, "last.dose.time")
  vax.age.group <- vax_age_group_for(dat)

  ## Find infected nodes ##
  idsInf <- which(active == 1 & status %in% c("a", "ic", "ip"))

  ## Common Parameters ##
  inf.prob.a.rr <- get_param(dat, "inf.prob.a.rr")
  inf.add <- get_param(dat, "inf.add")
  inf.sub <- get_param(dat, "inf.sub")
  inf.boost.start <- get_param(dat, "inf.boost.start")
  inf.boost.stop <- get_param(dat, "inf.boost.stop")
  inf.supp.start <- get_param(dat, "inf.supp.start")
  inf.supp.stop <- get_param(dat, "inf.supp.stop")
  inf.prob.mask.rr <- get_param(dat, "inf.prob.mask.rr")
  act.rate.iso.inter.time <- get_param(dat, "act.rate.iso.inter.time")
  act.rate.iso.inter.rr <- get_param(dat, "act.rate.iso.inter.rr")
  vax.schedule <- build_vax_schedule(dat, get_pathogen(dat))

  # Seasonal forcing of transmission (issue #26): a sinusoidal annual multiplier on
  # inf.prob. Defaults to no forcing (amplitude 0) when the params are unset.
  seasonal.amp <- get_param(dat, "seasonal.amp", override.null.error = TRUE)
  seasonal.phase <- get_param(dat, "seasonal.phase", override.null.error = TRUE)
  if (is.null(seasonal.amp)) seasonal.amp <- 0
  if (is.null(seasonal.phase)) seasonal.phase <- 0
  season_mult <- 1 + seasonal.amp * cos(2 * pi * (at - seasonal.phase) / 364)

  nLayers <- dat$num.nw
  nInf <- rep(0, nLayers)
  allNewInf <- integer(0)  # newly infected across all layers + imports, for stratified incidence
 
  if (length(idsInf) > 0) {
    for (layer in seq_len(nLayers)) {
      ## Look up discordant edgelist ##
      del <- discord_edgelist(dat, at, network = layer,
                              infstat = c("a", "ic", "ip"))

      ## If any discordant pairs, proceed ##
      if (!(is.null(del))
      ) {

        ## Parameters ##
        inf.prob <- get_param(dat, "inf.prob")[layer]
        act.rate <- get_param(dat, "exposure.rate")[layer]
        inf.prob.inter.rr <- get_param(dat, "inf.prob.inter.rr")[layer]
        inf.prob.inter.time <- get_param(dat, "inf.prob.inter.time")[layer]
        act.rate.inter.rr <- get_param(dat, "act.rate.inter.rr")[layer]
        act.rate.inter.time <- get_param(dat, "act.rate.inter.time")[layer]

        # Update inf.prob to account for seasonal variation
        for (i in seq_along(inf.boost.start)) {
          if (at >= inf.boost.start[i] & at <= inf.boost.stop[i]){
            inf.prob <- inf.prob + inf.add[i]
          }
        }

        for (i in seq_along(inf.supp.start)) {
          if (at >= inf.supp.start[i] & at <= inf.supp.stop[i])
            inf.prob <- inf.prob - inf.sub[i]
        }

        # Set parameters on discordant edgelist data frame (seasonally forced)
        del$transProb <- inf.prob * season_mult

        # Vaccine effect on susceptibility/transmission
        # Dose allocation is handled in mod-vax.R; this infection module only
        # applies vaccine-derived protection when calculating per-edge transmission.
        #
        # compute_ve() uses each susceptible node's current vaccine dose,
        # dose timing, age group, and the VE parameters stored in vax.schedule
        # to return the current relative risk for infection.
        #
        # rr = 1 means no vaccine-derived protection.
        # rr < 1 reduces the susceptible node's transmission probability.
        vax_eff <-
          compute_ve(
            at = at,
            ids = del$sus,
            vax = vax,
            last.dose.time = last.dose.time,
            vax.schedule = vax.schedule,
            outcome = "infect",
            vax.age.group = vax.age.group)
        # Store vaccination status and time since latest dose
        del$vaxSus <- vax[del$sus]
        # Apply vaccine-derived susceptibility reduction to transmission probability for each discordant edge.
        del$transProb <- del$transProb * vax_eff$rr


        # Asymptomatic infection
        del$stat <- status[del$inf]
        del$transProb[del$stat == "a"] <- del$transProb[del$stat == "a"] *
          inf.prob.a.rr

        # Generic inf.prob and act.rate interventions
        if (at >= inf.prob.inter.time) {
          del$transProb <- del$transProb * inf.prob.inter.rr
        }
        del$actRate <- act.rate
        if (at >= act.rate.inter.time) {
          del$actRate <- del$actRate * act.rate.inter.rr
        }

        # Case isolation for those in isolation process
        del$iso <- ifelse(!is.na(isolate[del$inf]),isolate[del$inf],0)
        if (at >= act.rate.iso.inter.time) {
          del$actRate[del$iso %in% c(1,2)] <- del$actRate[del$iso %in% c(1,2)] *
            act.rate.iso.inter.rr
        }

        # Masking for those in isolation process
        if (at >= act.rate.iso.inter.time) {
          del$transProb[del$iso %in% c(1,2,3,4)] <- del$transProb[del$iso %in% c(1,2,3,4)] *
            inf.prob.mask.rr
        }

        del$finalProb <- 1 - (1 - del$transProb)^del$actRate
        
        # debug warning message: In rbinom(nrow(del), 1, del$finalProb) : NAs produced
        bad <- which(is.na(del$finalProb) | del$finalProb < 0 | del$finalProb > 1)
        if (length(bad) > 0) {
          print(table(is.na(del$transProb), useNA = "ifany"))
          print(table(is.na(del$actRate), useNA = "ifany"))
          print(summary(del$transProb))
          print(summary(del$actRate))
          print(summary(del$finalProb))
          
          #browser()
        }
        
        # Stochastic transmission process
        transmit <- rbinom(nrow(del), 1, del$finalProb)

        # Keep rows where transmission occurred
        del <- del[which(transmit == 1), , drop = FALSE]

        # Look up new ids if any transmissions occurred
        idsNewInf <- unique(del$sus)
        nInf[layer] <- length(idsNewInf)

        # Set new attributes for those newly infected
        if (nInf[layer] > 0) {
          dat <- set_attr(dat, "status", "e", idsNewInf)
          dat <- set_attr(dat, "infTime", at, idsNewInf)
          dat <- set_attr(dat, "statusTime", at, idsNewInf)
          allNewInf <- c(allNewInf, idsNewInf)
        }
      }
    }
  }

  ## External re-importation (seasonal): a small per-step hazard of infection from
  ## outside the modeled population, so seasonal epidemics can re-ignite and zero
  ## prevalence is not absorbing. Defaults to 0 (off) when import.rate is unset.
  import.rate <- get_param(dat, "import.rate", override.null.error = TRUE)
  nImport <- 0L
  if (!is.null(import.rate) && import.rate > 0) {
    status_now <- get_attr(dat, "status")
    idsSus <- which(active == 1 & status_now == "s")
    if (length(idsSus) > 0) {
      import_prob <- min(1, import.rate * season_mult)
      idsImp <- idsSus[runif(length(idsSus)) < import_prob]
      if (length(idsImp) > 0) {
        dat <- set_attr(dat, "status", "e", idsImp)
        dat <- set_attr(dat, "infTime", at, idsImp)
        dat <- set_attr(dat, "statusTime", at, idsImp)
        nImport <- length(idsImp)
        allNewInf <- c(allNewInf, idsImp)
      }
    }
  }
  allNewInf <- unique(allNewInf)

  ## Summary statistics for incidence (network by layer + imports)
  dat <- set_epi(dat, "se.flow", at, sum(nInf) + nImport)
  dat <- set_epi(dat, "se.flow.l1", at, nInf[1])
  dat <- set_epi(dat, "se.flow.l2", at, nInf[2])
  dat <- set_epi(dat, "se.flow.l3", at, nInf[3])
  dat <- set_epi(dat, "se.flow.l4", at, if (length(nInf) >= 4) nInf[4] else 0)  # household layer
  dat <- set_epi(dat, "se.import.flow", at, nImport)

  ## Stratified incidence: by age group, vaccination status, and network position
  age <- get_attr(dat, "age")
  age.breaks <- get_param(dat, "age.breaks")
  ag <- cut(age[allNewInf], age.breaks, labels = FALSE, right = FALSE)
  for (g in seq_len(length(age.breaks) - 1)) {
    dat <- set_epi(dat, paste0("se.flow.age", g), at, sum(ag == g, na.rm = TRUE))
  }
  vaxed_new <- vax[allNewInf] >= 1
  dat <- set_epi(dat, "se.flow.vax", at, sum(vaxed_new, na.rm = TRUE))
  dat <- set_epi(dat, "se.flow.unvax", at, sum(!vaxed_new, na.rm = TRUE))
  degree_quartile <- get_attr(dat, "degree_quartile")
  if (!is.null(degree_quartile)) {
    dq_new <- degree_quartile[allNewInf]
    for (q in 1:4) dat <- set_epi(dat, paste0("se.flow.degQ", q), at, sum(dq_new == q, na.rm = TRUE))
  }
  is_bridge <- get_attr(dat, "is_bridge")
  if (!is.null(is_bridge)) {
    br_new <- is_bridge[allNewInf]
    dat <- set_epi(dat, "se.flow.bridge", at, sum(br_new, na.rm = TRUE))
    dat <- set_epi(dat, "se.flow.nonbridge", at, sum(!br_new, na.rm = TRUE))
  }

  return(dat)
}


infect_covid_corporate <- function(dat, at) { #maria's repo
  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  dxStatus <- get_attr(dat, "dxStatus")
  isolate <- get_attr(dat, "isolate")
  
  vax <- get_attr(dat, "vax")
  vax1Time <- get_attr(dat, "vax1Time")
  vax2Time <- get_attr(dat, "vax2Time")
  vax3Time <- get_attr(dat, "vax3Time")
  vax4Time <- get_attr(dat, "vax4Time")
  
  ## Find infected nodes ##
  idsInf <- which(active == 1 & status %in% c("a", "ic", "ip"))
  
  ## Common Parameters ##
  inf.prob.a.rr <- get_param(dat, "inf.prob.a.rr")
  vax1.rr.infect <- get_param(dat, "vax1.rr.infect")
  vax2.rr.infect <- get_param(dat, "vax2.rr.infect")
  vax3.rr.infect <- get_param(dat, "vax3.rr.infect")
  vax4.rr.infect <- get_param(dat, "vax4.rr.infect")
  half.life <- get_param(dat, "half.life")
  inf.add <- get_param(dat, "inf.add")
  inf.sub <- get_param(dat, "inf.sub")
  inf.boost.start <- get_param(dat, "inf.boost.start")
  inf.boost.stop <- get_param(dat, "inf.boost.stop")
  inf.supp.start <- get_param(dat, "inf.supp.start")
  inf.supp.stop <- get_param(dat, "inf.supp.stop")
  inf.prob.mask.rr <- get_param(dat, "inf.prob.mask.rr")
  act.rate.iso.inter.time <- get_param(dat, "act.rate.iso.inter.time")
  act.rate.iso.inter.rr <- get_param(dat, "act.rate.iso.inter.rr")
  
  nLayers <- dat$num.nw
  nInf <- rep(0, nLayers)
  
  if (length(idsInf) > 0) {
    for (layer in seq_len(nLayers)) {
      ## Look up discordant edgelist ##
      del <- discord_edgelist(dat, at, network = layer,
                              infstat = c("a", "ic", "ip"))
      
      ## If any discordant pairs, proceed ##
      if (!(is.null(del))
          ) {
        
        ## Parameters ##
        inf.prob <- get_param(dat, "inf.prob")[layer]
        act.rate <- get_param(dat, "act.rate")[layer]
        inf.prob.inter.rr <- get_param(dat, "inf.prob.inter.rr")[layer]
        inf.prob.inter.time <- get_param(dat, "inf.prob.inter.time")[layer]
        act.rate.inter.rr <- get_param(dat, "act.rate.inter.rr")[layer]
        act.rate.inter.time <- get_param(dat, "act.rate.inter.time")[layer]
        
        # Update inf.prob to account for seasonal variation
        for (i in seq_along(inf.boost.start)) {
          if (at >= inf.boost.start[i] & at <= inf.boost.stop[i]){
            inf.prob <- inf.prob + inf.add[i]
          }
        }
        
        for (i in seq_along(inf.supp.start)) {
          if (at >= inf.supp.start[i] & at <= inf.supp.stop[i])
            inf.prob <- inf.prob - inf.sub[i]
        }
        
        # Set parameters on discordant edgelist data frame
        del$transProb <- inf.prob
        
        # Vaccination
        del$vaxSus <- vax[del$sus]
        del$transProb[del$vaxSus == 1] <- del$transProb[del$vaxSus == 1] * vax1.rr.infect
        del$transProb[del$vaxSus == 2] <- del$transProb[del$vaxSus == 2] * vax2.rr.infect
        del$transProb[del$vaxSus == 3] <- del$transProb[del$vaxSus == 3] * vax3.rr.infect
        del$transProb[del$vaxSus == 4] <- del$transProb[del$vaxSus == 4] * vax4.rr.infect
        
        #Waning vaccine immunity
        sinceVax1 <- at - vax1Time[del$sus] # days since 1st vax
        sinceVax2 <- at - vax2Time[del$sus] # days since 2nd vax
        sinceVax3 <- at - vax3Time[del$sus] # days since 3rd vax
        sinceVax4 <- at - vax4Time[del$sus] # days since 4th vax
        
        latest.vax <- pmin(sinceVax1, sinceVax2, sinceVax3, sinceVax4, na.rm = TRUE) # identify the most recent vax by taking the minimum time since vaccination across doses.
        latest.vax[is.na(latest.vax)] <- 0 #  individuals with no vax are assigned 0 - no waning 
        
        del$latest.vax <- latest.vax
        del$transProb <- pmin(del$transProb * (2 ^ (del$latest.vax / half.life)), inf.prob) # Vaccine protection as exponential waning: the transmission probability increases by a factor of 2 every half.life time units since the most recent dose, capped at the baseline transmission probability (inf.prob)
        
        # Asymptomatic infection
        del$stat <- status[del$inf]
        del$transProb[del$stat == "a"] <- del$transProb[del$stat == "a"] *
          inf.prob.a.rr
        
        # Generic inf.prob and act.rate interventions
        if (at >= inf.prob.inter.time) {
          del$transProb <- del$transProb * inf.prob.inter.rr
        }
        del$actRate <- act.rate
        if (at >= act.rate.inter.time) {
          del$actRate <- del$actRate * act.rate.inter.rr
        }
        
        # Case isolation for those in isolation process
        del$iso <- ifelse(!is.na(isolate[del$inf]),isolate[del$inf],0)
        if (at >= act.rate.iso.inter.time) {
          del$actRate[del$iso %in% c(1,2)] <- del$actRate[del$iso %in% c(1,2)] *
            act.rate.iso.inter.rr
        }
        
        # Masking for those in isolation process
        if (at >= act.rate.iso.inter.time) {
          del$transProb[del$iso %in% c(1,2,3,4)] <- del$transProb[del$iso %in% c(1,2,3,4)] *
            inf.prob.mask.rr
        }
        
        del$finalProb <- 1 - (1 - del$transProb)^del$actRate
        
        # Stochastic transmission process
        transmit <- rbinom(nrow(del), 1, del$finalProb)
        
        # Keep rows where transmission occurred
        del <- del[which(transmit == 1), , drop = FALSE]
        
        # Look up new ids if any transmissions occurred
        idsNewInf <- unique(del$sus)
        nInf[layer] <- length(idsNewInf)
        
        # Set new attributes for those newly infected
        if (nInf[layer] > 0) {
          dat <- set_attr(dat, "status", "e", idsNewInf)
          dat <- set_attr(dat, "infTime", at, idsNewInf)
          dat <- set_attr(dat, "statusTime", at, idsNewInf)
        }
      }
    }
  }
  
  ## Summary statistics for incidence
  dat$epi$se.flow[at] <- sum(nInf)
  dat$epi$se.flow.l1[at] <- nInf[1]
  dat$epi$se.flow.l2[at] <- nInf[2]
  dat$epi$se.flow.l3[at] <- nInf[3]
  
  return(dat)
}

