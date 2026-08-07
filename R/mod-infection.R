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
  # Fetched here rather than after the transmission loop because the issue #45
  # age-specific susceptibility knob needs it per discordant edge.
  age <- get_attr(dat, "age")

  # Direct infant product (issue #41): when TRUE (the RSV infant_direct /
  # infant_direct_cocoon arms) a nirsevimab/maternal-like product is severity-
  # dominant, so dosed infants receive NO infection-VE (little effect on
  # acquisition) while their severity VE still applies downstream in mod-progress
  # and mod-departure. Default FALSE leaves every other arm unchanged.
  is_infant <- get_attr(dat, "is_infant")
  vax.infant.product <- get_param(dat, "vax.infant.product", override.null.error = TRUE)
  if (is.null(vax.infant.product)) vax.infant.product <- FALSE

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
        # Direct infant product (issue #41): strip infection-VE for dosed infants,
        # so a severity-dominant monoclonal/maternal product does not spuriously
        # cut infant acquisition; adult recipients (e.g. the cocooning co-residents)
        # keep their normal infection VE.
        if (isTRUE(vax.infant.product) && !is.null(is_infant)) {
          inf_sus <- which(is_infant[del$sus])
          if (length(inf_sus) > 0) vax_eff$rr[inf_sus] <- 1
        }
        # Apply vaccine-derived susceptibility reduction to transmission probability for each discordant edge.
        del$transProb <- del$transProb * vax_eff$rr

        # Age-specific susceptibility (issue #45 trajectory knob). A relative
        # risk on the SUSCEPTIBLE node's age group, applied multiplicatively to
        # the per-contact transmission probability. Absent, all-NA or all-ones
        # means no effect, which is the current behaviour, so this stays inert
        # unless a scenario sets it.
        #
        # It exists because the model's influenza age gradient is emergent from
        # network structure alone: rural flu incidence peaks in 10-19y because
        # that is where the school layer concentrates contact, whereas Krishnan
        # 2018, Sullender 2019 and PHIRST all give the opposite ordering. There
        # was no other lever on that gradient, since inf.prob is per-layer and
        # carries no age dimension. Sweeping susc.age.rr from all-ones to a
        # monotone decline is the "corrected" arm of that axis.
        susc.age.rr <- get_param(dat, "susc.age.rr", override.null.error = TRUE)
        if (!is.null(susc.age.rr) && !all(is.na(susc.age.rr)) &&
            any(susc.age.rr != 1, na.rm = TRUE)) {
          ag_sus <- cut(age[del$sus], get_param(dat, "age.breaks"),
                        labels = FALSE, right = FALSE)
          rr_sus <- susc.age.rr[ag_sus]
          rr_sus[is.na(rr_sus)] <- 1
          del$transProb <- del$transProb * rr_sus
        }

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
  # Age x vaccination-status cross-tab. The paper's indirect-protection claim
  # (vaccinating network hubs cuts infection among the UNVACCINATED vulnerable)
  # cannot be read from the age and vax marginals alone, so emit incidence among
  # the unvaccinated (and vaccinated) within each age group.
  for (g in seq_len(length(age.breaks) - 1)) {
    dat <- set_epi(dat, paste0("se.flow.unvax.age", g), at, sum(ag == g & !vaxed_new, na.rm = TRUE))
    dat <- set_epi(dat, paste0("se.flow.vax.age", g),   at, sum(ag == g &  vaxed_new, na.rm = TRUE))
  }
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
  # Infant incidence (issue #23): reported separately (and by vax status) so the
  # infant burden and the indirect protection of infants by parent-targeting are
  # not diluted in the coarse youngest age band.
  is_infant <- get_attr(dat, "is_infant")
  if (!is.null(is_infant)) {
    inf_new <- is_infant[allNewInf]
    dat <- set_epi(dat, "se.flow.infant", at, sum(inf_new, na.rm = TRUE))
    dat <- set_epi(dat, "se.flow.infant.unvax", at, sum(inf_new & !vaxed_new, na.rm = TRUE))
    dat <- set_epi(dat, "se.flow.infant.vax", at, sum(inf_new & vaxed_new, na.rm = TRUE))
    # Infant severity is sharply front-loaded, so the RSV infant contrast has to
    # be readable at the 0-6 / 6-12 month split rather than pooled under 1y
    # (issue #46). Deaths and hospitalisations already split this way; incidence
    # did not, which left the infant comparison undecomposable into its
    # attack-rate and severity parts.
    a_new <- age[allNewInf]
    dat <- set_epi(dat, "se.flow.infant.young", at, sum(inf_new & a_new <  0.5, na.rm = TRUE))
    dat <- set_epi(dat, "se.flow.infant.old",   at, sum(inf_new & a_new >= 0.5, na.rm = TRUE))
  }

  ## Person-level cumulative infection (issue #46) ------------------------------
  ## `se.flow` counts S->E transitions, i.e. EPISODES. With reinfection on, a
  ## cumulative per-100 built from it can and does exceed 100, which is not an
  ## attack rate and cannot be set against any seroprevalence or cohort estimate.
  ## `ever.inf` makes the person the unit: set once, never cleared, so summing it
  ## over the living population is the number of people ever infected.
  ever.inf <- get_attr(dat, "ever.inf")
  if (!is.null(ever.inf) && length(allNewInf) > 0) {
    ever.inf[allNewInf] <- 1
    dat <- set_attr(dat, "ever.inf", ever.inf)
  }

  ## Household secondary attack rate (issue #46) --------------------------------
  ## The one external anchor that survived the validation review is the PHIRST
  ## cross-pathogen household ordering, which is a per-person risk among exposed
  ## household contacts. `se.flow.l4` is the numerator already (layer 4 is the
  ## household layer), but with no denominator it is a count, not a rate.
  ##
  ## Denominator: susceptible nodes sharing a household with at least one
  ## infectious node at this step. Summing numerator and denominator over the
  ## epidemic and dividing gives a household SAR directly comparable to a cohort
  ## study's. Emitted as two series rather than a ratio so the analysis can
  ## aggregate before dividing; a per-step ratio would be undefined whenever no
  ## household is exposed.
  hh.ids <- get_attr(dat, "hh.ids")
  if (!is.null(hh.ids)) {
    status_now <- get_attr(dat, "status")
    hh_infectious <- unique(hh.ids[active == 1 & status_now %in% c("a", "ip", "ic", "h")])
    exposed_sus <- active == 1 & status_now == "s" & hh.ids %in% hh_infectious
    dat <- set_epi(dat, "hh.sar.den", at, sum(exposed_sus, na.rm = TRUE))
    dat <- set_epi(dat, "hh.sar.num", at, if (length(nInf) >= 4) nInf[4] else 0)
    dat <- set_epi(dat, "hh.exposed.n", at, length(hh_infectious))
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

