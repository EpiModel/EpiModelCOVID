
#' @rdname moduleset-corporate
#' @export
infect_covid_corporate <- function(dat, at) {

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
  inf.prob.mask.rr <- get_param(dat, "inf.prob.mask.rr")
  act.rate.iso.inter.time <- get_param(dat, "act.rate.iso.inter.time")
  act.rate.iso.inter.rr <- get_param(dat, "act.rate.iso.inter.rr")

  nLayers <- length(dat$el)
  nInf <- rep(0, nLayers)

  if (length(idsInf) > 0) {
    for (layer in seq_len(nLayers)) {
      ## Look up discordant edgelist ##
      del <- discord_edgelist(dat, at, network = layer,
                              infstat = c("a", "ic", "ip"))

      ## If any discordant pairs, proceed ##
      if (!(is.null(del))) {

        ## Parameters ##
        inf.prob <- get_param(dat, "inf.prob")[layer]
        act.rate <- get_param(dat, "act.rate")[layer]
        inf.prob.inter.rr <- get_param(dat, "inf.prob.inter.rr")[layer]
        inf.prob.inter.time <- get_param(dat, "inf.prob.inter.time")[layer]
        act.rate.inter.rr <- get_param(dat, "act.rate.inter.rr")[layer]
        act.rate.inter.time <- get_param(dat, "act.rate.inter.time")[layer]

        # Set parameters on discordant edgelist data frame
        del$transProb <- inf.prob

        # Vaccination
        del$vaxSus <- vax[del$sus]
        del$transProb[del$vaxSus == 1] <- del$transProb[del$vaxSus == 1] * vax1.rr.infect
        del$transProb[del$vaxSus == 2] <- del$transProb[del$vaxSus == 2] * vax2.rr.infect
        del$transProb[del$vaxSus == 3] <- del$transProb[del$vaxSus == 3] * vax3.rr.infect
        del$transProb[del$vaxSus == 4] <- del$transProb[del$vaxSus == 4] * vax4.rr.infect

        #Waning vaccine immunity
        sinceVax1 <- at - vax1Time[del$sus]
        sinceVax2 <- at - vax2Time[del$sus]
        sinceVax3 <- at - vax3Time[del$sus]
        sinceVax4 <- at - vax4Time[del$sus]

        latest.vax <- pmin(sinceVax1, sinceVax2, sinceVax3, sinceVax4, na.rm = TRUE)
        latest.vax[is.na(latest.vax)] <- 0

        del$latest.vax <- latest.vax
        del$transProb <- pmin(del$transProb * (2 ^ (del$latest.vax / half.life)), inf.prob)

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
