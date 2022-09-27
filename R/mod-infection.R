
#' @rdname moduleset-vaxDecisions
#' @export
infect_covid_vax_decisions <- function(dat, at) {

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  dxStatus <- get_attr(dat, "dxStatus")
  vax <- get_attr(dat, "vax")
  
  vax1Time <- get_attr(dat, "vax1Time")
  vax2Time <- get_attr(dat, "vax2Time")
  vax3Time <- get_attr(dat, "vax3Time")
  vax4Time <- get_attr(dat, "vax4Time")

  ## Find infected nodes ##
  idsInf <- which(active == 1 & status %in% c("a", "ic", "ip"))

  ## Common Parameters ##
  inf.prob.a.rr <- get_param(dat, "inf.prob.a.rr")
  act.rate.dx.inter.rr <- get_param(dat, "act.rate.dx.inter.rr")
  act.rate.dx.inter.time <- get_param(dat, "act.rate.dx.inter.time")
  act.rate.sympt.inter.rr <- get_param(dat, "act.rate.sympt.inter.rr")
  act.rate.sympt.inter.time <- get_param(dat, "act.rate.sympt.inter.time")
  vax1.rr.infect <- get_param(dat, "vax1.rr.infect")
  vax2.rr.infect <- get_param(dat, "vax2.rr.infect")
  vax3.rr.infect <- get_param(dat, "vax3.rr.infect")
  vax4.rr.infect <- get_param(dat, "vax4.rr.infect")
  half.life <- get_param(dat, "half.life")

  nLayers <- length(dat$el)
  nInf <- rep(0, nLayers)

  if (length(idsInf) > 0) {
    for (layer in seq_len(nLayers)) {
      ## Look up discordant edgelist ##
      del <- discord_edgelist(dat, at, network = layer,
                              infstat = c("a", "ic", "ip"))
      del$dx <- dxStatus[del$inf]

      ## If any discordant pairs, proceed ##
      if (!(is.null(del))) {

        ## Parameters ##
        inf.prob <- get_param(dat, "inf.prob")[layer]
        act.rate <- get_param(dat, "act.rate")[layer]
        
        # Set parameters on discordant edgelist data frame
        del$transProb <- inf.prob
        del$actRate <- act.rate

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
        del$transProb <- del$transProb * (2 ^ (del$latest.vax / half.life))

        # Asymptomatic infection
        del$stat <- status[del$inf]
        del$transProb[del$stat == "a"] <- del$transProb[del$stat == "a"] *
                                          inf.prob.a.rr

        # Case isolation with diagnosed or symptomatic infection
        if (any(is.na(del$actRate[del$dx == 2]))) {browser()}
        if (at >= act.rate.dx.inter.time) {
          del$actRate[del$dx == 2] <- del$actRate[del$dx == 2] *
                                      act.rate.dx.inter.rr
        }
        if (at >= act.rate.sympt.inter.time) {
          del$actRate[del$stat == "ic"] <- del$actRate[del$stat == "ic"] *
                                           act.rate.sympt.inter.rr
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

  return(dat)
}
