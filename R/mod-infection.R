

#' @rdname moduleset-netjail
#' @export
infect_covid_netjail <- function(dat, at) {

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")
  statusTime <- get_attr(dat, "statusTime")
  quar <- get_attr(dat, "quar")
  tracedTime <- get_attr(dat, "tracedTime")
  quarEnd <- get_attr(dat, "quarEnd")
  transmissions <- dat$attr$transmissions
  age.grp <- get_attr(dat, "age.grp")
  vax <- get_attr(dat, "vax")

  ## Find infected nodes ##
  idsInf <- which(active == 1 & status %in% c("a", "ic", "ip"))
  nActive <- sum(active == 1)
  nElig <- length(idsInf)

  ## Common Parameters ##
  inf.prob.a.rr <- get_param(dat, "inf.prob.a.rr")
  act.rate.dx.inter.rr <- get_param(dat, "act.rate.dx.inter.rr")
  act.rate.dx.inter.time <- get_param(dat, "act.rate.dx.inter.time")
  act.rate.sympt.inter.rr <- get_param(dat, "act.rate.sympt.inter.rr")
  act.rate.sympt.inter.time <- get_param(dat, "act.rate.sympt.inter.time")
  act.rate.quar.inter.rr <- get_param(dat, "act.rate.quar.inter.rr")
  act.rate.quar.inter.time <- get_param(dat, "act.rate.quar.inter.time")
  vax.rr.infect.partial <- get_param(dat, "vax.rr.infect.partial")
  vax.rr.infect.full <- get_param(dat, "vax.rr.infect.full")
  com.inf.rate <- get_param(dat, "com.inf.rate")

  nLayers <- length(dat$el)
  nInf <- rep(0, nLayers)
  nAge <- rep(0, 9)

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
        # act.rate.quar.inter.rr <- get_param(dat, "act.rate.quar.inter.rr")[layer]
        # act.rate.quar.inter.time <- get_param(dat, "act.rate.quar.inter.time")[layer]
        
        # Set parameters on discordant edgelist data frame
        del$transProb <- inf.prob

        # Vaccination
        del$vaxSus <- vax[del$sus]
        del$transProb[del$vaxSus == 2] <- del$transProb[del$vaxSus == 2] *
                                          vax.rr.infect.partial
        del$transProb[del$vaxSus == 3] <- del$transProb[del$vaxSus == 3] *
                                          vax.rr.infect.full

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

        # Case isolation with diagnosed or symptomatic infection
        if (at >= act.rate.dx.inter.time) {
          del$actRate[del$dx == 2] <- del$actRate[del$dx == 2] *
            act.rate.dx.inter.rr
        }
        if (at >= act.rate.sympt.inter.time) {
          del$actRate[del$stat == "ic"] <- del$actRate[del$stat == "ic"] *
            act.rate.sympt.inter.rr
        }
        
        # Contact quarantine with tracing
        del$quar <- quar[del$sus]
        del$tracedTime <- tracedTime[del$sus]
        del$quarEnd <- quarEnd[del$sus]
        
        if (at >= act.rate.quar.inter.time) {
          del$actRate[del$quar == 1 & !is.na(del$quar) & 
                      at >= del$tracedTime & at <= del$quarEnd] <- del$actRate[del$quar == 1 & 
                                                                                 !is.na(del$quar) &
                                                                                 at >= del$tracedTime & 
                                                                                 at <= del$quarEnd] *
            act.rate.quar.inter.rr
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
  
  idsSus <- which(active == 1 & status %in% "s")
  idsInffromstaff <- which(rbinom(idsSus, 1, com.inf.rate) == 1)
  nInffromstaff <- length(idsInffromstaff)
  
  dat <- set_attr(dat, "status", "e", idsInffromstaff)
  dat <- set_attr(dat, "infTime", at, idsInffromstaff)
  dat <- set_attr(dat, "statusTime", at, idsInffromstaff)
  
  ## Summary statistics for incidence
  
  dat$epi$se.flow[at] <- sum(nInf) + nInffromstaff
  dat$epi$se.flow.staff[at] <- nInffromstaff
  dat$epi$se.flow.cell[at] <- nInf[1]
  dat$epi$se.flow.block[at] <- nInf[2]
  
  return(dat)
}
