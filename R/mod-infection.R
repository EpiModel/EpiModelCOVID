
#' @rdname moduleset-ship
#' @export
infect_covid_ship <- function(dat, at) {

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  transmissions <- get_attr(dat, "transmissions")
  type <- get_attr(dat, "type")

  ## Find infected nodes ##
  idsInf <- which(active == 1 & status %in% c("a", "ic", "ip"))

  ## Common Parameters ##
  inf.prob.a.rr <- get_param(dat, "inf.prob.a.rr")

  act.rate.dx.inter.rr <- get_param(dat, "act.rate.dx.inter.rr")
  act.rate.dx.inter.time <- get_param(dat, "act.rate.dx.inter.time")
  act.rate.sympt.inter.rr <- get_param(dat, "act.rate.sympt.inter.rr")
  act.rate.sympt.inter.time <- get_param(dat, "act.rate.sympt.inter.time")

  network.lockdown.time <- get_param(dat, "network.lockdown.time")
  if (at < network.lockdown.time) {
    nets <- 1:3
  } else {
    nets <- 4:6
  }

  ## Initialize default incidences at 0 ##
  nInf.PtoP <- 0
  nInf.PtoC <- 0
  nInf.CtoP <- 0
  nInf.CtoC <- 0

  # Pass/Pass Contacts
  if (length(idsInf) > 0) {

    ## Look up discordant edgelist ##
    del.PP <- discord_edgelist_covid_ship(dat, nw = nets[1])

    ## If any discordant pairs, proceed ##
    if (!(is.null(del.PP))) {

      ## Parameters ##
      inf.prob <- get_param(dat, "inf.prob.pp")
      act.rate <- get_param(dat, "act.rate.pp")
      inf.prob.pp.inter.rr <- get_param(dat, "inf.prob.pp.inter.rr")
      inf.prob.pp.inter.time <- get_param(dat, "inf.prob.pp.inter.time")
      act.rate.pp.inter.rr <- get_param(dat, "act.rate.pp.inter.rr")
      act.rate.pp.inter.time <- get_param(dat, "act.rate.pp.inter.time")

      # Set parameters on discordant edgelist data frame
      del.PP$transProb <- inf.prob
      del.PP$transProb[del.PP$stat == "a"] <- del.PP$transProb[del.PP$stat == "a"] *
                                              inf.prob.a.rr
      if (at >= inf.prob.pp.inter.time) {
        del.PP$transProb <- del.PP$transProb * inf.prob.pp.inter.rr
      }
      del.PP$actRate <- act.rate
      if (at >= act.rate.pp.inter.time) {
        del.PP$actRate <- del.PP$actRate * act.rate.pp.inter.rr
      }
      if (at >= act.rate.dx.inter.time) {
        del.PP$actRate[del.PP$dx == 2] <- del.PP$actRate[del.PP$dx == 2] *
                                          act.rate.dx.inter.rr
      }
      if (at >= act.rate.sympt.inter.time) {
        del.PP$actRate[del.PP$stat == "ic"] <- del.PP$actRate[del.PP$stat == "ic"] *
                                               act.rate.sympt.inter.rr
      }
      del.PP$finalProb <- 1 - (1 - del.PP$transProb)^del.PP$actRate

      # Stochastic transmission process
      transmit <- rbinom(nrow(del.PP), 1, del.PP$finalProb)

      # Keep rows where transmission occurred
      del.PP <- del.PP[which(transmit == 1), ]

      # Look up new ids if any transmissions occurred
      idsNewInf.PtoP <- unique(del.PP$sus)
      nInf.PtoP <- length(idsNewInf.PtoP)
      transIds <- del.PP$inf

      # Set new attributes for those newly infected
      if (nInf.PtoP > 0) {
        dat <- set_attr(dat, "status", "e", posit_ids = idsNewInf.PtoP)
        dat <- set_attr(dat, "infTime", at, posit_ids = idsNewInf.PtoP)
        dat <- set_attr(dat, "statusTime", at, posit_ids = idsNewInf.PtoP)
        trans_count <- get_attr(dat, "transmissions")
        for (tt in seq_along(transIds)) {
          trans_count[transIds[tt]] <- trans_count[transIds[tt]] + 1
        }
        dat <- set_attr(dat, "transmissions", trans_count)
      }
    }

    # Crew/Crew Contacts
    del.CC <- discord_edgelist_covid_ship(dat, nw = nets[2])
    if (!(is.null(del.CC))) {

      ## Parameters ##
      inf.prob <- get_param(dat, "inf.prob.cc")
      act.rate <- get_param(dat, "act.rate.cc")
      inf.prob.cc.inter.rr <- get_param(dat, "inf.prob.cc.inter.rr")
      inf.prob.cc.inter.time <- get_param(dat, "inf.prob.cc.inter.time")
      act.rate.cc.inter.rr <- get_param(dat, "act.rate.cc.inter.rr")
      act.rate.cc.inter.time <- get_param(dat, "act.rate.cc.inter.time")

      # Set parameters on discordant edgelist data frame
      del.CC$transProb <- inf.prob
      del.CC$transProb[del.CC$stat == "a"] <- del.CC$transProb[del.CC$stat == "a"] *
                                              inf.prob.a.rr
      if (at >= inf.prob.cc.inter.time) {
        del.CC$transProb <- del.CC$transProb * inf.prob.cc.inter.rr
      }
      del.CC$actRate <- act.rate
      if (at >= act.rate.cc.inter.time) {
        del.CC$actRate <- del.CC$actRate * act.rate.cc.inter.rr
      }
      if (at >= act.rate.dx.inter.time) {
        del.CC$actRate[del.CC$dx == 2] <- del.CC$actRate[del.CC$dx == 2] *
                                          act.rate.dx.inter.rr
      }
      if (at >= act.rate.sympt.inter.time) {
        del.CC$actRate[del.CC$stat == "ic"] <- del.CC$actRate[del.CC$stat == "ic"] *
                                               act.rate.sympt.inter.rr
      }
      del.CC$finalProb <- 1 - (1 - del.CC$transProb)^del.CC$actRate

      # Stochastic transmission process
      transmit <- rbinom(nrow(del.CC), 1, del.CC$finalProb)

      # Keep rows where transmission occurred
      del.CC <- del.CC[which(transmit == 1), ]

      # Look up new ids if any transmissions occurred
      idsNewInf.CtoC <- unique(del.CC$sus)
      nInf.CtoC <- length(idsNewInf.CtoC)
      transIds <- del.CC$inf

      # Set new attributes for those newly infected
      if (nInf.CtoC > 0) {
        dat <- set_attr(dat, "status", "e", posit_ids = idsNewInf.CtoC)
        dat <- set_attr(dat, "infTime", at, posit_ids = idsNewInf.CtoC)
        dat <- set_attr(dat, "statusTime", at, posit_ids = idsNewInf.CtoC)
        trans_count <- get_attr(dat, "transmissions")
        for (tt in seq_along(transIds)) {
          trans_count[transIds[tt]] <- trans_count[transIds[tt]] + 1
        }
        dat <- set_attr(dat, "transmissions", trans_count)
      }
    }

    # Pass/Crew Contacts
    del.PC <- discord_edgelist_covid_ship(dat, nw = nets[3])
    if (!(is.null(del.PC))) {

      ## Parameters ##
      inf.prob <- get_param(dat, "inf.prob.pc")
      act.rate <- get_param(dat, "act.rate.pc")
      inf.prob.pc.inter.rr <- get_param(dat, "inf.prob.pc.inter.rr")
      inf.prob.pc.inter.time <- get_param(dat, "inf.prob.pc.inter.time")
      act.rate.pc.inter.rr <- get_param(dat, "act.rate.pc.inter.rr")
      act.rate.pc.inter.time <- get_param(dat, "act.rate.pc.inter.time")

      # Set parameters on discordant edgelist data frame
      del.PC$transProb <- inf.prob
      del.PC$transProb[del.PC$stat == "a"] <- del.PC$transProb[del.PC$stat == "a"] *
                                              inf.prob.a.rr
      if (at >= inf.prob.pc.inter.time) {
        del.PC$transProb <- del.PC$transProb * inf.prob.pc.inter.rr
      }
      del.PC$actRate <- act.rate
      if (at >= act.rate.pc.inter.time) {
        del.PC$actRate <- del.PC$actRate * act.rate.pc.inter.rr
      }
      if (at >= act.rate.dx.inter.time) {
        del.PC$actRate[del.PC$dx == 2] <- del.PC$actRate[del.PC$dx == 2] *
                                          act.rate.dx.inter.rr
      }
      if (at >= act.rate.sympt.inter.time) {
        del.PC$actRate[del.PC$stat == "ic"] <- del.PC$actRate[del.PC$stat == "ic"] *
                                               act.rate.sympt.inter.rr
      }
      del.PC$finalProb <- 1 - (1 - del.PC$transProb)^del.PC$actRate

      # Stochastic transmission process
      transmit <- rbinom(nrow(del.PC), 1, del.PC$finalProb)

      # Keep rows where transmission occurred
      del.PC <- del.PC[which(transmit == 1), ]

      # New transmissions by direction
      idsNewInf.CtoP <- unique(del.PC$sus[type[del.PC$sus] == "p"])
      nInf.CtoP <- length(idsNewInf.CtoP)

      idsNewInf.PtoC <- unique(del.PC$sus[type[del.PC$sus] == "c"])
      nInf.PtoC <- length(idsNewInf.PtoC)

      # Either direction
      idsNewInf.PC <- union(idsNewInf.CtoP, idsNewInf.PtoC)
      transIds <- del.PC$inf

      # Set new attributes for those newly infected
      if ((nInf.CtoP + nInf.PtoC) > 0) {
        dat <- set_attr(dat, "status", "e", posit_ids = idsNewInf.PC)
        dat <- set_attr(dat, "infTime", at, posit_ids = idsNewInf.PC)
        dat <- set_attr(dat, "statusTime", at, posit_ids = idsNewInf.PC)
        trans_count <- get_attr(dat, "transmissions")
        for (tt in seq_along(transIds)) {
          trans_count[transIds[tt]] <- trans_count[transIds[tt]] + 1
        }
        dat <- set_attr(dat, "transmissions", trans_count)
      }
    }
  }

  ## Save summary statistics for S->E flow
  dat <- set_epi(dat, "se.flow", nInf.PtoP + nInf.PtoC + nInf.CtoP + nInf.CtoC, at = at)
  dat <- set_epi(dat, "se.pp.flow", nInf.PtoP, at = at)
  dat <- set_epi(dat, "se.pc.flow", nInf.PtoC, at = at)
  dat <- set_epi(dat, "se.cp.flow", nInf.CtoP, at = at)
  dat <- set_epi(dat, "se.cc.flow", nInf.CtoC, at = at)
  dat <- set_epi(dat, "Rt", mean(transmissions[status %in% c("a", "ip", "ic", "r")]), at = at)
  return(dat)
}


discord_edgelist_covid_ship <- function(dat, nw = 1) {

  status <- get_attr(dat, "status")
  dxStatus <- get_attr(dat, "dxStatus")

  el <- get_edgelist(dat, nw)

  del <- NULL
  if (nrow(el) > 0) {
    el <- el[sample(seq_len(nrow(el))), , drop = FALSE]
    stat <- matrix(status[el], ncol = 2)
    isInf <- matrix(stat %in% c("a", "ic", "ip"), ncol = 2)
    isSus <- matrix(stat %in% "s", ncol = 2)
    SIpairs <- el[isSus[, 1] * isInf[, 2] == 1, , drop = FALSE]
    ISpairs <- el[isSus[, 2] * isInf[, 1] == 1, , drop = FALSE]
    pairs <- rbind(SIpairs, ISpairs[, 2:1])
    if (nrow(pairs) > 0) {
      sus <- pairs[, 1]
      inf <- pairs[, 2]
      del <- data.frame(sus = sus, inf = inf, stat = status[inf], dx = dxStatus[inf])
    }
  }

  return(del)
}


#' @rdname moduleset-corporate
#' @export
infect_covid_corporate <- function(dat, at) {

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  vax <- get_attr(dat, "vax")

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

  nLayers <- dat$num.nw
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
        del$transProb[del$vaxSus %in% 2:3] <- del$transProb[del$vaxSus %in% 2:3] *
                                          vax1.rr.infect
        del$transProb[del$vaxSus == 4] <- del$transProb[del$vaxSus == 4] *
                                          vax2.rr.infect

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
