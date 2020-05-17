
#' @rdname moduleset-ship
#' @export
infect_covid_ship <- function(dat, at) {

  ## Attributes ##
  active <- dat$attr$active
  status <- dat$attr$status
  infTime <- dat$attr$infTime
  statusTime <- dat$attr$statusTime
  transmissions <- dat$attr$transmissions

  ## Find infected nodes ##
  idsInf <- which(active == 1 & status %in% c("a", "ic", "ip"))
  nActive <- sum(active == 1)
  nElig <- length(idsInf)

  ## Common Parameters ##
  inf.prob.a.rr <- dat$param$inf.prob.a.rr

  act.rate.dx.inter.rr <- dat$param$act.rate.dx.inter.rr
  act.rate.dx.inter.time <- dat$param$act.rate.dx.inter.time
  act.rate.sympt.inter.rr <- dat$param$act.rate.sympt.inter.rr
  act.rate.sympt.inter.time <- dat$param$act.rate.sympt.inter.time

  network.lockdown.time <- dat$param$network.lockdown.time
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
      inf.prob <- dat$param$inf.prob.pp
      act.rate <- dat$param$act.rate.pp
      inf.prob.pp.inter.rr <- dat$param$inf.prob.pp.inter.rr
      inf.prob.pp.inter.time <- dat$param$inf.prob.pp.inter.time
      act.rate.pp.inter.rr <- dat$param$act.rate.pp.inter.rr
      act.rate.pp.inter.time <- dat$param$act.rate.pp.inter.time

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
        del.PP$actRate[del.PP$dx == 1] <- del.PP$actRate[del.PP$dx == 1] *
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
        dat$attr$status[idsNewInf.PtoP] <- "e"
        dat$attr$infTime[idsNewInf.PtoP] <- at
        dat$attr$statusTime[idsNewInf.PtoP] <- at
        dat$attr$transmissions[transIds] <- transmissions[transIds] + 1
      }
    }

    # Crew/Crew Contacts
    del.CC <- discord_edgelist_covid_ship(dat, nw = nets[2])
    if (!(is.null(del.CC))) {

      ## Parameters ##
      inf.prob <- dat$param$inf.prob.cc
      act.rate <- dat$param$act.rate.cc
      inf.prob.cc.inter.rr <- dat$param$inf.prob.cc.inter.rr
      inf.prob.cc.inter.time <- dat$param$inf.prob.cc.inter.time
      act.rate.cc.inter.rr <- dat$param$act.rate.cc.inter.rr
      act.rate.cc.inter.time <- dat$param$act.rate.cc.inter.time

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
        del.CC$actRate[del.CC$dx == 1] <- del.CC$actRate[del.CC$dx == 1] *
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
        dat$attr$status[idsNewInf.CtoC] <- "e"
        dat$attr$infTime[idsNewInf.CtoC] <- at
        dat$attr$statusTime[idsNewInf.CtoC] <- at
        dat$attr$transmissions[transIds] <- transmissions[transIds] + 1
      }
    }

    # Pass/Crew Contacts
    del.PC <- discord_edgelist_covid_ship(dat, nw = nets[3])
    if (!(is.null(del.PC))) {

      ## Parameters ##
      inf.prob <- dat$param$inf.prob.pc
      act.rate <- dat$param$act.rate.pc
      inf.prob.pc.inter.rr <- dat$param$inf.prob.pc.inter.rr
      inf.prob.pc.inter.time <- dat$param$inf.prob.pc.inter.time
      act.rate.pc.inter.rr <- dat$param$act.rate.pc.inter.rr
      act.rate.pc.inter.time <- dat$param$act.rate.pc.inter.time

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
        del.PC$actRate[del.PP$dx == 1] <- del.PC$actRate[del.PC$dx == 1] *
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
      idsNewInf.CtoP <- unique(del.PC$sus[dat$attr$type[del.PC$sus] == "p"])
      nInf.CtoP <- length(idsNewInf.CtoP)

      idsNewInf.PtoC <- unique(del.PC$sus[dat$attr$type[del.PC$sus] == "c"])
      nInf.PtoC <- length(idsNewInf.PtoC)

      # Either direction
      idsNewInf.PC <- union(idsNewInf.CtoP, idsNewInf.PtoC)
      transIds <- del.PC$inf

      # Set new attributes for those newly infected
      if ((nInf.CtoP + nInf.PtoC) > 0) {
        dat$attr$status[idsNewInf.PC] <- "e"
        dat$attr$infTime[idsNewInf.PC] <- at
        dat$attr$statusTime[idsNewInf.PC] <- at
        dat$attr$transmissions[transIds] <- transmissions[transIds] + 1
      }
    }
  }


  ## Save summary statistics for S->E flow
  dat$epi$se.flow[at] <- nInf.PtoP + nInf.PtoC + nInf.CtoP + nInf.CtoC
  dat$epi$se.pp.flow[at] <- nInf.PtoP
  dat$epi$se.pc.flow[at] <- nInf.PtoC
  dat$epi$se.cp.flow[at] <- nInf.CtoP
  dat$epi$se.cc.flow[at] <- nInf.CtoC
  dat$epi$Rt[at] <- mean(transmissions[status %in% c("a", "ip", "ic", "r")])
  return(dat)
}

discord_edgelist_covid_ship <- function(dat, nw = 1) {

  status <- dat$attr$status
  dxStatus <- dat$attr$dxStatus
  type <- dat$attr$type

  el <- dat$el[[nw]]

  del <- NULL
  if (nrow(el) > 0) {
    el <- el[sample(1:nrow(el)), , drop = FALSE]
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
