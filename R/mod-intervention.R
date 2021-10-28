#' @rdname moduleset-contacttrace
#' @export

intervention_covid_contacttrace <- function(dat, at) {
  
  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")
  statusTime.Ic <- get_attr(dat, "statusTime.Ic")
  statusTime <- get_attr(dat, "statusTime")
  dxTime <- get_attr(dat, "dxTime")
  dxStatus <- get_attr(dat, "dxStatus")
  eligible.case <- get_attr(dat, "eligible.case")
  
  ## Identify pool of eligible cases ##
  idsEligCI <- which(active == 1 & status %in% c("a", "ic", "ip") & 
                       dxStatus == 2 & is.na(eligible.case))
  nEligCI <- length(idsEligCI)
  
  ## Common Parameters ##
  inf.prob.a.rr <- get_param(dat, "inf.prob.a.rr")
  act.rate.dx.inter.rr <- get_param(dat, "act.rate.dx.inter.rr")
  act.rate.dx.inter.time <- get_param(dat, "act.rate.dx.inter.time")
  act.rate.sympt.inter.rr <- get_param(dat, "act.rate.sympt.inter.rr")
  act.rate.sympt.inter.time <- get_param(dat, "act.rate.sympt.inter.time")
  
  if (length(idsEligCI) > 0) {
      
      ## Assign eligible case attribute for tracking later on ##
      eligible.case[idsEligCI] <- 1
    
      ## Look up discordant edgelist ##
      del_ct <- get_partners(dat, idsEligCI, max.age = 2,
                             only.active = TRUE)
      
      ## If any discordant pairs, proceed ##
      if (!(is.null(del_ct))) {
        
        
        
        
        
        
        
        ## Parameters ##
        inf.prob <- get_param(dat, "inf.prob")[layer]
        act.rate <- get_param(dat, "act.rate")[layer]
        inf.prob.inter.rr <- get_param(dat, "inf.prob.inter.rr")[layer]
        inf.prob.inter.time <- get_param(dat, "inf.prob.inter.time")[layer]
        act.rate.inter.rr <- get_param(dat, "act.rate.inter.rr")[layer]
        act.rate.inter.time <- get_param(dat, "act.rate.inter.time")[layer]
        
        # Set parameters on discordant edgelist data frame
        del$transProb <- inf.prob
      
        
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
        
        # Set new attributes 
        if (nInf[layer] > 0) {
          dat <- set_attr(dat, "status", "e", idsNewInf)
          dat <- set_attr(dat, "infTime", at, idsNewInf)
          dat <- set_attr(dat, "statusTime", at, idsNewInf)
        }
      }
  }
  
}
