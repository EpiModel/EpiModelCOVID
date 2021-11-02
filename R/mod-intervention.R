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
  iso.end <- get_attr(dat, "iso.end")
  eligible.cc <- get_attr(dat, "eligible.cc")
  
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
      del_ct <- get_partners(dat, idsEligCI)
      
      ## If any discordant pairs, proceed ##
      if (!(is.null(del_ct))) {
        
        dxTime <- get_param(dat, "dxTime")[idsEligCI]
        statusTime.Ic <- get_param(dat, "statusTime.Ic")[idsEligCI]
        symendTime <- get_param(dat, "symendTime")[idsEligCI]
        
        # Set parameters on discordant edgelist data frame
        del_ct$dxTime <- dxTime ## do I need to do a left_join here in case some 
                                ## of the index cases aren't in the discordant edgelist bc no partners?
        del_ct$statusTime.Ic <- statusTime.Ic
        del_ct$symendTime <- symendTime
        del_ct$status <- status[del_ct$index]
        
        # Assign new parameters to discordant edgelist data frame
        del_ct$iso.end <- del_ct$dxTime + 10
        
        # Filter discordant edgelist for eligible contacts
        
        
        
        
        
        
      

        
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
        
        # Save updated attributes 
          dat <- set_attr(dat, "eligible.case", eligible.case )
          dat <- set_attr(dat, "iso.end", iso.end)
          dat <- set_attr(dat, "eligible.cc", eligible.cc)

      }
  }
  
}
