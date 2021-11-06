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
  traced.cc <- get_attr(dat, "traced.cc")
  quar <- get_attr(dat, "quar")
  
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
  
  if (nEligCI > 0) {
      
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
                                ## same question for next 2 parameters
        del_ct$statusTime.Ic <- statusTime.Ic
        del_ct$symendTime <- symendTime
        del_ct$status <- status[del_ct$index]
        
        # Assign new isolattion end attribute to discordant edgelist data frame
        del_ct$iso.end[del_ct$status == 'a'] <- del_ct$dxTime[del_ct$status == 'a'] + 10
        
        del_ct$iso.end[del_ct$status == 'ip'] <- del_ct$dxTime[del_ct$status == 'ip'] + 10 
        
        del_ct$iso.end[del_ct$status == 'ic'] <- max((del_ct$statusTime.Ic[del_ct$status == 'ic'] + 10),
                                                     del_ct$symendTime)
        # comparing 10 days after symptom onset to symptom resolution to find max
        
        # Filter discordant edgelist for eligible contacts, assign new attribute for eligibility
        del_ct$eligible.cc[del_ct$status == 'a'] <- ifelse(
          start <= iso.end | stop >= (dxTime - 2), 1, 0
        )
        
        del_ct$eligible.cc[del_ct$status == 'ic'] <- ifelse(
          start <= iso.end | stop >= (dxTime - 2), 1, 0
        )
        
        del_ct$eligible.cc[del_ct$status == 'ip'] <- ifelse(
          start <= iso.end | stop >= (statusTime.Ic - 2), 1, 0
        )
        
        # Keep only eligible close contacts
        del_ct <- del_ct[which(eligible.cc == 1), , drop = FALSE]
        nEligCT <- nrow(del_ct)
        
        ## Intervention 1: Varying fraction of traced contacts
        # Sample pool of eligible close contacts
        if (nEligCT > 0) {
          
        
        
        
        
        }
        
        
        ## Intervention 2: Varying time to index case/close contact interview
        # Sample pool of eligible close contacts
        if (nEligCT > 0) {
          
          
          
          
          
          
          # Stochastic transmission process
          transmit <- rbinom(nrow(del), 1, del$finalProb)
          
          # Keep rows where transmission occurred
          del <- del[which(transmit == 1), , drop = FALSE]
          
          
          if (num.Ich > 0) {
            vec.new.H <- which(rbinom(num.Ich, 1, ich.rate) == 1)
            if (length(vec.new.H) > 0) {
              ids.new.H <- ids.Ich[vec.new.H]
              num.new.IctoH <- length(ids.new.H)
              status[ids.new.H] <- "h"
              statusTime[ids.new.H] <- at
              symendTime[ids.new.H] <- at
          
          
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
          dat <- set_attr(dat, "eligible.case", eligible.case)
          dat <- set_attr(dat, "traced.cc", traced.cc)
          dat <- set_attr(dat, "quar", quar)

        # Set new attributes for those newly traced
          if (nInf[layer] > 0) {
            dat <- set_attr(dat, "traced.cc", 1, idsNewInf)
            dat <- set_attr(dat, "infTime", at, idsNewInf)
            dat <- set_attr(dat, "statusTime", at, idsNewInf)
          }
          
          
        ## Summary statistics
          dat <- set_epi(dat, "nDx", at, length(idsDx.sympt) + length(idsDx.other))
          dat <- set_epi(dat, "nDx.pos", at, length(idsDx.sympt.pos) +
                           length(idsDx.other.pos.true))
          dat <- set_epi(dat, "nDx.pos.sympt", at, length(idsDx.sympt.pos))
          dat <- set_epi(dat, "nDx.pos.fn", at, length(idsDx.sympt.neg) +
                           length(idsDx.other.pos.false))
          
      }
  }
  
}
