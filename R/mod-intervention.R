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
  tracedTime <- get_attr(dat, "tracedTime")
  quarEnd <- get_attr(dat, "quarEnd")
  
  
  ## Identify pool of eligible cases ##
  idsEligCI <- which(active == 1 & status %in% c("a", "ic", "ip") & 
                       dxStatus == 2 & is.na(eligible.case))
  nEligCI <- length(idsEligCI)
  
  ## Common Parameters ##
  prop.traced.1 <- get_param(dat, "prop.traced.1")
  prop.traced.2 <- get_param(dat, "prop.traced.2")
  time.lag <- get_param(dat, "time.lag")
  intervention <- get_param(dat, "intervention")
  
  if (nEligCI > 0) {
      
      ## Assign eligible case attribute for tracking later on ##
      eligible.case[idsEligCI] <- 1
    
      ## Look up discordant edgelist ##
      del_ct <- get_partners(dat, idsEligCI)
      
      ## If any discordant pairs, proceed ##
      if (!(is.null(del_ct))) {
        
        dxTime <- get_param(dat, "dxTime")
        statusTime.Ic <- get_param(dat, "statusTime.Ic")
        symendTime <- get_param(dat, "symendTime")
        
        # Set parameters on discordant edgelist data frame
        del_ct$dxTime <- dxTime[del_ct$index]
        del_ct$statusTime.Ic <- statusTime.Ic[del_ct$index]
        del_ct$symendTime <- symendTime[del_ct$index]
        del_ct$status <- status[del_ct$index]
        
        # Assign new isolation end attribute to discordant edgelist data frame
        del_ct$iso.end[del_ct$status %in% c('a', 'ip')] <- del_ct$dxTime[del_ct$status %in% c('a', 'ip')] + 10
        
        del_ct$iso.end[del_ct$status == 'ic'] <- max((del_ct$statusTime.Ic[del_ct$status == 'ic'] + 10),
                                                     del_ct$symendTime)
        # comparing 10 days after symptom onset to symptom resolution to find max
        
        # Filter discordant edgelist for eligible contacts, assign new attribute for eligibility
        del_ct$eligible.cc <- 0
        
        del_ct$eligible.cc[del_ct$status %in% c("a", "ip") & 
                             (del_ct$start <= del_ct$iso.end | 
                                del_ct$stop >= (del_ct$dxTime - 2))] = 1
        
        del_ct$eligible.cc[del_ct$status == "ip" & 
                             (del_ct$start <= del_ct$iso.end | 
                                del_ct$stop >= (del_ct$dxTime - 2))] = 1
      

        del_ct$eligible.cc[del_ct$status == 'ip'] <- ifelse(
          start <= iso.end | stop >= (statusTime.Ic - 2), 1, 0
        )
        
        # Keep only eligible close contacts
        del_ct <- del_ct[which(eligible.cc == 1), , drop = FALSE]
        nEligCT <- nrow(del_ct)
        
        ## Intervention 1: Varying fraction of traced contacts
        # Sample pool of eligible close contacts
        
        if (nEligCT > 0 & intervention == 1) {
          # Only sample group that has not already been traced
          ids.not.traced <- which(traced.cc != 1)
          num.not.traced <- length(ids.not.traced)
          if (num.not.traced > 0) {
            vec.traced.status <- rbinom(num.not.traced, 1, prop.traced.1)
            traced.cc[ids.not.traced] <- vec.traced.status
          }
          
          ids.traced <- which(traced.cc == 1)
          
          # Apply contact tracing attributes to close contacts
          ids.missing.quar <- which(is.na(quar))
          num.missing.quar <- length(ids.missing.quar)
          if (num.missing.quar > 0) {
            tracedTime[ids.missing.quar] <- at
            quarEnd[ids.missing.quar] <- tracedTime[ids.missing.quar] + 14
            
            # Selecting sample of individuals to actually complete quarantine
            vec.quar.status <- rbinom(num.missing.quar, 1, 0.8)
            quar[ids.missing.quar] <- vec.quar.status
          }
          
          # Check quarantine windows for those with quarantine attributes
          ids.with.quar <- which(!is.na(quar))
          num.with.quar <- length(ids.with.quar)
          if (num.with.quar > 0) {
            # for those still quarantining (quar = 1 or 0) keep current attributes
            ids.quar.cont <- which(quarEnd >= at)
            
            # for those finished with quarantine, transition back to missing quar
            ids.quar.discont <- which(quarEnd < at)
            num.quar.discont <- length(ids.quar.discont)
            if (num.quar.discont > 0) {
              quar[ids.quar.discont] <- NA
            }
          }
        }
        
        ## Intervention 2: Varying time to index case/close contact interview
        # Sample pool of eligible close contacts
        if (nEligCT > 0 & intervention == 2) {
          # Only sample group that has not already been traced
          ids.not.traced <- which(traced.cc != 1)
          num.not.traced <- length(ids.not.traced)
          if (num.not.traced > 0) {
            vec.traced.status <- rbinom(nrow(del_ct), 1, prop.traced.2)
            traced.cc[ids.not.traced] <- vec.traced.status
          }
          
          # Apply contact tracing attributes to close contacts
          ids.missing.quar <- which(is.na(quar))
          num.missing.quar <- length(ids.missing.quar)
          if (num.missing.quar > 0) {
            tracedTime[ids.missing.quar] <- at + time.lag
            quarEnd[ids.missing.quar] <- tracedTime[ids.missing.quar] + 14
            
            # Selecting sample of individuals to actually complete quarantine
            vec.quar.status <- rbinom(num.missing.quar, 1, 0.8)
            quar[ids.missing.quar] <- vec.quar.status
          }
          
          # Check quarantine windows for those with quarantine attributes
          ids.with.quar <- which(!is.na(quar))
          num.with.quar <- length(ids.with.quar)
          if (num.with.quar > 0) {
            # for all those still quarantining (quar = 1 or 0) keep current attributes
            ids.quar.cont <- which(quarEnd >= at & tracedTime <= at)
            
            # for those not started or already finished with quarantine, transition to missing quar
            ids.quar.discont <- which(at < tracedTime | quarEnd < at)
            num.quar.discont <- length(ids.quar.discont)
            if (num.quar.discont > 0) {
              quar[ids.quar.discont] <- NA
            }
          }
        }

        # Save updated attributes 
          dat <- set_attr(dat, "eligible.case", eligible.case)
          dat <- set_attr(dat, "traced.cc", traced.cc)
          dat <- set_attr(dat, "quar", quar)
          dat <- set_attr(dat, "iso.end", iso.end)
          dat <- set_attr(dat, "tracedTime", tracedTime)
          dat <- set_attr(dat, "quarEnd", quarEnd)

        ## Summary statistics
          dat <- set_epi(dat, "nQuar", at, length(ids.with.quar)) # number theoretically quarantining
          dat <- set_epi(dat, "nTraced", at, length(ids.traced)) # number of contacts traced
          dat <- set_epi(dat, "nElig.CC", at, nEligCT) # number of contacts eligible for tracing

    }
  }
  
  return(dat)
  
}
