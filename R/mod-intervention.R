#' @rdname moduleset-contacttrace
#' @export

intervention_covid_contacttrace <- function(dat, at) {
  
  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")
  statusTime.Ic <- get_attr(dat, "statusTime.Ic")
  statusTime <- get_attr(dat, "statusTime")
  symendTime <- get_attr(dat, "symendTime")
  dxTime <- get_attr(dat, "dxTime")
  dxStatus <- get_attr(dat, "dxStatus")
  eligible.case <- get_attr(dat, "eligible.case")
  traced.cc <- get_attr(dat, "traced.cc")
  quar <- get_attr(dat, "quar")
  tracedTime <- get_attr(dat, "tracedTime")
  quarEnd <- get_attr(dat, "quarEnd")
  
  
  ## Identify pool of eligible cases ##
  if (length(active) != length(status) & length(active) != length(dxStatus) & 
      length(active)!= length(eligible.case)) browser()
  
  
  idsEligCI <- which(active == 1 & status %in% c("a", "ic", "ip") & 
                       dxStatus == 2 & is.na(eligible.case))
  nEligCI <- length(idsEligCI)
  
  ## Common Parameters ##
  prop.traced.1 <- get_param(dat, "prop.traced.1")
  prop.traced.2 <- get_param(dat, "prop.traced.2")
  time.lag <- get_param(dat, "time.lag")
  baseline.lag <- get_param(dat, "baseline.lag")
  intervention <- get_param(dat, "intervention")
  inter.start.time <- get_param(dat, "inter.start.time")
  
  if (nEligCI > 0) {
      
      #if (nEligCI > 0) browser()
    
      ## Assign eligible case attribute for tracking later on ##
      eligible.case[idsEligCI] <- 1
    
      ## Look up discordant edgelist ##
      del_ct <- get_partners(dat, get_posit_ids(dat, idsEligCI))
      
      ## If any discordant pairs, proceed ##
      if (!(is.null(del_ct))) {
        
        # dxTime <- get_param(dat, "dxTime")
        # statusTime.Ic <- get_param(dat, "statusTime.Ic")
        # symendTime <- get_param(dat, "symendTime")
        
        # Set parameters on discordant edgelist data frame
        # initialize columns to add from dat
        del_ct$dxTime <- 0
        del_ct$statusTime.Ic <- 0
        del_ct$symendTime <- 0
        del_ct$traced.cc <- 0
        del_ct$tracedTime <- 0
        del_ct$quar <- 0
        del_ct$quarEnd <- 0
        # consider putting in initialization for status column as 's'
        
        del_ct$dxTime <- dxTime[del_ct$index]
        del_ct$statusTime.Ic <- statusTime.Ic[del_ct$index]
        del_ct$symendTime <- symendTime[del_ct$index]
        del_ct$status <- status[del_ct$index]
        
        del_ct$traced.cc <- traced.cc[del_ct$partner]
        del_ct$tracedTime <- tracedTime[del_ct$partner]
        del_ct$quar <- quar[del_ct$partner]
        del_ct$quarEnd <- quarEnd[del_ct$partner]
        
        # Assign new isolation end attribute to discordant edgelist data frame
        # initialize iso.end column
        del_ct$iso.end <- 0
       
        del_ct$iso.end[del_ct$status %in% c('a', 'ip')] <- del_ct$dxTime[del_ct$status %in% c('a', 'ip')] + 10
        
        if (any(del_ct$status == 'ic')) {
          del_ct$iso.end[del_ct$status == 'ic'] <- max((del_ct$statusTime.Ic[del_ct$status == 'ic'] + 10),
                                                       del_ct$symendTime[del_ct$status == 'ic'],
                                                       na.rm = TRUE)
        } ## this is not working for ic individuals 
        # comparing 10 days after symptom onset to symptom resolution to find max
        
        # Filter discordant edgelist for eligible contacts, assign new attribute for eligibility
        del_ct$eligible.cc <- 0
        
        del_ct$eligible.cc[del_ct$status %in% c("a", "ip") & 
                             (del_ct$stop >= (del_ct$dxTime - 2) | is.na(del_ct$stop)) &
                             (del_ct$start <= del_ct$iso.end | 
                                del_ct$stop >= (del_ct$dxTime - 2))] <- 1
        
        del_ct$eligible.cc[del_ct$status == "ic" & 
                             (del_ct$stop >= (del_ct$statusTime.Ic - 2) | is.na(del_ct$stop)) &
                             (del_ct$start <= del_ct$iso.end | 
                                del_ct$stop >= (del_ct$statusTime.Ic - 2))] <- 1
        
        # Keep only eligible close contacts
        del_ct <- del_ct[which(del_ct$eligible.cc == 1), , drop = FALSE]
        nEligCT <- length(unique(del_ct$partner))
        
        ## Intervention 1: Varying fraction of traced contacts
        # Sample pool of eligible close contacts
        
        if (nEligCT > 0 & intervention == 1 & at >= inter.start.time) {
          
          # Only sample group that has not already been traced
          ids.not.traced <- which(del_ct$traced.cc == 0 | is.na(del_ct$traced.cc))
          num.not.traced <- length(ids.not.traced)
          if (num.not.traced > 0) {
            vec.traced.status <- rbinom(num.not.traced, 1, prop.traced.1)
            del_ct$traced.cc[ids.not.traced] <- vec.traced.status
          }
          
          # Apply contact tracing attributes to close contacts
          ids.missing.quar <- which(is.na(del_ct$quar))
          num.missing.quar <- length(ids.missing.quar)
          if (num.missing.quar > 0) {
            del_ct$tracedTime[ids.missing.quar] <- at + baseline.lag
            del_ct$quarEnd[ids.missing.quar] <- del_ct$tracedTime[ids.missing.quar] + 14
            
            # Selecting sample of individuals to actually complete quarantine
            vec.quar.status <- rbinom(num.missing.quar, 1, 0.8)
            del_ct$quar[ids.missing.quar] <- vec.quar.status
          }
        }
        
        ## Intervention 2: Varying time to index case/close contact interview
        # Sample pool of eligible close contacts
        if (nEligCT > 0 & intervention == 2 & at >= inter.start.time) {
          # Only sample group that has not already been traced
          ids.not.traced <- which(del_ct$traced.cc == 0 | is.na(del_ct$traced.cc))
          num.not.traced <- length(ids.not.traced)
          if (num.not.traced > 0) {
            vec.traced.status <- rbinom(num.not.traced, 1, prop.traced.2)
            del_ct$traced.cc[ids.not.traced] <- vec.traced.status
          }
          
          # Apply contact tracing attributes to close contacts
          ids.missing.quar <- which(is.na(del_ct$quar))
          cc.missing.quar <- del_ct$partner[is.na(del_ct$quar)]
          num.missing.quar <- length(ids.missing.quar)
          if (num.missing.quar > 0) {
            del_ct$tracedTime[ids.missing.quar] <- at + time.lag
            del_ct$quarEnd[ids.missing.quar] <- del_ct$tracedTime[ids.missing.quar] + 14
            
            # Selecting sample of individuals to actually complete quarantine
            vec.quar.status <- rbinom(num.missing.quar, 1, 0.8)
            del_ct$quar[ids.missing.quar] <- vec.quar.status
          }
        }

        ids.traced <- del_ct$partner[del_ct$traced.cc == 1 & !is.na(del_ct$traced.cc)]
        cc.not.traced <- del_ct$partner[del_ct$traced.cc == 0 & !is.na(del_ct$traced.cc)]
        
        nTraced <- length(which(!is.na(del_ct$partner[del_ct$traced.cc == 1])))
        nNotTraced <- length(which(!is.na(del_ct$partner[del_ct$traced.cc == 0])))
        
        # Check quarantine windows for those with quarantine attributes
        ids.with.quar <- which(!is.na(del_ct$quar))
        num.with.quar <- length(ids.with.quar)
        if (num.with.quar > 0) {
          # for those still quarantining (quar = 1 or 0) keep current attributes
          ids.quar.current <- which(at <= quarEnd & at >= tracedTime)

          # for those who haven't started quarantining yet, keep quarantine status
          ids.quar.future <- which(at < del_ct$tracedTime)
          num.quar.future <- length(ids.quar.future)
          
          # for those finished with quarantine, transition back to missing quar
          ids.quar.discont <- which(del_ct$quarEnd < at)
          num.quar.discont <- length(ids.quar.discont)
          if (num.quar.discont > 0) {
            del_ct$quar[ids.quar.discont] <- NA
          }
        }
        
        ids.quar <- del_ct$partner[del_ct$quar == 1 & !is.na(del_ct$quar)]
        ids.not.quar <- del_ct$partner[del_ct$quar == 0 & !is.na(del_ct$quar)]
        
        
        # Save updated attributes 
          dat <- set_attr(dat, "eligible.case", eligible.case)
          dat <- set_attr(dat, "traced.cc", 0, get_posit_ids(dat, cc.not.traced))
          dat <- set_attr(dat, "traced.cc", 1, get_posit_ids(dat, ids.traced))
          dat <- set_attr(dat, "quar", 0, get_posit_ids(dat, ids.not.quar))
          dat <- set_attr(dat, "quar", 1, get_posit_ids(dat, ids.quar))
          dat <- set_attr(dat, "tracedTime", at + time.lag, get_posit_ids(dat, cc.missing.quar))
          dat <- set_attr(dat, "quarEnd", tracedTime + 14, get_posit_ids(dat, cc.missing.quar))

        ## Summary statistics
          dat <- set_epi(dat, "nQuar", at, length(ids.with.quar)) # number theoretically quarantining
          dat <- set_epi(dat, "nTraced", at, nTraced) # number of contacts traced
          dat <- set_epi(dat, "nElig.CC", at, nEligCT) # number of contacts eligible for tracing

    }
  }
  
  return(dat)
  
}
