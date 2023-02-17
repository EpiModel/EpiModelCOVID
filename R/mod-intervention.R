
#' @rdname moduleset-common
#' @export
contact_trace_covid <- function(dat, at) {

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
  # if (length(active) != length(status) & length(active) != length(dxStatus) &
  #     length(active)!= length(eligible.case)) browser()


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
  
  if (nEligCI > 0 & at >= inter.start.time) {
      
      ## Assign eligible case attribute for tracking later on ##
      eligible.case[idsEligCI] <- 1
      
      if (intervention == 1) {
        ## Look up discordant edgelist ##
        del_ct <- get_partners(dat, idsEligCI, only.active.nodes = TRUE, networks = 1)
        del_ct$index <- get_posit_ids(dat, del_ct$index)
        del_ct$partner <- get_posit_ids(dat, del_ct$partner)
      }
      
      if (intervention == 2) {
        ## Look up discordant edgelist ##
        del_ct <- get_partners(dat, idsEligCI, only.active.nodes = TRUE)
        del_ct$index <- get_posit_ids(dat, del_ct$index)
        del_ct$partner <- get_posit_ids(dat, del_ct$partner)
      }
   
      ## If any discordant pairs, proceed ##
      if (!(is.null(del_ct))) {

        # Set parameters on discordant edgelist data frame
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
          del_ct$iso.end[del_ct$status == 'ic'] <- pmax((del_ct$statusTime.Ic[del_ct$status == 'ic'] + 10),
                                                      del_ct$symendTime[del_ct$status == 'ic'],
                                                      na.rm = TRUE)
        }
        # comparing 10 days after symptom onset to symptom resolution to find max


        # Filter now for only contacts that have NOT been traced and for traced contacts
        #   FINISHED with their quarantine
        # Filter discordant edgelist for eligible contacts, assign new attribute for eligibility

        del_ct_traced <- del_ct[which(del_ct$traced.cc == 1), ,
                         drop = FALSE]
        del_ct <- del_ct[which(del_ct$traced.cc == 0 | is.na(del_ct$traced.cc)), ,
                         drop = FALSE]
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
        
        nIndex <- length(unique(del_ct$index))
        avg.partners <- nEligCT/nIndex

      ## Intervention 1: Tracing contacts at the cell level (tracing greater proportion of contacts & quarantining 100% of close contacts)
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
        ids.missing.quar <- which(del_ct$traced.cc == 1 & is.na(del_ct$quar))
        # cc.missing.quar <- del_ct$partner[ids.missing.quar]
        # if (length(ids.missing.quar < 1)) browser()
        num.missing.quar <- length(ids.missing.quar)
        if (num.missing.quar > 0) {
          del_ct$tracedTime[ids.missing.quar] <- at + baseline.lag
          del_ct$quarEnd[ids.missing.quar] <- del_ct$tracedTime[ids.missing.quar] + 14

          # Selecting sample of individuals to actually complete quarantine
          vec.quar.status <- rbinom(num.missing.quar, 1, 1)
          del_ct$quar[ids.missing.quar] <- vec.quar.status
        }
      }

      ## Intervention 2: Tracing contacts at both the cell and block levels (tracing fewer proportion of contacts & quarantining 50% of close contacts)
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
        ids.missing.quar <- which(del_ct$traced.cc == 1 & is.na(del_ct$quar))
        # cc.missing.quar <- del_ct$partner[ids.missing.quar]
        num.missing.quar <- length(ids.missing.quar)
        if (num.missing.quar > 0) {
          del_ct$tracedTime[ids.missing.quar] <- at + baseline.lag + time.lag
          del_ct$quarEnd[ids.missing.quar] <- del_ct$tracedTime[ids.missing.quar] + 14

          # Selecting sample of individuals to actually complete quarantine
          vec.quar.status <- rbinom(num.missing.quar, 1, 0.5)
          del_ct$quar[ids.missing.quar] <- vec.quar.status
        }
      }

      # Supplying empty set of contacts before interventions begin
      if (at < inter.start.time) {
        cc.missing.quar <- del_ct$partner[!is.na(del_ct$traced.cc) & is.na(del_ct$quar)]
      }

      if (at >= inter.start.time) {
        cc.missing.quar <- del_ct$partner[del_ct$traced.cc == 1 & is.na(del_ct$quar)]
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
      nQuar <- length(which(!is.na(del_ct$partner[del_ct$quar == 1])))
      ids.not.quar <- del_ct$partner[del_ct$quar == 0 & !is.na(del_ct$quar)]
      ids.quar.disc <- del_ct$partner[is.na(del_ct$quar)]

      # Check tracing and quarantining attributes for those already traced
      # Supplying empty set of contacts before interventions begin
      if (at < inter.start.time | length(del_ct_traced$partner) < 1) {
        ids.finished.quar <- del_ct_traced$partner[del_ct_traced$quarEnd < at]
      }

      if (length(del_ct_traced$partner) > 0 & at >= inter.start.time) {
        ids.still.quar <- del_ct_traced$partner[del_ct_traced$quarEnd >= at]

        ids.finished.quar <- del_ct_traced$partner[del_ct_traced$quarEnd < at]
      }
    }


      # if (length(ids.quar) < 1) browser()

      # Save updated attributes
        dat <- set_attr(dat, "eligible.case", eligible.case)
        dat <- set_attr(dat, "traced.cc", NA, ids.finished.quar)
        dat <- set_attr(dat, "traced.cc", 0, cc.not.traced)
        dat <- set_attr(dat, "traced.cc", 1, ids.traced)
        dat <- set_attr(dat, "quar", NA, ids.quar.disc)
        dat <- set_attr(dat, "quar", NA, ids.finished.quar)
        dat <- set_attr(dat, "quar", 0, ids.not.quar)
        dat <- set_attr(dat, "quar", 1, ids.quar)
        dat <- set_attr(dat, "tracedTime", at + baseline.lag, ids.traced)
        dat <- set_attr(dat, "tracedTime", NA, ids.finished.quar)
        dat <- set_attr(dat, "quarEnd", at + baseline.lag + 14, ids.traced)
        dat <- set_attr(dat, "quarEnd", NA, ids.finished.quar)

      ## Summary statistics
        dat <- set_epi(dat, "nQuar", at, nQuar) # number actually quarantining
        dat <- set_epi(dat, "nTraced", at, nTraced) # number of contacts traced
        dat <- set_epi(dat, "nElig.CC", at, nEligCT) # number of contacts eligible for tracing
        dat <- set_epi(dat, "avg.partners", at, avg.partners)

      }
    
  
  return(dat)

}
