#' @rdname moduleset-corporate
#' @export
isolate_covid_corporate <- function(dat, at) { # placed after the progression
  
  ## Attributes
  active     <- get_attr(dat, "active")
  status     <- get_attr(dat, "status")
  statusTime <- get_attr(dat, "statusTime")
  
  isolate <- get_attr(dat, "isolate")
  isoTime <- get_attr(dat, "isoTime")
  
  dxStatus <- get_attr(dat, "dxStatus")
  dxTime   <- get_attr(dat, "dxTime")
  
  hh.ids   <- get_attr(dat, "hh.ids")# labeled as household in Maria's work
  deg_work <- get_attr(dat, "deg_work")
  non.office <- ifelse(deg_work > 0, 0, 1)  # 0=office, 1=non-office
  
  ## Parameters
  iso.prob.office <- get_param(dat, "iso.prob.office")
  iso.prob.other  <- get_param(dat, "iso.prob.other")
  notif.prob      <- get_param(dat, "notif.prob")  # length 3: work/community/hh
  
  # ------------------------------------------------------------
  # Exposure(diagnostic status)-based isolation
  # ------------------------------------------------------------
  num.elig.exp <- 0
  num.new.iso4.other <- 0
  num.new.iso4.office <- 0
  
  ids.elig.exp <- which(active == 1 &
                          status %in% c("s","a","e","ip","r") &
                          is.na(isolate) &
                          dxStatus %in% 0:1 # 0 - not diagnosed, 1- negative
                        ) # ids for non-symptomatic, not diagnosed or neg, not isolating
  num.elig.exp <- length(ids.elig.exp)
  
  if (num.elig.exp > 0) {
    
    # pull partners of ids identified above 
    ids.exp <- get_partners(dat, ids.elig.exp, only.active.nodes = TRUE)
    ids.exp$index_posit_ids   <- get_posit_ids(dat, ids.exp$index)
    ids.exp$partner_posit_ids <- get_posit_ids(dat, ids.exp$partner)
    ids.exp$hh.index          <- hh.ids[ids.exp$index_posit_ids]
    
    # identify people who have been diagnosed pos in the past 3 time steps
    ids.dx <- which(active == 1 & dxStatus == 2 & dxTime >= at - 2)
    ids.dx.time <- data.frame(partner_posit_ids = ids.dx,
                              dxTime = dxTime[ids.dx],
                              hh.partner = hh.ids[ids.dx],
                              non.office = non.office[ids.dx]) 
    ids.exp.dx <- merge(ids.dx.time,ids.exp,by = "partner_posit_ids")
    
    # pull household members of those diagnosed pos in past 3 time steps
    ids.exp.dx.hh <- data.frame(
      partner_posit_ids = numeric(),
      dxTime = numeric(),
      hh.partner = numeric(),
      non.office = numeric(),
      index = numeric(),
      partner = numeric(),
      start = numeric(),
      stop = numeric(),
      network = numeric(),
      index_posit_ids = numeric(),
      hh.index = numeric()
    )
    for (i in ids.dx) { # each l represent exposed household member in a houshold
      l <- length(which(hh.ids==hh.ids[i]))
      new_row <- list(
        partner_posit_ids = rep(i,l),
        dxTime = rep(dxTime[i],l),
        hh.partner = rep(hh.ids[i],l),
        non.office = rep(non.office[i],l),
        index = get_unique_ids(dat,which(hh.ids==hh.ids[i])),
        partner = rep(get_unique_ids(dat,i),l),
        start = rep(1,l),
        stop = rep(NA,l),
        network = rep(3,l),
        index_posit_ids = which(hh.ids==hh.ids[i]), # household posit id where the positive case reside in
        hh.index = hh.ids[which(hh.ids==hh.ids[i])]
      )
      ids.exp.dx.hh <- rbind(ids.exp.dx.hh,new_row) # all household members include those not eligible (pos)
    }
    ids.exp.dx.hh <- ids.exp.dx.hh[ids.exp.dx.hh$index_posit_ids %in% ids.exp$index_posit_ids,] # eligible exposed houshold members
    
    # add household edges
    ids.exp.dx <- rbind(ids.exp.dx,ids.exp.dx.hh) # positive (left column), exposed (right column)

    # keep only contacts that occurred in recent window before the diagnosis, by network
    ids.exp.dx2 <- subset(ids.exp.dx,
                          (network == 3 & ids.exp.dx$dxTime == at-1) |
                            (network == 1 & ids.exp.dx$dxTime >= at-3) |
                            (network == 2 &
                               ids.exp.dx$dxTime >= ids.exp.dx$stop &
                               ids.exp.dx$dxTime <= ids.exp.dx$stop + 3))
        
    if (length(unique(ids.exp.dx2$index_posit_ids)) > 0) {
      # loop through networks
      for (j in 1:3) {
        ids.index <- unique(ids.exp.dx2$index_posit_ids[ids.exp.dx2$network == j])
        
        # identify those who were notified of the exposure by their contact
        vec.new.notif <- which(rbinom(length(ids.index),1,notif.prob[j]) == 1) # identify layer-specific contacts that contacted the recent positive case
        
        if (length(vec.new.notif) > 0) {
          ids.new.notif <- ids.index[vec.new.notif]
          
          # identify those who decide to follow isolation guidelines
          ### office workers
          ids.new.notif.office <- intersect(ids.new.notif, which(is.na(isolate) & non.office==0)) # posit ids of those need notification, not isolated and at office
          vec.new.iso4.office <- which(rbinom(length(ids.new.notif.office), 1, iso.prob.office) == 1)
          if (length(vec.new.iso4.office) > 0) {
            ids.new.iso4.office <- ids.new.notif.office[vec.new.iso4.office]
            num.new.iso4.office <- num.new.iso4.office + length(ids.new.iso4.office)
            isolate[ids.new.iso4.office] <- 4 # masking due to exposure notification
            
            partner.dx.time <- subset(ids.exp.dx2, ids.exp.dx2$index_posit_ids %in% ids.new.iso4.office) # recent exposed contacts who will mask because they were notified
            if (j == 2) {
              # for community contacts, start isolation at time of contact
              partner.dx.time <- partner.dx.time[,c("start",
                                                    "index_posit_ids")]
              partner.dx.time <- partner.dx.time[ave(partner.dx.time$start, # pick one isolate start time per exposed person then assign as isoTime
                                                     partner.dx.time$index_posit_ids,
                                                     FUN = function(x) x == min(x)) == 1, ]
              partner.dx.time <- unique(partner.dx.time)
              # if (length(ids.new.iso4.office) != length(partner.dx.time$start)) browser()
              isoTime[ids.new.iso4.office] <- partner.dx.time$start
            } else {
              # for HH and office contacts, start isolation at time of diagnosis
              # since contacts are 'permanent'
              partner.dx.time <- partner.dx.time[,c("dxTime",
                                                    "index_posit_ids")]
              partner.dx.time <- partner.dx.time[ave(partner.dx.time$dxTime,
                                                     partner.dx.time$index_posit_ids,
                                                     FUN = function(x) x == min(x)) == 1, ]
              partner.dx.time <- unique(partner.dx.time)
              # if (length(ids.new.iso4.office) != length(partner.dx.time$dxTime)) browser()
              isoTime[ids.new.iso4.office] <- partner.dx.time$dxTime
            }
            
          }
          ### non-office workers
          ids.new.notif.other <- intersect(ids.new.notif, which(is.na(isolate) & non.office==1))
          vec.new.iso4.other <- which(rbinom(length(ids.new.notif.other), 1, iso.prob.other) == 1)
          if (length(vec.new.iso4.other) > 0) {
            ids.new.iso4.other <- ids.new.notif.other[vec.new.iso4.other]
            num.new.iso4.other <- num.new.iso4.other + length(ids.new.iso4.other)
            isolate[ids.new.iso4.other] <- 4 # masking due to exposure notification
            
            partner.dx.time <- subset(ids.exp.dx2, ids.exp.dx2$index_posit_ids %in% ids.new.iso4.other)
            if (j == 2) {
              # for community contacts, start isolation at time of contact
              partner.dx.time <- partner.dx.time[,c("start",
                                                    "index_posit_ids")]
              partner.dx.time <- partner.dx.time[ave(partner.dx.time$start,
                                                     partner.dx.time$index_posit_ids,
                                                     FUN = function(x) x == min(x)) == 1, ]
              partner.dx.time <- unique(partner.dx.time)
              # if (length(ids.new.iso4.other) != length(partner.dx.time$start)) browser()
              isoTime[ids.new.iso4.other] <- partner.dx.time$start
            } else {
              # for HH and office contacts, start isolation at time of diagnosis
              # since contacts are 'permanent'
              partner.dx.time <- partner.dx.time[,c("dxTime",
                                                    "index_posit_ids")]
              partner.dx.time <- partner.dx.time[ave(partner.dx.time$dxTime,
                                                     partner.dx.time$index_posit_ids,
                                                     FUN = function(x) x == min(x)) == 1, ]
              partner.dx.time <- unique(partner.dx.time)
              # if (length(ids.new.iso4.other) != length(partner.dx.time$dxTime)) browser()
              isoTime[ids.new.iso4.other] <- partner.dx.time$dxTime
            }
            
          }
        }
      }
    }
    }
  
  # ------------------------------------------------------------
  # Clinical-progression-based isolation
  # ------------------------------------------------------------
  num.new.iso1.office <- 0
  num.new.iso1.other  <- 0
  
  # Start isolation for mild infection when status becomes ic
  ids.new.Ic <- which(active == 1 & status == "ic" & statusTime == at & is.na(isolate)) # Ic at current time point at 
  if (length(ids.new.Ic) > 0) {
    
    ids.new.Ic.other <- intersect(ids.new.Ic,which(non.office == 1))
    ids.new.Ic.office <- intersect(ids.new.Ic,which(non.office == 0))
    
    ### office workers
    vec.new.iso.office <- which(rbinom(length(ids.new.Ic.office), 1, iso.prob.office) == 1)
    if (length(vec.new.iso.office) > 0) {
      ids.new.iso1.office <- ids.new.Ic.office[vec.new.iso.office]
      num.new.iso1.office <- length(ids.new.iso1.office)
      isolate[ids.new.iso1.office] <- 1 # isolation for mild infection
      isoTime[ids.new.iso1.office] <- at
    }
    ### non-office workers
    vec.new.iso.other <- which(rbinom(length(ids.new.Ic.other), 1, iso.prob.other) == 1)
    if (length(vec.new.iso.other) > 0) {
      ids.new.iso1.other <- ids.new.Ic.other[vec.new.iso.other]
      num.new.iso1.other <- length(ids.new.iso1.other)
      isolate[ids.new.iso1.other] <- 1 # isolation for mild infection
      isoTime[ids.new.iso1.other] <- at
    }
  }
  
  # Escalate isolation for mild (1) to severe infection (2) for newly hospitalized (status becomes h)
  num.new.iso2  <- 0
  num.new.iso2.w <- 0
  
  ids.new.H <- which(active == 1 & status == "h" & statusTime == at)
  if (length(ids.new.H) > 0) {
    ids.new.iso2 <- intersect(ids.new.H, which(isolate == 1)) # for those who are hospitalized and isolated for mild infection update the isolation to 2
    num.new.iso2 <- length(ids.new.iso2)
    if (num.new.iso2 > 0) {
      isolate[ids.new.iso2] <- 2 # isolation for severe infection
      ids.new.iso2.w <- intersect(ids.new.iso2,which(non.office==0))
      num.new.iso2.w <- length(ids.new.iso2.w)
    }
  }
  
  # Move mild isolation (isolation =1) that started more than 5 day ago to post-isolation masking among R
  num.new.iso3 <- 0
  num.new.iso3.w <- 0
  
  ids.new.iso3 <- which(active == 1 & status == "r" & isolate == 1 & (at - isoTime) > 5)
  num.new.iso3 <- length(ids.new.iso3)
  if (num.new.iso3 > 0) {
    isolate[ids.new.iso3] <- 3 # post-isolation masking
    ids.new.iso3.w <- intersect(ids.new.iso3,which(non.office==0))
    num.new.iso3.w <- length(ids.new.iso3.w)
  }
  
  # For those who 1) recovered and isolated for severe infection or post-isolation masking, or 2) masked due to exposure notification and isolated more than 10 days, end isolation
  num.new.iso.end     <- 0
  num.new.iso.end.w   <- 0
  ids.new.iso.end <- which(active == 1 &
                             ((status == "r" & isolate %in% c(2,3)) | isolate == 4) &
                             (at - isoTime) > 10)
  num.new.iso.end <- length(ids.new.iso.end)
  if (num.new.iso.end > 0) {
    isolate[ids.new.iso.end] <- NA # end isolation pathway
    isoTime[ids.new.iso.end] <- NA
    ids.new.iso.end.w <- intersect(ids.new.iso.end,which(non.office==0))
    num.new.iso.end.w <- length(ids.new.iso.end.w)
  }
  
  ## Save updated attributes
  dat <- set_attr(dat, "isolate", isolate)
  dat <- set_attr(dat, "isoTime", isoTime)
  
  ## Save summary stats
  dat <- set_epi(dat, "iso1.flow", at, num.new.iso1.other + num.new.iso1.office)
  dat <- set_epi(dat, "iso2.flow", at, num.new.iso2)
  dat <- set_epi(dat, "iso3.flow", at, num.new.iso3)
  dat <- set_epi(dat, "iso4.flow", at, num.new.iso4.other + num.new.iso4.office)
  dat <- set_epi(dat, "isoend.flow", at, num.new.iso.end)

  dat <- set_epi(dat, "iso1.flow.w", at, num.new.iso1.office)
  dat <- set_epi(dat, "iso2.flow.w", at, num.new.iso2.w)
  dat <- set_epi(dat, "iso3.flow.w", at, num.new.iso3.w)
  dat <- set_epi(dat, "iso4.flow.w", at, num.new.iso4.office)
  dat <- set_epi(dat, "isoend.flow.w", at, num.new.iso.end.w)
  
  return(dat)
}
