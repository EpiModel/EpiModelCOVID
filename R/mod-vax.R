
#' @rdname moduleset-vaxDecisions
#' @export
vax_covid_vax_decisions <- function(dat, at) {
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  dxStatus <- get_attr(dat, "dxStatus")
  dxTime <- get_attr(dat, "dxTime")
  age <- get_attr(dat, "age")
  vax <- get_attr(dat, "vax")
  vaxSE <- get_attr(dat, "vaxSE")
  vax1Time <- get_attr(dat, "vax1Time")
  vax2Time <- get_attr(dat, "vax2Time")
  vax3Time <- get_attr(dat, "vax3Time")
  vax4Time <- get_attr(dat, "vax4Time")
  infTime <- get_attr(dat, "infTime")
  symptStartTime <- get_attr(dat, "symptStartTime")
  
  vax1.start <- get_param(dat, "vax1.start")
  vax2.interval <- get_param(dat, "vax2.interval")
  vax3.start <- get_param(dat, "vax3.start")
  vax3.interval <- get_param(dat, "vax3.interval")
  vax4.start <- get_param(dat, "vax4.start")
  vax4.interval <- get_param(dat, "vax4.interval")

  vax1.rate <- get_param(dat, "vax1.rate")
  vax2.rate <- get_param(dat, "vax2.rate")
  vax3.rate <- get_param(dat, "vax3.rate")
  vax4.rate <- get_param(dat, "vax4.rate")
  
  se.prob <- get_param(dat, "se.prob")
  se.rr <- get_param(dat, "se.rr")
  
  prevax.inf.rr <- get_param(dat, "prevax.inf.rr")
  postvax.inf.rr <- get_param(dat, "postvax.inf.rr")

  ## First vax
  nVax1 <- 0
  if (at >= vax1.start) {
    idsElig.vax1 <- which(active == 1 & !(status %in% c("ic", "h")) 
                          & !(dxStatus == 2 & (at - dxTime <= 10)) & vax == 0) 
    nElig.vax1 <- length(idsElig.vax1)
    if (nElig.vax1 > 0) {
      #set vax rate based on age and infection history
      age.group1 <- pmin((floor(age[idsElig.vax1] / 10)) + 1, 8)
      vax1.rate.vec <- vax1.rate[age.group1]
      locInf1 <- which(!is.na(symptStartTime[idsElig.vax1]))
      vax1.rate.vec[locInf1] <- vax1.rate.vec[locInf1] * prevax.inf.rr 

      #roll for initial vaccination
      vecVax1 <- which(rbinom(nElig.vax1, 1, vax1.rate.vec) == 1)
      idsVax1 <- idsElig.vax1[vecVax1]
      nVax1 <- length(idsVax1)
      
      #roll for side effects
      vecSE1 <- rbinom(nVax1, 1, se.prob[1]) 
      
      #update
      if (nVax1 > 0) {
        vax[idsVax1] <- 1
        vax1Time[idsVax1] <- at
        vaxSE[idsVax1] <- vecSE1
      }
    }
  }

  ## Second vax
  nVax2 <- 0
  idsElig.vax2 <- which(active == 1 & !(status %in% c("ic", "h")) 
                        & !(dxStatus == 2 & (at - dxTime <= 10)) & vax == 1 
                        & (at - vax1Time >= vax2.interval))
  nElig.vax2 <- length(idsElig.vax2)
  if (nElig.vax2 > 0) {
    #set vax rate based on age, infection history, and SE from dose 1
    age.group2 <- pmin((floor(age[idsElig.vax2] / 10)) + 1, 8)
    vax2.rate.vec <- vax2.rate[age.group2]
    
    locInf2 <- which(!is.na(symptStartTime[idsElig.vax2]))
    vax2.rate.vec[locInf2] <- vax2.rate.vec[locInf2] * prevax.inf.rr 
    
    locSE2 <- which(vaxSE[idsElig.vax2] == 1)
    vax2.rate.vec[locSE2] <- vax2.rate.vec[locSE2] * se.rr
    
    #roll for second dose
    vecVax2 <- which(rbinom(nElig.vax2, 1, vax2.rate.vec) == 1)
    idsVax2 <- idsElig.vax2[vecVax2]
    nVax2 <- length(idsVax2)
    
    #roll for side effects
    vecSE2 <- rbinom(nVax2, 1, se.prob[2]) 
    
    #update
    if (nVax2 > 0) {
      vax[idsVax2] <- 2
      vax2Time[idsVax2] <- at
      vaxSE[idsVax2] <- vecSE2
    }
  }
  
  ## Third vax / 1st booster
  nVax3 <- 0
  if (at >= vax3.start) {
    idsElig.vax3 <- which(active == 1 & !(status %in% c("ic", "h")) 
                          & !(dxStatus == 2 & (at - dxTime <= 10)) & vax == 2 
                          & (at - vax2Time >= vax3.interval)) 
    nElig.vax3 <- length(idsElig.vax3)
    if (nElig.vax3 > 0) {
      #set vax rate based on age, infection history, and SE from dose 2
      age.group3 <- pmin((floor(age[idsElig.vax3] / 10)) + 1, 8)
      vax3.rate.vec <- vax3.rate[age.group3]
      
      locInf3 <- which(symptStartTime[idsElig.vax3] >= 
                         pmax(vax2Time[idsElig.vax3] + vax3.interval, vax3.start))
      vax3.rate.vec[locInf3] <- vax3.rate.vec[locInf3] * prevax.inf.rr
      
      locBTInf3 <- which(symptStartTime[idsElig.vax3] < 
                         pmax(vax2Time[idsElig.vax3] + vax3.interval, vax3.start) 
                         & symptStartTime[idsElig.vax3] > vax2Time[idsElig.vax3])
      vax3.rate.vec[locBTInf3] <- vax3.rate.vec[locBTInf3] * postvax.inf.rr
      
      locSE3 <- which(vaxSE[idsElig.vax3] == 1)
      vax3.rate.vec[locSE3] <- vax3.rate.vec[locSE3] * se.rr
      
      #roll for 3rd dose
      vecVax3 <- which(rbinom(nElig.vax3, 1, vax3.rate.vec) == 1)
      idsVax3 <- idsElig.vax3[vecVax3]
      nVax3 <- length(idsVax3)
      
      #roll for side effects
      vecSE3 <- rbinom(nVax3, 1, se.prob[3]) 
      
      #update
      if (nVax3 > 0) {
        vax[idsVax3] <- 3
        vax3Time[idsVax3] <- at
        vaxSE[idsVax3] <- vecSE3
      }
    }
  }
  
  ## Fourth vax / 2nd booster
  nVax4 <- 0
  if (at >= vax4.start) {
    idsElig.vax4 <- which(active == 1 & !(status %in% c("ic", "h")) 
                          & !(dxStatus == 2 & (at - dxTime <= 10)) & vax == 3 
                          & (at - vax3Time >= vax4.interval)) 
    nElig.vax4 <- length(idsElig.vax4)
    if (nElig.vax4 > 0) {
      #set vax rate based on age, infection history, and SE from dose 3
      age.group4 <- pmin((floor(age[idsElig.vax4] / 10)) + 1, 8)
      vax4.rate.vec <- vax4.rate[age.group4]
      
      locInf4 <- which(symptStartTime[idsElig.vax4] >= 
                         pmax(vax3Time[idsElig.vax4] + vax4.interval, vax4.start))
      vax4.rate.vec[locInf4] <- vax4.rate.vec[locInf4] * prevax.inf.rr
      
      locBTInf4 <- which(symptStartTime[idsElig.vax4] < 
                           pmax(vax3Time[idsElig.vax4] + vax4.interval, vax4.start) 
                         & symptStartTime[idsElig.vax4] > vax3Time[idsElig.vax4])
      vax4.rate.vec[locBTInf4] <- vax4.rate.vec[locBTInf4] * postvax.inf.rr
      
      locSE4 <- which(vaxSE[idsElig.vax4] == 1)
      vax4.rate.vec[locSE4] <- vax4.rate.vec[locSE4] * se.rr
      
      #roll for fourth dose
      vecVax4 <- which(rbinom(nElig.vax4, 1, vax4.rate.vec) == 1)
      idsVax4 <- idsElig.vax4[vecVax4]
      nVax4 <- length(idsVax4)
      
      #roll for side effects
      vecSE4 <- rbinom(nVax4, 1, se.prob[4]) 
      
      #update
      if (nVax4 > 0) {
        vax[idsVax4] <- 4
        vax4Time[idsVax4] <- at
        vaxSE[idsVax4] <- vecSE4
      }
    }
  }

  ## Replace attr
  dat <- set_attr(dat, "vax", vax)
  dat <- set_attr(dat, "vaxSE", vaxSE)
  dat <- set_attr(dat, "vax1Time", vax1Time)
  dat <- set_attr(dat, "vax2Time", vax2Time)
  dat <- set_attr(dat, "vax3Time", vax3Time)
  dat <- set_attr(dat, "vax4Time", vax4Time)

  ## Summary statistics ##
  dat <- set_epi(dat, "nVax1", at, nVax1)
  dat <- set_epi(dat, "nVax2", at, nVax2)
  dat <- set_epi(dat, "nVax3", at, nVax3)
  dat <- set_epi(dat, "nVax4", at, nVax4)

  return(dat)
}
