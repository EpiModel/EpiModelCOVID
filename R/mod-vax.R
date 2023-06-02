
#' @rdname moduleset-vaxDecisions
#' @export
vax_covid_vax_decisions <- function(dat, at) {
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  dxStatus <- get_attr(dat, "dxStatus")
  dxTime <- get_attr(dat, "dxTime")
  age <- get_attr(dat, "age")
  vax <- get_attr(dat, "vax")
  vaxType <- get_attr(dat, "vaxType")
  vaxSE <- get_attr(dat, "vaxSE")
  vax1Time <- get_attr(dat, "vax1Time")
  vax2Time <- get_attr(dat, "vax2Time")
  vax3Time <- get_attr(dat, "vax3Time")
  vax4Time <- get_attr(dat, "vax4Time")
  symptStartTime <- get_attr(dat, "symptStartTime")
  vax.age.group <- get_attr(dat, "vax.age.group")
  
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
  
  hosp.nudge.prob <- get_param(dat, "hosp.nudge.prob")
  se.nudge.prob <- get_param(dat, "se.nudge.prob")
  bt.nudge.prob <- get_param(dat, "bt.nudge.prob")
  misc.nudge.prob.1 <- get_param(dat, "misc.nudge.prob.1")
  misc.nudge.prob.2 <- get_param(dat, "misc.nudge.prob.2")
  misc.nudge.prob.3 <- get_param(dat, "misc.nudge.prob.3")
  
  hosp.th <- get_param(dat, "hosp.th")

  num.new.w <- 0
  num.new.r <- 0
  hosp.flag <- 0
  
  ## Update vaccination types (if roll-out has begun)
  if (at > min(vax1.start)){
    
    # Check hospitalization trends
    h.7avg.num.prev <- get_epi(dat, "h.7avg.num", at - 2)
    h.7avg.num.curr <- get_epi(dat, "h.7avg.num", at - 1)
    prev.hosp.flag <- get_epi(dat, "hosp.flag")
    
    if (h.7avg.num.prev / sum(active == 1) < hosp.th &
        h.7avg.num.curr / sum(active == 1) >= hosp.th) {
      hosp.flag <- 1
    }
    
    # 1. Resistant -> Willing
    if (hosp.flag == 1 & sum(prev.hosp.flag, na.rm = TRUE) == 0) {
      idsElig1 <- which(active == 1 & vaxType == 0)
      
      # set prob rate based on age
      hosp.prob.vec <- hosp.nudge.prob[vax.age.group[idsElig1]]
      
      vaxType.new.1 <- rbinom(length(idsElig1), 1, hosp.prob.vec)
      num.new.w <- abs(sum(vaxType[idsElig1] - vaxType.new.1))
      vaxType[idsElig1] <- vaxType.new.1
    }

    # 2. Willing -> Resistant
    # 2a.1: 1 dose - side effects
    idsElig2a.1 <- which(active == 1 & vaxType == 1 & vax == 1
                       & vax1Time == (at - 1) & vaxSE == 1)
    se.prob.vec <- se.nudge.prob[vax.age.group[idsElig2a.1]]
    vaxType.new.2a.1 <- rbinom(length(idsElig2a.1), 1, (1 - se.prob.vec))
    num.new.r <- abs(sum(vaxType[idsElig2a.1] - vaxType.new.2a.1))
    vaxType[idsElig2a.1] <- vaxType.new.2a.1
    
    # 2a.2: 1 dose - misc.
    idsElig2a.2 <- which(active == 1 & vaxType == 1 & vax == 1 
                         & vax1Time == (at - 1))
    misc.prob.vec <- misc.nudge.prob.1[vax.age.group[idsElig2a.2]]
    vaxType.new.2a.2 <- rbinom(length(idsElig2a.2), 1, (1 - misc.prob.vec))
    num.new.r <- num.new.r + abs(sum(vaxType[idsElig2a.2] - vaxType.new.2a.2))
    vaxType[idsElig2a.2] <- vaxType.new.2a.2

    # 2b.1: 2 doses - side effects
    idsElig2b.1 <- which(active == 1 & vaxType == 1 & vax == 2
                         & vax2Time == (at - 1) & vaxSE == 1)
    se.prob.vec <- se.nudge.prob[vax.age.group[idsElig2b.1]]
    vaxType.new.2b.1 <- rbinom(length(idsElig2b.1), 1, (1 - se.prob.vec))
    num.new.r <- num.new.r + abs(sum(vaxType[idsElig2b.1] - vaxType.new.2b.1))
    vaxType[idsElig2b.1] <- vaxType.new.2b.1

    # 2b.2: 2 doses - breakthrough infections
    idsElig2b.2 <- which(active == 1 & vaxType == 1 & vax == 2
                         & symptStartTime == at
                         & at < pmax(vax2Time + vax3.interval,
                                     vax3.start[vax.age.group]))
    bt.prob.vec <- bt.nudge.prob[vax.age.group[idsElig2b.2]]
    vaxType.new.2b.2 <- rbinom(length(idsElig2b.2), 1, (1 - bt.prob.vec))
    num.new.r <- num.new.r + abs(sum(vaxType[idsElig2b.2] - vaxType.new.2b.2))
    vaxType[idsElig2b.2] <- vaxType.new.2b.2
    
    # 2b.3: 2 doses - misc.
    idsElig2b.3 <- which(active == 1 & vaxType == 1 & vax == 2
                         & vax2Time == (at - 1))
    misc.prob.vec <- misc.nudge.prob.2[vax.age.group[idsElig2b.3]]
    vaxType.new.2b.3 <- rbinom(length(idsElig2b.3), 1, (1 - misc.prob.vec))
    num.new.r <- num.new.r + abs(sum(vaxType[idsElig2b.3] - vaxType.new.2b.3))
    vaxType[idsElig2b.3] <- vaxType.new.2b.3

    # 2c.1: 3 doses - side effects
    idsElig2c.1 <- which(active == 1 & vaxType == 1 & vax == 3
                         & vax3Time == (at - 1) & vaxSE == 1)
    se.prob.vec <- se.nudge.prob[vax.age.group[idsElig2c.1]]
    vaxType.new.2c.1 <- rbinom(length(idsElig2c.1), 1, (1 - se.prob.vec))
    num.new.r <- num.new.r + abs(sum(vaxType[idsElig2c.1] - vaxType.new.2c.1))
    vaxType[idsElig2c.1] <- vaxType.new.2c.1

    # 2c.2: 3 doses - breakthrough infections
    idsElig2c.2 <- which(active == 1 & vaxType == 1 & vax == 3
                         & symptStartTime == at
                         & at < pmax(vax3Time + vax4.interval,
                                     vax4.start[vax.age.group]))
    bt.prob.vec <- bt.nudge.prob[vax.age.group[idsElig2c.2]]
    vaxType.new.2c.2 <- rbinom(length(idsElig2c.2), 1, (1 - bt.prob.vec))
    num.new.r <- num.new.r + abs(sum(vaxType[idsElig2c.2] - vaxType.new.2c.2))
    vaxType[idsElig2c.2] <- vaxType.new.2c.2

    # 2c.3: 3 doses - misc
    idsElig2c.3 <- which(active == 1 & vaxType == 1 & vax == 3
                         & vax3Time == (at - 1))
    misc.prob.vec <- misc.nudge.prob.3[vax.age.group[idsElig2c.3]]
    vaxType.new.2c.3 <- rbinom(length(idsElig2c.3), 1, (1 - misc.prob.vec))
    num.new.r <- num.new.r + abs(sum(vaxType[idsElig2c.3] - vaxType.new.2c.3))
    vaxType[idsElig2c.3] <- vaxType.new.2c.3
    
    # 2d.1: 4 doses - side effects
    idsElig2d.1 <- which(active == 1 & vaxType == 1 & vax == 4
                         & vax4Time == (at - 1) & vaxSE == 1)
        se.prob.vec <- se.nudge.prob[vax.age.group[idsElig2d.1]]
    
    vaxType.new.2d.1 <- rbinom(length(idsElig2d.1), 1, (1 - se.prob.vec))
    num.new.r <- num.new.r + abs(sum(vaxType[idsElig2d.1] - vaxType.new.2d.1))
    vaxType[idsElig2d.1] <- vaxType.new.2d.1

    # 2d.2: 4 doses - breakthrough infections
    idsElig2d.2 <- which(active == 1 & vaxType == 1 & vax == 4
                         & symptStartTime == at)
    bt.prob.vec <- bt.nudge.prob[vax.age.group[idsElig2d.2]]
    vaxType.new.2d.2 <- rbinom(length(idsElig2d.2), 1, (1 - bt.prob.vec))
    num.new.r <- num.new.r + abs(sum(vaxType[idsElig2d.2] - vaxType.new.2d.2))
    vaxType[idsElig2d.2] <- vaxType.new.2d.2
  }
  
  ## First vax
  nVax1 <- 0
  idsElig.vax1 <- which(active == 1 & !(status %in% c("ic", "h"))
                        & !(dxStatus == 2 & (at - dxTime <= 10)) & vax == 0
                        & (vaxType == 1 | is.na(vaxType)) 
                        & at >= vax1.start[vax.age.group])
  nElig.vax1 <- length(idsElig.vax1)
  if (nElig.vax1 > 0) {
    #set vax rate based on age
    vax1.rate.vec <- vax1.rate[vax.age.group[idsElig.vax1]]

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
  
  
  ## Second vax
  nVax2 <- 0
  idsElig.vax2 <- which(active == 1 & !(status %in% c("ic", "h"))
                        & !(dxStatus == 2 & (at - dxTime <= 10)) & vax == 1
                        & (at - vax1Time >= vax2.interval) & 
                          (vaxType == 1 | is.na(vaxType)))
  nElig.vax2 <- length(idsElig.vax2)
  if (nElig.vax2 > 0) {
    #set vax rate based on age
    vax2.rate.vec <- vax2.rate[vax.age.group[idsElig.vax2]]

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
  idsElig.vax3 <- which(active == 1 & !(status %in% c("ic", "h"))
                        & !(dxStatus == 2 & (at - dxTime <= 10)) & vax == 2
                        & (at - vax2Time >= vax3.interval) 
                        & (vaxType == 1 | is.na(vaxType))
                        & at >= vax3.start[vax.age.group])
  nElig.vax3 <- length(idsElig.vax3)
  if (nElig.vax3 > 0) {
    #set vax rate based on age
    vax3.rate.vec <- vax3.rate[vax.age.group[idsElig.vax3]]

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
  
  
  ## Fourth vax / 2nd booster
  nVax4 <- 0
  idsElig.vax4 <- which(active == 1 & !(status %in% c("ic", "h"))
                        & !(dxStatus == 2 & (at - dxTime <= 10)) & vax == 3
                        & (at - vax3Time >= vax4.interval) 
                        & (vaxType == 1 | is.na(vaxType))
                        & at >= vax4.start[vax.age.group])
  nElig.vax4 <- length(idsElig.vax4)
  if (nElig.vax4 > 0) {
    #set vax rate based on age
    vax4.rate.vec <- vax4.rate[vax.age.group[idsElig.vax4]]

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
  
  
  ## Replace attr
  dat <- set_attr(dat, "vax", vax)
  dat <- set_attr(dat, "vaxSE", vaxSE)
  dat <- set_attr(dat, "vaxType", vaxType)
  dat <- set_attr(dat, "vax1Time", vax1Time)
  dat <- set_attr(dat, "vax2Time", vax2Time)
  dat <- set_attr(dat, "vax3Time", vax3Time)
  dat <- set_attr(dat, "vax4Time", vax4Time)
  
  ## Summary statistics ##
  dat <- set_epi(dat, "nVax1", at, nVax1)
  dat <- set_epi(dat, "nVax2", at, nVax2)
  dat <- set_epi(dat, "nVax3", at, nVax3)
  dat <- set_epi(dat, "nVax4", at, nVax4)
  
  ## Summary stats -- vaccine coverage ##
  dat <- set_epi(dat, "cov_vax1_0to4", at, length(which(vax.age.group == 1 & vax >= 1)) / length(which(vax.age.group == 1)))
  dat <- set_epi(dat, "cov_vax1_5to17", at, length(which(vax.age.group == 2 & vax >= 1)) / length(which(vax.age.group == 2)))
  dat <- set_epi(dat, "cov_vax1_18to64", at, length(which((vax.age.group == 3 | vax.age.group == 4) & vax >= 1)) / length(which(vax.age.group == 3 | vax.age.group == 4)))
  dat <- set_epi(dat, "cov_vax1_65p", at, length(which(vax.age.group == 5 & vax >= 1)) / length(which(vax.age.group == 5)))

  dat <- set_epi(dat, "cov_vax2_0to4", at, length(which(vax.age.group == 1 & vax >= 2)) / length(which(vax.age.group == 1)))
  dat <- set_epi(dat, "cov_vax2_5to17", at, length(which(vax.age.group == 2 & vax >= 2)) / length(which(vax.age.group == 2)))
  dat <- set_epi(dat, "cov_vax2_18to64", at, length(which((vax.age.group == 3 | vax.age.group == 4) & vax >= 2)) / length(which(vax.age.group == 3 | vax.age.group == 4)))
  dat <- set_epi(dat, "cov_vax2_65p", at, length(which(vax.age.group == 5 & vax >= 2)) / length(which(vax.age.group == 5)))

  dat <- set_epi(dat, "cov_vax3_5to17", at, length(which(vax.age.group == 2 & vax >= 3)) / length(which(vax.age.group == 2)))
  dat <- set_epi(dat, "cov_vax3_18to49", at, length(which(vax.age.group == 3 & vax >= 3)) / length(which(vax.age.group == 3)))
  dat <- set_epi(dat, "cov_vax3_50to64", at, length(which(vax.age.group == 4 & vax >= 3)) / length(which(vax.age.group == 4)))
  dat <- set_epi(dat, "cov_vax3_65p", at, length(which(vax.age.group == 5 & vax >= 3)) / length(which(vax.age.group == 5)))

  dat <- set_epi(dat, "cov_vax4_50to64", at, length(which(vax.age.group == 4 & vax >= 4)) / length(which(vax.age.group == 4)))
  dat <- set_epi(dat, "cov_vax4_65p", at, length(which(vax.age.group == 5 & vax >= 4)) / length(which(vax.age.group == 5)))

  dat <- set_epi(dat, "num.new.w", at, num.new.w)
  dat <- set_epi(dat, "num.new.r", at, num.new.r)
  dat <- set_epi(dat, "prop.willing.1", at, length(which(vax.age.group == 3 & vaxType == 1)) / length(which(vax.age.group == 3)))
  dat <- set_epi(dat, "prop.willing.2", at, length(which(vax.age.group == 4 & vaxType == 1)) / length(which(vax.age.group == 4)))
  dat <- set_epi(dat, "prop.willing.3", at, length(which(vax.age.group == 5 & vaxType == 1)) / length(which(vax.age.group == 5)))
  
  dat <- set_epi(dat, "hosp.flag", at, hosp.flag)
  
  return(dat)
}
