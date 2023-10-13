
#' @rdname moduleset-common
#' @export
vax_covid <- function(dat, at) {
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  age <- get_attr(dat, "age")
  vax <- get_attr(dat, "vax")
  vax1Time <- get_attr(dat, "vax1Time")
  vax2Time <- get_attr(dat, "vax2Time")
  vax3Time <- get_attr(dat, "vax3Time")
  vax4Time <- get_attr(dat, "vax4Time")
  vax.age.group <- get_attr(dat, "vax.age.group")
  dxStatus <- get_attr(dat, "dxStatus")
  dxTime <- get_attr(dat, "dxTime")

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
  
  vax1.boost <- get_param(dat, "vax1.boost")
  vax2.boost <- get_param(dat, "vax2.boost")
  vax3.boost <- get_param(dat, "vax3.boost")
  vax4.boost <- get_param(dat, "vax4.boost")
  vax1.boost.start <- get_param(dat, "vax1.boost.start")
  vax2.boost.start <- get_param(dat, "vax2.boost.start")
  vax3.boost.start <- get_param(dat, "vax3.boost.start")
  vax4.boost.start <- get_param(dat, "vax4.boost.start")
  
  if (any(at == vax1.boost.start)) {
    idsElig.vax1.boost <- which(active == 1 & !(status %in% c("ic", "h"))
                            & !(dxStatus == 2 & (at - dxTime <= 10)) & vax == 0
                            & at == vax1.boost.start[vax.age.group])
    nElig.vax1.boost <- length(idsElig.vax1.boost)
    
    if (nElig.vax1.boost > 0) {
      vax1.boost.age <- vax1.boost[vax.age.group[idsElig.vax1.boost]]
      vec.vax1.boost <- which(rbinom(nElig.vax1.boost,1,vax1.boost.age) == 1)
      ids.vax1.boost <- idsElig.vax1.boost[vec.vax1.boost]
      
      if (length(ids.vax1.boost) > 0) {
        vax[ids.vax1.boost] <- 1
        vax1Time[ids.vax1.boost] <- at
      }
    }
  }
  
  if (any(at == vax2.boost.start)) {
    idsElig.vax2.boost <- which(active == 1 & !(status %in% c("ic", "h"))
                                & !(dxStatus == 2 & (at - dxTime <= 10)) & vax == 1
                                & at == vax2.boost.start[vax.age.group])
    nElig.vax2.boost <- length(idsElig.vax2.boost)
    
    if (nElig.vax2.boost > 0) {
      vax2.boost.age <- vax2.boost[vax.age.group[idsElig.vax2.boost]]
      vec.vax2.boost <- which(rbinom(nElig.vax2.boost,1,vax2.boost.age) == 1)
      ids.vax2.boost <- idsElig.vax2.boost[vec.vax2.boost]
      
      if (length(ids.vax2.boost) > 0) {
        vax[ids.vax2.boost] <- 1
        vax2Time[ids.vax2.boost] <- at
      }
    }
  }
  
  if (any(at == vax3.boost.start)) {
    idsElig.vax3.boost <- which(active == 1 & !(status %in% c("ic", "h"))
                                & !(dxStatus == 2 & (at - dxTime <= 10)) & vax == 2
                                & at == vax3.boost.start[vax.age.group])
    nElig.vax3.boost <- length(idsElig.vax3.boost)
    
    if (nElig.vax3.boost > 0) {
      vax3.boost.age <- vax3.boost[vax.age.group[idsElig.vax3.boost]]
      vec.vax3.boost <- which(rbinom(nElig.vax3.boost,1,vax3.boost.age) == 1)
      ids.vax3.boost <- idsElig.vax3.boost[vec.vax3.boost]
      
      if (length(ids.vax3.boost) > 0) {
        vax[ids.vax3.boost] <- 1
        vax3Time[ids.vax3.boost] <- at
      }
    }
  }
  
  if (any(at == vax4.boost.start)) {
    idsElig.vax4.boost <- which(active == 1 & !(status %in% c("ic", "h"))
                                & !(dxStatus == 2 & (at - dxTime <= 10)) & vax == 3
                                & at == vax4.boost.start[vax.age.group])
    nElig.vax4.boost <- length(idsElig.vax4.boost)
    
    if (nElig.vax4.boost > 0) {
      vax4.boost.age <- vax4.boost[vax.age.group[idsElig.vax4.boost]]
      vec.vax4.boost <- which(rbinom(nElig.vax4.boost,1,vax4.boost.age) == 1)
      ids.vax4.boost <- idsElig.vax4.boost[vec.vax4.boost]
      
      if (length(ids.vax4.boost) > 0) {
        vax[ids.vax4.boost] <- 1
        vax4Time[ids.vax4.boost] <- at
      }
    }
  }

  ## First vax
  nVax1 <- 0
  idsElig.vax1 <- which(active == 1 & !(status %in% c("ic", "h"))
                        & !(dxStatus == 2 & (at - dxTime <= 10)) & vax == 0
                        & at >= vax1.start[vax.age.group])
  nElig.vax1 <- length(idsElig.vax1)
  if (nElig.vax1 > 0) {
    #set vax rate based on age
    vax1.rate.vec <- vax1.rate[vax.age.group[idsElig.vax1]]

    vecVax1 <- which(rbinom(nElig.vax1, 1, vax1.rate.vec) == 1)
    idsVax1 <- idsElig.vax1[vecVax1]
    nVax1 <- length(idsVax1)
    if (nVax1 > 0) {
      vax[idsVax1] <- 1
      vax1Time[idsVax1] <- at
    }
  }

  ## Second vax
  nVax2 <- 0
  idsElig.vax2 <- which(active == 1 & !(status %in% c("ic", "h"))
                        & !(dxStatus == 2 & (at - dxTime <= 10)) & vax == 1
                        & (at - vax1Time >= vax2.interval))
  nElig.vax2 <- length(idsElig.vax2)
  if (nElig.vax2 > 0) {
    #set vax rate based on age
    vax2.rate.vec <- vax2.rate[vax.age.group[idsElig.vax2]]

    #roll for second dose
    vecVax2 <- which(rbinom(nElig.vax2, 1, vax2.rate.vec) == 1)
    idsVax2 <- idsElig.vax2[vecVax2]
    nVax2 <- length(idsVax2)

    #update
    if (nVax2 > 0) {
      vax[idsVax2] <- 2
      vax2Time[idsVax2] <- at
    }
  }

  ## Third vax / 1st booster
  nVax3 <- 0
  idsElig.vax3 <- which(active == 1 & !(status %in% c("ic", "h"))
                        & !(dxStatus == 2 & (at - dxTime <= 10)) & vax == 2
                        & (at - vax2Time >= vax3.interval)
                        & at >= vax3.start[vax.age.group])
  nElig.vax3 <- length(idsElig.vax3)
  if (nElig.vax3 > 0) {
    #set vax rate based on age
    vax3.rate.vec <- vax3.rate[vax.age.group[idsElig.vax3]]

    #roll for 3rd dose
    vecVax3 <- which(rbinom(nElig.vax3, 1, vax3.rate.vec) == 1)
    idsVax3 <- idsElig.vax3[vecVax3]
    nVax3 <- length(idsVax3)

    #update
    if (nVax3 > 0) {
      vax[idsVax3] <- 3
      vax3Time[idsVax3] <- at
    }
  }


  ## Fourth vax / 2nd booster
  nVax4 <- 0
  idsElig.vax4 <- which(active == 1 & !(status %in% c("ic", "h"))
                        & !(dxStatus == 2 & (at - dxTime <= 10)) & vax == 3
                        & (at - vax3Time >= vax4.interval)
                        & at >= vax4.start[vax.age.group])
  nElig.vax4 <- length(idsElig.vax4)
  if (nElig.vax4 > 0) {
    #set vax rate based on age
    vax4.rate.vec <- vax4.rate[vax.age.group[idsElig.vax4]]

    #roll for fourth dose
    vecVax4 <- which(rbinom(nElig.vax4, 1, vax4.rate.vec) == 1)
    idsVax4 <- idsElig.vax4[vecVax4]
    nVax4 <- length(idsVax4)

    #update
    if (nVax4 > 0) {
      vax[idsVax4] <- 4
      vax4Time[idsVax4] <- at
    }
  }

  ## Replace attr
  dat <- set_attr(dat, "vax", vax)
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


  return(dat)
}
