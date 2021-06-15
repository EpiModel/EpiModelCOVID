
#' @rdname moduleset-common
#' @export
dx_covid <- function(dat, at) {

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  dxStatus <- get_attr(dat, "dxStatus")

  dx.rate.sympt <- get_param(dat, "dx.rate.sympt")
  if (length(dx.rate.sympt) > 1) {
    dx.rate.sympt <- dx.rate.sympt[at]
  }
  dx.rate.other <- get_param(dat, "dx.rate.other")
  if (length(dx.rate.other) > 1) {
    dx.rate.other <- dx.rate.other[at]
  }
  allow.rescreen <- get_param(dat, "allow.rescreen")
  pcr.sens <- get_param(dat, "pcr.sens")

  idsDx.sympt <- idsDx.other <- NULL
  idsDx.sympt.pos <- idsDx.other.pos.true <- NULL
  idsDx.sympt.neg <- idsDx.other.pos.false <- NULL

  idsElig.sympt <- which(active == 1 & dxStatus %in% 0:1 & status == "ic")
  if (allow.rescreen == TRUE) {
    idsElig.other <- which(active == 1 & dxStatus %in% 0:1 &
                             status %in% c("s", "e", "a", "ip", "r"))
  } else {
    idsElig.other <- which(active == 1 & dxStatus == 0 &
                             status %in% c("s", "e", "a", "ip", "r"))
  }

  nElig.sympt <- length(idsElig.sympt)
  if (nElig.sympt > 0) {
    vecDx.sympt <- which(rbinom(nElig.sympt, 1, dx.rate.sympt) == 1)
    idsDx.sympt <- idsElig.sympt[vecDx.sympt]
    nDx.sympt <- length(idsDx.sympt)
    if (nDx.sympt > 0) {
      vecDx.sympt.pos <- rbinom(nDx.sympt, 1, pcr.sens)
      idsDx.sympt.pos <- idsDx.sympt[which(vecDx.sympt.pos == 1)]
      idsDx.sympt.neg <- idsDx.sympt[which(vecDx.sympt.pos == 0)]
      dxStatus[idsDx.sympt.pos] <- 2
      dxStatus[idsDx.sympt.neg] <- 1
    }
  }

  nElig.other <- length(idsElig.other)
  if (nElig.other > 0) {
    vecDx.other <- which(rbinom(nElig.other, 1, dx.rate.other) == 1)
    idsDx.other <- idsElig.other[vecDx.other]
    nDx.other <- length(idsDx.other)
    if (nDx.other > 0) {
      idsDx.other.neg <- intersect(idsDx.other, which(status == "s"))
      idsDx.other.pos.all <- intersect(idsDx.other,
                                       which(status %in% c("e", "a", "ip", "r")))
      vecDx.other.pos <- rbinom(length(idsDx.other.pos.all), 1, pcr.sens)
      idsDx.other.pos.true <- idsDx.other.pos.all[which(vecDx.other.pos == 1)]
      idsDx.other.pos.false <- idsDx.other.pos.all[which(vecDx.other.pos == 0)]
      dxStatus[idsDx.other.neg] <- 1
      dxStatus[idsDx.other.pos.false] <- 1
      dxStatus[idsDx.other.pos.true] <- 2
    }
  }

  ## Replace attr
  dat <- set_attr(dat, "dxStatus", dxStatus)

  ## Summary statistics ##
  dat <- set_epi(dat, "nDx", at, length(idsDx.sympt) + length(idsDx.other))
  dat <- set_epi(dat, "nDx.pos", at, length(idsDx.sympt.pos) +
                   length(idsDx.other.pos.true))
  dat <- set_epi(dat, "nDx.pos.sympt", at, length(idsDx.sympt.pos))
  dat <- set_epi(dat, "nDx.pos.fn", at, length(idsDx.sympt.neg) +
                   length(idsDx.other.pos.false))

  return(dat)
}


#' @rdname moduleset-common
#' @export
dx_covid_contacttrace <- function(dat, at) {
  
  ## Pull attributes
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  dxStatus <- get_attr(dat, "dxStatus")
  
  ## Pull parameters
  dx.rate.sympt <- get_param(dat, "dx.rate.sympt")
  if (length(dx.rate.sympt) > 1) {
    dx.rate.sympt <- dx.rate.sympt[at]
  }
  dx.rate.other <- get_param(dat, "dx.rate.other")
  if (length(dx.rate.other) > 1) {
    dx.rate.other <- dx.rate.other[at]
  }
  allow.rescreen <- get_param(dat, "allow.rescreen")
  pcr.sens <- get_param(dat, "pcr.sens")
  
  ## Initialize trackers
  idsDx.sympt <- idsDx.other <- NULL
  idsDx.sympt.pos <- idsDx.other.pos.true <- NULL
  idsDx.sympt.neg <- idsDx.other.pos.false <- NULL
  
  ## Determine screening eligibility
  idsElig.sympt <- which(active == 1 & dxStatus %in% 0:1 & status == "ic")
  if (allow.rescreen == TRUE) {
    idsElig.other <- which(active == 1 & dxStatus %in% 0:1 &
                             status %in% c("s", "e", "a", "ip", "r"))
  } else {
    idsElig.other <- which(active == 1 & dxStatus == 0 &
                             status %in% c("s", "e", "a", "ip", "r"))
  }
  
  ## Symptomatic testing
  nElig.sympt <- length(idsElig.sympt)
  if (nElig.sympt > 0) {
    vecDx.sympt <- which(rbinom(nElig.sympt, 1, dx.rate.sympt) == 1)
    idsDx.sympt <- idsElig.sympt[vecDx.sympt]
    nDx.sympt <- length(idsDx.sympt)
    if (nDx.sympt > 0) {
      vecDx.sympt.pos <- rbinom(nDx.sympt, 1, pcr.sens)
      idsDx.sympt.pos <- idsDx.sympt[which(vecDx.sympt.pos == 1)]
      idsDx.sympt.neg <- idsDx.sympt[which(vecDx.sympt.pos == 0)]
      dxStatus[idsDx.sympt.pos] <- 2
      dxStatus[idsDx.sympt.neg] <- 1
    }
  }
  
  ## Asymptomatic screening
  nElig.other <- length(idsElig.other)
  if (nElig.other > 0) {
    vecDx.other <- which(rbinom(nElig.other, 1, dx.rate.other) == 1)
    idsDx.other <- idsElig.other[vecDx.other]
    nDx.other <- length(idsDx.other)
    if (nDx.other > 0) {
      idsDx.other.neg <- intersect(idsDx.other, which(status == "s"))
      idsDx.other.pos.all <- intersect(idsDx.other,
                                       which(status %in% c("e", "a", "ip", "r")))
      vecDx.other.pos <- rbinom(length(idsDx.other.pos.all), 1, pcr.sens)
      idsDx.other.pos.true <- idsDx.other.pos.all[which(vecDx.other.pos == 1)]
      idsDx.other.pos.false <- idsDx.other.pos.all[which(vecDx.other.pos == 0)]
      dxStatus[idsDx.other.neg] <- 1
      dxStatus[idsDx.other.pos.false] <- 1
      dxStatus[idsDx.other.pos.true] <- 2
    }
  }
  
  ## Set attributes
  dat <- set_attr(dat, "dxStatus", dxStatus)
  
  ## Summary statistics
  dat <- set_epi(dat, "nDx", at, length(idsDx.sympt) + length(idsDx.other))
  dat <- set_epi(dat, "nDx.pos", at, length(idsDx.sympt.pos) +
                   length(idsDx.other.pos.true))
  dat <- set_epi(dat, "nDx.pos.sympt", at, length(idsDx.sympt.pos))
  dat <- set_epi(dat, "nDx.pos.fn", at, length(idsDx.sympt.neg) +
                   length(idsDx.other.pos.false))
  
  return(dat)
}
