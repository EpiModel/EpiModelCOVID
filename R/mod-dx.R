
#' @rdname moduleset-ship
#' @export
dx_covid_ship <- function(dat, at) {

  active <- dat$attr$active
  status <- dat$attr$status
  dxStatus <- dat$attr$dxStatus

  dx.rate.sympt <- dat$param$dx.rate.sympt[at]
  dx.rate.other <- dat$param$dx.rate.other[at]
  allow.rescreen <- dat$param$allow.rescreen
  pcr.sens <- dat$param$pcr.sens

  idsDx.sympt <- idsDx.other <- NULL
  idsDx.sympt.pos <- idsDx.other.pos.true <- NULL
  idsDx.sympt.neg <- idsDx.other.pos.false <- NULL

  idsElig.sympt <- which(active == 1 & dxStatus %in% 0:1 & status == "ic")
  if (allow.rescreen == TRUE) {
    idsElig.other <- which(active == 1 & dxStatus %in% 0:1 & status %in% c("s", "e", "a", "ip", "r"))
  } else {
    idsElig.other <- which(active == 1 & dxStatus == 0 & status %in% c("s", "e", "a", "ip", "r"))
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
      idsDx.other.pos.all <- intersect(idsDx.other, which(status %in% c("e", "a", "ip", "r")))
      vecDx.other.pos <- rbinom(length(idsDx.other.pos.all), 1, pcr.sens)
      idsDx.other.pos.true <- idsDx.other.pos.all[which(vecDx.other.pos == 1)]
      idsDx.other.pos.false <- idsDx.other.pos.all[which(vecDx.other.pos == 0)]
      dxStatus[idsDx.other.neg] <- 1
      dxStatus[idsDx.other.pos.false] <- 1
      dxStatus[idsDx.other.pos.true] <- 2
    }
  }

  ## Replace attr
  dat$attr$dxStatus <- dxStatus

  ## Summary statistics ##
  dat$epi$nDx[at] <- length(idsDx.sympt) + length(idsDx.other)
  dat$epi$nDx.pos[at] <- length(idsDx.sympt.pos) + length(idsDx.other.pos.true)
  dat$epi$nDx.pos.sympt[at] <- length(idsDx.sympt.pos)
  dat$epi$nDx.pos.fn[at] <- length(idsDx.sympt.neg) + length(idsDx.other.pos.false)

  return(dat)
}


#' @rdname moduleset-corporate
#' @export
dx_covid_corporate <- function(dat, at) {

  active <- dat$attr$active
  status <- dat$attr$status
  dxStatus <- dat$attr$dxStatus

  dx.rate.sympt <- dat$param$dx.rate.sympt[at]
  dx.rate.other <- dat$param$dx.rate.other[at]
  allow.rescreen <- dat$param$allow.rescreen
  pcr.sens <- dat$param$pcr.sens

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
      idsDx.other.pos.all <- intersect(idsDx.other, which(status %in% c("e", "a", "ip", "r")))
      vecDx.other.pos <- rbinom(length(idsDx.other.pos.all), 1, pcr.sens)
      idsDx.other.pos.true <- idsDx.other.pos.all[which(vecDx.other.pos == 1)]
      idsDx.other.pos.false <- idsDx.other.pos.all[which(vecDx.other.pos == 0)]
      dxStatus[idsDx.other.neg] <- 1
      dxStatus[idsDx.other.pos.false] <- 1
      dxStatus[idsDx.other.pos.true] <- 2
    }
  }

  ## Replace attr
  dat$attr$dxStatus <- dxStatus

  ## Summary statistics ##
  dat$epi$nDx[at] <- length(idsDx.sympt) + length(idsDx.other)
  dat$epi$nDx.pos[at] <- length(idsDx.sympt.pos) + length(idsDx.other.pos.true)
  dat$epi$nDx.pos.sympt[at] <- length(idsDx.sympt.pos)
  dat$epi$nDx.pos.fn[at] <- length(idsDx.sympt.neg) + length(idsDx.other.pos.false)

  return(dat)
}
