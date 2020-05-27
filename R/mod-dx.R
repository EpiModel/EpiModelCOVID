
#' @rdname moduleset-ship
#' @export
dx_covid_ship <- function(dat, at) {

  active <- dat$attr$active
  status <- dat$attr$status
  dxStatus <- dat$attr$dxStatus

  dx.rate.sympt <- dat$param$dx.rate.sympt[at]
  dx.rate.other <- dat$param$dx.rate.other[at]

  nDx.sympt <- nDx.other <- 0
  idsDx.other.pos <- NULL

  idsElig.sympt <- which(active == 1 & dxStatus %in% 0:1 & status == "ic")
  idsElig.other <- which(active == 1 & dxStatus %in% 0:1 & status %in% c("s", "e", "a", "ip", "r"))

  nElig.sympt <- length(idsElig.sympt)
  if (nElig.sympt > 0) {
    vecDx.sympt <- which(rbinom(nElig.sympt, 1, dx.rate.sympt) == 1)
    idsDx.sympt <- idsElig.sympt[vecDx.sympt]
    nDx.sympt <- length(idsDx.sympt)
    if (nDx.sympt > 0) {
      dxStatus[idsDx.sympt] <- 2
    }
  }

  nElig.other <- length(idsElig.other)
  if (nElig.other > 0) {
    vecDx.other <- which(rbinom(nElig.other, 1, dx.rate.other) == 1)
    idsDx.other <- idsElig.other[vecDx.other]
    nDx.other <- length(idsDx.other)
    if (nDx.other > 0) {
      idsDx.other.neg <- intersect(idsDx.other, which(status == "s"))
      idsDx.other.pos <- intersect(idsDx.other, which(status %in% c("e", "a", "ip", "r")))
      dxStatus[idsDx.other.neg] <- 1
      dxStatus[idsDx.other.pos] <- 2
    }
  }

  ## Replace attr
  dat$attr$dxStatus <- dxStatus

  ## Summary statistics ##
  dat$epi$nDx[at] <- nDx.sympt + nDx.other
  dat$epi$nDx.pos[at] <- nDx.sympt + length(idsDx.other.pos)
  dat$epi$nDx.pos.sympt[at] <- nDx.sympt

  return(dat)
}
