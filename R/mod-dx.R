
#' @rdname moduleset-ship
#' @export
dx_covid_ship <- function(dat, at) {

  active <- dat$attr$active
  status <- dat$attr$status
  type <- dat$attr$type
  dxStatus <- dat$attr$dxStatus

  dx.start <- dat$param$dx.start
  dx.rate.pass <- dat$param$dx.rate.pass
  dx.rate.crew <- dat$param$dx.rate.crew
  dx.elig.status <- dat$param$dx.elig.status

  idsElig <- which(active == 1 & dxStatus == 0 & status %in% dx.elig.status)
  nElig <- length(idsElig)
  nDx <- 0
  nDx.pos <- 0
  nDx.pos.sympt <- 0

  if (at >= dx.start & nElig > 0) {
    dx.rates <- ifelse(type[idsElig] == "p", dx.rate.pass, dx.rate.crew)
    vecDx <- which(rbinom(nElig, 1, dx.rates) == 1)
    idsDx <- idsElig[vecDx]
    nDx <- length(idsDx)
    nDx.pos <- length(intersect(idsDx, which(status %in% c("e", "a", "ip", "ic"))))
    nDx.pos.sympt <- length(intersect(idsDx, which(status == "ic")))
    if (nDx > 0) {
      dxStatus[idsDx] <- 1
    }
  }

  ## Replace attr
  dat$attr$dxStatus <- dxStatus

  ## Summary statistics ##
  dat$epi$nDx[at] <- nDx
  dat$epi$nDx.pos[at] <- nDx.pos
  dat$epi$nDx.pos.sympt[at] <- nDx.pos.sympt

  return(dat)
}
