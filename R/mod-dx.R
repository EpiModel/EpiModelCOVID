
#' @rdname moduleset-common
#' @export
dx_covid <- function(dat, at) {

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  dxStatus <- get_attr(dat, "dxStatus")
  dxTime <- get_attr(dat, "dxTime")
  isolate <- get_attr(dat, "isolate")
  isoTime <- get_attr(dat, "isoTime")
  non.office <- get_attr(dat, "non.office")

  iso.prob.office <- get_param(dat, "iso.prob.office")
  iso.prob.other <- get_param(dat, "iso.prob.other")

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

  #symptomatic IDs not on isolation pathway (use symptomatic test rate)
  idsElig.sympt <- which(active == 1 & dxStatus %in% 0:1 & status == "ic" &
                           is.na(isolate))
  #asymptomatic IDs not on isolation pathway (use symptomatic test rate)
  if (allow.rescreen == TRUE) {
    idsElig.other <- which(active == 1 & dxStatus %in% 0:1 &
                             status %in% c("s", "e", "a", "ip", "r") &
                             is.na(isolate))
  } else {
    idsElig.other <- which(active == 1 & dxStatus == 0 &
                             status %in% c("s", "e", "a", "ip", "r") &
                             is.na(isolate))
  }
  #IDs on isolation pathway --> all test
  # including those with positive dx so they can exit isolation if negative
  idsElig.iso <- which(active == 1 &
                         (isolate == 3 | (isolate == 4 & (at - isoTime) > 5)))

  # Testing for symptomatic (not currently isolating)
  nElig.sympt <- length(idsElig.sympt)
  num.new.dx.sympt.iso.office <- 0
  num.new.dx.sympt.iso.other <- 0
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
      dxTime[idsDx.sympt.pos] <- at
      dxTime[idsDx.sympt.neg] <- at

      # start isolation pathway for positive tests not already isolating
      ### office workers
      ids.new.dx.sympt.not.iso.office <- intersect(idsDx.sympt.pos, which(is.na(isolate) & non.office==0))
      vec.dx.sympt.iso.office <- which(rbinom(length(ids.new.dx.sympt.not.iso.office), 1, iso.prob.office) == 1)
      if (length(vec.dx.sympt.iso.office) > 0) {
        ids.new.dx.sympt.iso.office <- ids.new.dx.sympt.not.iso.office[vec.dx.sympt.iso.office]
        num.new.dx.sympt.iso.office <- length(ids.new.dx.sympt.iso.office)
        isolate[ids.new.dx.sympt.iso.office] <- 1
        isoTime[ids.new.dx.sympt.iso.office] <- at
      }
      ### non-office workers
      ids.new.dx.sympt.not.iso.other <- intersect(idsDx.sympt.pos, which(is.na(isolate) & non.office==1))
      vec.dx.sympt.iso.other <- which(rbinom(length(ids.new.dx.sympt.not.iso.other), 1, iso.prob.other) == 1)
      if (length(vec.dx.sympt.iso.other) > 0) {
        ids.new.dx.sympt.iso.other <- ids.new.dx.sympt.not.iso.other[vec.dx.sympt.iso.other]
        num.new.dx.sympt.iso.other <- length(ids.new.dx.sympt.iso.other)
        isolate[ids.new.dx.sympt.iso.other] <- 1
        isoTime[ids.new.dx.sympt.iso.other] <- at
      }
    }
  }

  # Testing for asymptomatic (not currently isolating)
  nElig.other <- length(idsElig.other)
  num.new.dx.iso.office <- 0
  num.new.dx.iso.other <- 0
  if (nElig.other > 0) {
    vecDx.other <- which(rbinom(nElig.other, 1, dx.rate.other) == 1)
    idsDx.other <- idsElig.other[vecDx.other]
    nDx.other <- length(idsDx.other)
    if (nDx.other > 0) {
      idsDx.other.neg <- intersect(idsDx.other, which(status == "s", "e", "r"))
      idsDx.other.pos.all <- intersect(idsDx.other,
                                       which(status %in% c("a", "ip")))
      vecDx.other.pos <- rbinom(length(idsDx.other.pos.all), 1, pcr.sens)
      idsDx.other.pos.true <- idsDx.other.pos.all[which(vecDx.other.pos == 1)]
      idsDx.other.pos.false <- idsDx.other.pos.all[which(vecDx.other.pos == 0)]
      dxStatus[idsDx.other.neg] <- 1
      dxStatus[idsDx.other.pos.false] <- 1
      dxStatus[idsDx.other.pos.true] <- 2
      dxTime[idsDx.other.neg] <- at
      dxTime[idsDx.other.pos.false] <- at
      dxTime[idsDx.other.pos.true] <- at

      # start isolation pathway for positive tests not already isolating
      ### office workers
      ids.new.dx.not.iso.office <- intersect(idsDx.other.pos.true, which(is.na(isolate) & non.office==0))
      vec.dx.iso.office <- which(rbinom(length(ids.new.dx.not.iso.office), 1, iso.prob.office) == 1)
      if (length(vec.dx.iso.office) > 0) {
        ids.new.dx.iso.office <- ids.new.dx.not.iso.office[vec.dx.iso.office]
        num.new.dx.iso.office <- length(ids.new.dx.iso.office)
        isolate[ids.new.dx.iso.office] <- 1
        isoTime[ids.new.dx.iso.office] <- at
      }
      ### non-office workers
      ids.new.dx.not.iso.other <- intersect(idsDx.other.pos.true, which(is.na(isolate) & non.office==1))
      vec.dx.iso.other <- which(rbinom(length(ids.new.dx.not.iso.other), 1, iso.prob.other) == 1)
      if (length(vec.dx.iso.other) > 0) {
        ids.new.dx.iso.other <- ids.new.dx.not.iso.other[vec.dx.iso.other]
        num.new.dx.iso.other <- length(ids.new.dx.iso.other)
        isolate[ids.new.dx.iso.other] <- 1
        isoTime[ids.new.dx.iso.other] <- at
      }
    }
  }

  # Testing for those on isolation pathway
  nElig.iso <- length(idsElig.iso)
  num.iso.end.iso <- 0
  if (nElig.iso > 0) {
    vecDx.iso <- which(rbinom(nElig.iso, 1, 1) == 1) #all test (for now)
    idsDx.iso <- idsElig.iso[vecDx.iso]
    nDx.iso <- length(idsDx.iso)
    if (nDx.iso > 0) {
      idsDx.iso.neg <- intersect(idsDx.iso, which(status == "s", "e", "r"))
      idsDx.iso.pos.all <- intersect(idsDx.iso,
                                       which(status %in% c("a", "ip", "ic")))
      vecDx.iso.pos <- rbinom(length(idsDx.iso.pos.all), 1, pcr.sens)
      idsDx.iso.pos.true <- idsDx.iso.pos.all[which(vecDx.iso.pos == 1)]
      idsDx.iso.pos.false <- idsDx.iso.pos.all[which(vecDx.iso.pos == 0)]
      dxStatus[idsDx.iso.neg] <- 1
      dxStatus[idsDx.iso.pos.false] <- 1
      dxStatus[idsDx.iso.pos.true] <- 2
      dxTime[idsDx.iso.neg] <- at
      dxTime[idsDx.iso.pos.false] <- at
      dxTime[idsDx.iso.pos.true] <- at

      # end isolation pathway for those who test negative
      ids.iso.end.iso <- intersect(c(idsDx.iso.neg,idsDx.iso.pos.false), which(!is.na(isolate)))
      if (length(ids.iso.end.iso) > 0) {
        num.iso.end.iso <- length(ids.iso.end.iso)
        isolate[ids.iso.end.iso] <- NA
        isoTime[ids.iso.end.iso] <- NA
      }
    }
  }

  ## Replace attr
  dat <- set_attr(dat, "dxStatus", dxStatus)
  dat <- set_attr(dat, "dxTime", dxTime)
  dat <- set_attr(dat, "isolate", isolate)
  dat <- set_attr(dat, "isoTime", isoTime)

  ## Summary statistics ##
  dat <- set_epi(dat, "nDx", at, length(idsDx.sympt) + length(idsDx.other))
  dat <- set_epi(dat, "nDx.pos", at, length(idsDx.sympt.pos) +
                   length(idsDx.other.pos.true))
  dat <- set_epi(dat, "nDx.pos.sympt", at, length(idsDx.sympt.pos))
  dat <- set_epi(dat, "nDx.pos.fn", at, length(idsDx.sympt.neg) +
                   length(idsDx.other.pos.false))
  dat <- set_epi(dat, "nDx.new.iso", at, num.new.dx.sympt.iso.office + num.new.dx.sympt.iso.other + num.new.dx.iso.office + num.new.dx.iso.other)
  dat <- set_epi(dat, "nDx.end.iso", at, num.iso.end.iso)

  return(dat)
}
