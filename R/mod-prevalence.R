
#' @rdname moduleset-ship
#' @export
prevalence_covid_ship <- function(dat, at) {

  active <- dat$attr$active
  status <- dat$attr$status
  type <- dat$attr$type
  dxStatus <- dat$attr$dxStatus

  nsteps <- dat$control$nsteps

  # Initialize Outputs
  var.names <- c("num", "s.num", "e.num", "a.num", "ip.num", "ic.num", "r.num",
                 "i.pass.num", "i.crew.num",
                 "dx.pos.num",
                 "se.flow", "ea.flow", "ar.flow", "Rt",
                 "eip.flow", "ipic.flow", "icr.flow",
                 "d.flow", "d.ic.flow", "exit.flow",
                 "nDx", "nDx.pos", "nDx.pos.sympt", "nDx.pos.fn",
                 "se.pp.flow", "se.pc.flow", "se.cp.flow", "se.cc.flow",
                 "meanAge", "meanClinic")
  if (at == 1) {
    for (i in seq_along(var.names)) {
      dat$epi[[var.names[i]]] <- rep(0, nsteps)
    }
  }

  # Update Outputs
  dat$epi$num[at] <- sum(active == 1)

  dat$epi$s.num[at] <- sum(active == 1 & status == "s")
  dat$epi$e.num[at] <- sum(active == 1 & status == "e")
  dat$epi$a.num[at] <- sum(active == 1 & status == "a")
  dat$epi$ip.num[at] <- sum(active == 1 & status == "ip")
  dat$epi$ic.num[at] <- sum(active == 1 & status == "ic")
  dat$epi$r.num[at] <- sum(active == 1 & status == "r")

  dat$epi$dx.pos.num[at] <- sum(active == 1 & dxStatus == 2, na.rm = TRUE)

  dat$epi$i.pass.num[at] <- sum(active == 1 & status %in% c("ip", "ic", "a") & type == "p")
  dat$epi$i.crew.num[at] <- sum(active == 1 & status %in% c("ip", "ic", "a") & type == "c")

  dat$epi$meanAge[at] <- mean(dat$attr$age, na.rm = TRUE)
  dat$epi$meanClinic[at] <- mean(dat$attr$clinical, na.rm = TRUE)

  return(dat)
}

#' @rdname moduleset-corporate
#' @export
prevalence_covid_corporate <- function(dat, at) {

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  vax <- get_attr(dat, "vax")

  # Initialize Outputs
  var.names <- c("num", "s.num", "e.num", "a.num", "ip.num", "ic.num", "r.num",
                 "h.num", "v1.num", "v2.num")
  if (at == 1) {
    for (i in seq_along(var.names)) {
      dat <- add_epi(dat, var.names[i])
    }
  }

  # Update Outputs
  dat <- set_epi(dat, "num", at, sum(active == 1))

  dat <- set_epi(dat, "s.num", at, sum(active == 1 & status == "s"))
  dat <- set_epi(dat, "e.num", at, sum(active == 1 & status == "e"))
  dat <- set_epi(dat, "a.num", at, sum(active == 1 & status == "a"))
  dat <- set_epi(dat, "ip.num", at, sum(active == 1 & status == "ip"))
  dat <- set_epi(dat, "ic.num", at, sum(active == 1 & status == "ic"))
  dat <- set_epi(dat, "r.num", at, sum(active == 1 & status == "r"))
  dat <- set_epi(dat, "h.num", at, sum(active == 1 & status == "h"))
  dat <- set_epi(dat, "v1.num", at, sum(active == 1 & status == "s" & vax == 1))
  dat <- set_epi(dat, "v2.num", at, sum(active == 1 & status == "s" & vax == 3))

  ## Person-denominated outputs (issue #46) --------------------------------------
  ##
  ## The model emitted only episode flows, which is why several validation
  ## comparisons could not be made at all: `cum_inf` derives from `se.flow`, so
  ## with reinfection on it counts episodes rather than people and a per-100
  ## above 100 cannot be set against any seroprevalence. These series are
  ## prevalences of a cumulative person-level flag, so they are attack rates by
  ## construction and are directly comparable to cohort and serosurvey numbers
  ## (CORES for COVID in Vellore is the local anchor).
  ##
  ## Emitted as counts with their own denominators rather than as rates, because
  ## departures shrink the living population over the run and the analysis needs
  ## to choose its denominator explicitly rather than inherit one.
  alive <- active == 1
  ever.inf <- get_attr(dat, "ever.inf")
  ever.sympt <- get_attr(dat, "ever.sympt")
  age <- get_attr(dat, "age")
  age.breaks <- get_param(dat, "age.breaks")
  nAge <- length(age.breaks) - 1
  ag <- cut(age, age.breaks, labels = FALSE, right = FALSE)

  if (!is.null(ever.inf)) {
    dat <- set_epi(dat, "ever.inf.num", at, sum(alive & ever.inf == 1, na.rm = TRUE))
    for (g in seq_len(nAge)) {
      dat <- set_epi(dat, paste0("ever.inf.age", g), at,
                     sum(alive & ag == g & ever.inf == 1, na.rm = TRUE))
    }
  }
  if (!is.null(ever.sympt)) {
    dat <- set_epi(dat, "ever.sympt.num", at, sum(alive & ever.sympt == 1, na.rm = TRUE))
    for (g in seq_len(nAge)) {
      dat <- set_epi(dat, paste0("ever.sympt.age", g), at,
                     sum(alive & ag == g & ever.sympt == 1, na.rm = TRUE))
    }
  }
  ## Per-band denominators. Without these the ever.inf.age* counts are not attack
  ## rates, and the age structure drifts over a 365-day run through ageing,
  ## arrivals and deaths, so a fixed t=1 denominator would be wrong.
  for (g in seq_len(nAge)) {
    dat <- set_epi(dat, paste0("n.age", g), at, sum(alive & ag == g, na.rm = TRUE))
  }

  ## Infant person-time, split at 6 months (issue #46). Infant N was never
  ## stored, so the infant series had no denominator at all and infant outcomes
  ## could only be reported per 100 of the WHOLE population, which buries the
  ## estimand. Summing these over time gives person-days; RSV severity is
  ## sharply front-loaded, hence the 0-6 / 6-12 split that matches
  ## dis.death.rate.infant / .infant6 and prop.hospit.infant / .infant6.
  is_infant <- get_attr(dat, "is_infant")
  if (!is.null(is_infant)) {
    inf_alive <- alive & is_infant
    dat <- set_epi(dat, "n.infant", at, sum(inf_alive, na.rm = TRUE))
    dat <- set_epi(dat, "n.infant.young", at, sum(inf_alive & age <  0.5, na.rm = TRUE))
    dat <- set_epi(dat, "n.infant.old",   at, sum(inf_alive & age >= 0.5, na.rm = TRUE))
    if (!is.null(ever.inf)) {
      dat <- set_epi(dat, "ever.inf.infant", at, sum(inf_alive & ever.inf == 1, na.rm = TRUE))
    }
  }

  return(dat)
}
