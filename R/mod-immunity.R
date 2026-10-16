# Episodes, post-infection immunity, prior immunity and per-infection severity
# (GlobalMix India #69, #70, #71)
#
# Every infection is an episode with its own severity. start_episode() clears an
# episode's clinical, hospit and fatal draws, so a reinfection draws them again
# with the person's current age, vaccine and immunity; they used to be drawn
# once per person and replayed on every reinfection (#69).
#
# Immunity after an infection in the run (#69). A person who recovers is fully
# protected for immune.full days (status r; Inf keeps them there for the rest
# of the run), then returns to s with partial protection that wanes from its
# value at recovery like the vaccine's: floor + (peak - floor) x 0.5^(d /
# halflife), d the days since recovery. The protections are marginal, against
# infection, symptomatic disease, admission and death, as studies of
# reinfection report them, and are converted to the per-stage relative risks
# the model draws with: infection per contact 1 - P_infect, symptomatic given
# infection (1 - P_clinical) / (1 - P_infect), and admission and death given
# symptomatic disease (1 - P_hosp) / (1 - P_clinical) and (1 - P_death) /
# (1 - P_clinical).
#
# Prior immunity (#70). init.immune is the share of each severity band immune
# at the start of the run (status r for the whole run, prior.immune = 1). Who
# is immune within a band is drawn with probability min(1, c x w), w the
# person's propensity in the attribute named by init.immune.attr and c set so
# the band's expected share is its target; with no attribute, everyone in the
# band has the same probability.
#
# Severity per infection (#71). sev.age.breaks are the severity bands, finer
# than the network bands. prop.clinical, ihr and ifr are, by severity band, the
# probabilities per infection of symptomatic disease, admission and death of an
# unvaccinated person with no immunity. A symptomatic case is admitted with
# probability ihr / prop.clinical and dies with probability ifr /
# prop.clinical, drawn once at onset, so a change to prop.clinical leaves the
# IFR and IHR per infection where they are. death.hosp.rr is the relative risk
# of death of an admitted case against one not admitted (1: independent of
# admission), with the death probability given symptomatic disease unchanged.

# Severity bands --------------------------------------------------------------------

severity_params <- c("prop.clinical", "ihr", "ifr")

# Severity band of each age, 1 to the number of bands (right-open bands)
sev_band <- function(dat, age) {
  findInterval(age, get_param(dat, "sev.age.breaks"))
}

check_severity_params <- function(dat) {
  br <- get_param(dat, "sev.age.breaks", override.null.error = TRUE)
  if (is.null(br) || !is.numeric(br) || length(br) < 2 || br[1] != 0 || !is.infinite(br[length(br)]) ||
      any(diff(br) <= 0)) {
    stop("sev.age.breaks must be increasing lower limits of the severity bands, starting at 0 ",
         "and closed by Inf (GlobalMix India #71)", call. = FALSE)
  }
  nb <- length(br) - 1
  for (p in severity_params) {
    v <- get_param(dat, p, override.null.error = TRUE)
    if (is.null(v) || !is.numeric(v) || length(v) != nb || anyNA(v) || any(v < 0 | v > 1)) {
      stop("`", p, "` needs one probability per severity band (", nb, ", sev.age.breaks)",
           call. = FALSE)
    }
  }
  pc <- get_param(dat, "prop.clinical")
  for (p in c("ihr", "ifr")) {
    v <- get_param(dat, p)
    if (any(v > pc * (1 + 1e-9))) {
      stop("`", p, "` is a probability per infection and cannot exceed prop.clinical, the ",
           "symptomatic share; it does in severity band(s) ", toString(which(v > pc * (1 + 1e-9))),
           call. = FALSE)
    }
  }
  rr <- param_or(dat, "death.hosp.rr", 1)
  if (!(is.numeric(rr) && length(rr) == 1 && is.finite(rr) && rr >= 0)) {
    stop("death.hosp.rr must be one non-negative relative risk", call. = FALSE)
  }
  invisible(TRUE)
}

# Admission and death probabilities of new symptomatic cases, given whether
# they are admitted: P(death | symptomatic) is ifr / prop.clinical, split by
# admission with relative risk death.hosp.rr so that its mean is unchanged.
# p_hosp is the admission probability the case was drawn with.
death_given_admission <- function(p_death, p_hosp, hospit, rr) {
  if (rr == 1) return(pmin(1, p_death))
  base <- p_death / (p_hosp * rr + 1 - p_hosp)
  pmin(1, ifelse(hospit == 1, rr * base, base))
}

# Immunity after an infection --------------------------------------------------------

immunity_outcomes <- c("infect", "clinical", "hosp", "death")

# The names are written out so the parameter scan in test-gmc19-params.R sees
# every one. Unset, there is no partial protection and no return to s.
immunity_pars <- function(dat) {
  pr <- list(
    infect = c(peak = param_or(dat, "immune.peak.infect", 0),
               halflife = param_or(dat, "immune.halflife.infect", Inf),
               floor = param_or(dat, "immune.floor.infect", 0)),
    clinical = c(peak = param_or(dat, "immune.peak.clinical", 0),
                 halflife = param_or(dat, "immune.halflife.clinical", Inf),
                 floor = param_or(dat, "immune.floor.clinical", 0)),
    hosp = c(peak = param_or(dat, "immune.peak.hosp", 0),
             halflife = param_or(dat, "immune.halflife.hosp", Inf),
             floor = param_or(dat, "immune.floor.hosp", 0)),
    death = c(peak = param_or(dat, "immune.peak.death", 0),
              halflife = param_or(dat, "immune.halflife.death", Inf),
              floor = param_or(dat, "immune.floor.death", 0)))
  list(full = param_or(dat, "immune.full", Inf), protection = pr)
}

check_immunity_params <- function(dat) {
  ip <- immunity_pars(dat)
  if (!(is.numeric(ip$full) && length(ip$full) == 1 && !is.na(ip$full) && ip$full >= 1)) {
    stop("immune.full, the days of full protection after recovery, must be one value of at ",
         "least 1 (Inf: no reinfection in the run)", call. = FALSE)
  }
  pr <- ip$protection
  for (o in immunity_outcomes) {
    v <- pr[[o]]
    if (any(is.na(v)) || any(v[c("peak", "floor")] < 0 | v[c("peak", "floor")] > 1) ||
        v[["halflife"]] <= 0) {
      stop("immune.peak.", o, " and immune.floor.", o, " must be in [0, 1] and immune.halflife.",
           o, " positive", call. = FALSE)
    }
  }
  # Protection against a stage cannot be below protection against the stage
  # before it, or the conditional relative risk would exceed 1
  for (k in c("peak", "floor")) {
    p <- vapply(pr, function(v) v[[k]], numeric(1))
    if (p[["infect"]] > p[["clinical"]] || p[["clinical"]] > p[["hosp"]] ||
        p[["clinical"]] > p[["death"]]) {
      stop("post-infection protection must not fall along the stages (infect <= clinical <= ",
           "hosp, and clinical <= death); immune.", k, ".* are ", toString(signif(p, 3)),
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

# Marginal protection against outcome `o` at `days` since recovery; 0 for NA
# (no infection earlier in the run)
immune_protection <- function(ip, days, o) {
  v <- ip$protection[[o]]
  p <- v[["floor"]] + (v[["peak"]] - v[["floor"]]) * 0.5^(days / v[["halflife"]])
  p[is.na(days)] <- 0
  p
}

# Conditional relative risk of a stage given the stage before it
cond_rr <- function(p_out, p_prev) ifelse(p_prev < 1, pmin(1, (1 - p_out) / (1 - p_prev)), 0)

# Relative risk from immunity at stage `stage` for people `ids`: per contact for
# a susceptible ("infect", days since recovery at `at`), or for an episode's
# draws, at the days from the previous recovery to this infection
immune_rr <- function(dat, ids, stage, at = NULL) {
  recTime <- get_attr(dat, "recTime", override.null.error = TRUE)
  if (is.null(recTime) || length(ids) == 0) return(rep(1, length(ids)))
  ref <- if (stage == "infect") rep(at, length(ids)) else get_attr(dat, "infTime")[ids]
  days <- ref - recTime[ids]
  if (all(is.na(days))) return(rep(1, length(ids)))
  ip <- immunity_pars(dat)
  P <- function(o) immune_protection(ip, days, o)
  switch(stage,
         infect = 1 - P("infect"),
         clinical = cond_rr(P("clinical"), P("infect")),
         hosp = cond_rr(P("hosp"), P("clinical")),
         death = cond_rr(P("death"), P("clinical")),
         stop("unknown stage ", stage))
}

# Relative per-contact acquisition probability of susceptibles `ids`: the
# vaccine (with no infection VE for infants on the direct infant product,
# #41), the age band (susc.age.rr), young infants (susc.infant.rr) and immunity
# from an earlier infection (#69). Network contacts, the background hazard and
# imports all take it; imports used to take none of it (#70).
susceptibility_rr <- function(dat, ids, at, vax.schedule = build_vax_schedule(dat, get_pathogen(dat))) {
  if (length(ids) == 0) return(numeric(0))
  rr <- compute_ve(at = at, ids = ids, vax = get_attr(dat, "vax"),
                   last.dose.time = get_attr(dat, "last.dose.time"), vax.schedule = vax.schedule,
                   outcome = "infect", vax.age.group = vax_age_group_for(dat))$rr
  is_infant <- get_attr(dat, "is_infant", override.null.error = TRUE)
  if (isTRUE(as.logical(param_or(dat, "vax.infant.product", FALSE))) && !is.null(is_infant)) {
    rr[is_infant[ids]] <- 1
  }
  age <- get_attr(dat, "age")
  susc.age.rr <- get_param(dat, "susc.age.rr", override.null.error = TRUE)
  if (!is.null(susc.age.rr) && !all(is.na(susc.age.rr)) && any(susc.age.rr != 1, na.rm = TRUE)) {
    rs <- susc.age.rr[cut(age[ids], get_param(dat, "age.breaks"), labels = FALSE, right = FALSE)]
    rs[is.na(rs)] <- 1
    rr <- rr * rs
  }
  rr_young <- infant_susc_rr(dat, age[ids])
  if (!is.null(rr_young)) rr <- rr * rr_young
  rr * immune_rr(dat, ids, "infect", at)
}

# Episodes -------------------------------------------------------------------------

# Start an episode for `ids` at step `at`. Its severity draws are cleared so
# progress_general() draws them for this episode (#69), and inf.count counts the
# person's infections in the run.
start_episode <- function(dat, ids, at) {
  if (length(ids) == 0) return(dat)
  dat <- set_attr(dat, "status", "e", ids)
  dat <- set_attr(dat, "infTime", at, ids)
  dat <- set_attr(dat, "statusTime", at, ids)
  dat <- set_attr(dat, "clinical", NA, ids)
  dat <- set_attr(dat, "hospit", NA, ids)
  dat <- set_attr(dat, "fatal", NA, ids)
  n <- get_attr(dat, "inf.count")
  set_attr(dat, "inf.count", n[ids] + 1, ids)
}

# The introduction pulse (#70): import.pulse expected introductions over
# import.pulse.days days from day import.pulse.start, as a daily probability per
# susceptible that sums to import.pulse in a fully susceptible population
import_pulse_prob <- function(dat, at) {
  n <- param_or(dat, "import.pulse", 0)
  if (n <= 0) return(0)
  start <- param_or(dat, "import.pulse.start", 1)
  days <- param_or(dat, "import.pulse.days", 1)
  if (at < start || at >= start + days) return(0)
  n / days / sum(get_attr(dat, "active") == 1)
}

check_import_params <- function(dat) {
  vals <- c(import.rate = param_or(dat, "import.rate", 0), import.pulse = param_or(dat, "import.pulse", 0),
            import.pulse.start = param_or(dat, "import.pulse.start", 1),
            import.pulse.days = param_or(dat, "import.pulse.days", 1))
  if (any(!is.finite(vals)) || any(vals < 0) || vals[["import.rate"]] > 1 ||
      vals[["import.pulse.days"]] < 1 || vals[["import.pulse.start"]] < 1) {
    stop("import.rate must be a daily probability, import.pulse a non-negative count, and ",
         "import.pulse.start and import.pulse.days at least 1", call. = FALSE)
  }
  invisible(TRUE)
}

# Outcomes by severity band (#71 acceptance), with sev.diagnostics = TRUE:
# per step and band, the events `what` (inf, sympt, hosp, death) among `ids`,
# all episodes and first episodes. With the infections, the living and the
# prior-immune by band, so realized IFR and IHR per infection and the starting
# immunity can be read by band.
sev_diagnostics <- function(dat) isTRUE(as.logical(param_or(dat, "sev.diagnostics", FALSE)))

sev_flows <- function(dat, at, what, ids, first) {
  nb <- length(get_param(dat, "sev.age.breaks")) - 1
  b <- sev_band(dat, get_attr(dat, "age")[ids])
  n_all <- tabulate(b, nb)
  n_first <- tabulate(b[first], nb)
  for (k in seq_len(nb)) {
    dat <- set_epi(dat, paste0("sev.", what, ".", k), n_all[k], at = at)
    dat <- set_epi(dat, paste0("sev.", what, ".first.", k), n_first[k], at = at)
  }
  if (what == "inf") {
    live <- get_attr(dat, "active") == 1
    bl <- sev_band(dat, get_attr(dat, "age")[live])
    prior <- get_attr(dat, "prior.immune")[live] == 1
    n_live <- tabulate(bl, nb)
    n_prior <- tabulate(bl[prior], nb)
    for (k in seq_len(nb)) {
      dat <- set_epi(dat, paste0("sev.n.", k), n_live[k], at = at)
      dat <- set_epi(dat, paste0("sev.prior.", k), n_prior[k], at = at)
    }
  }
  dat
}

# Prior immunity and the seeds at initialization ----------------------------------------

check_prior_immunity <- function(dat) {
  f <- param_or(dat, "init.immune", 0)
  nb <- length(get_param(dat, "sev.age.breaks")) - 1
  if (length(f) == 1 && f == 0) return(invisible(TRUE))
  if (!is.numeric(f) || length(f) != nb || anyNA(f) || any(f < 0 | f > 1)) {
    stop("init.immune needs one share in [0, 1] per severity band (", nb, ")", call. = FALSE)
  }
  attr <- immune_weight_attr(dat)
  if (!is.null(attr) && !attr %in% names(get_attr_list(dat))) {
    stop("init.immune.attr names `", attr, "`, which the population does not carry as a ",
         "vertex attribute", call. = FALSE)
  }
  invisible(TRUE)
}

immune_weight_attr <- function(dat) {
  a <- param_or(dat, "init.immune.attr", "")
  if (is.null(a) || is.na(a) || identical(as.character(a), "")) NULL else as.character(a)
}

# Inclusion probabilities min(1, c w) with mean `share`, c found by bisection
weighted_inclusion <- function(w, share) {
  n <- length(w)
  if (share <= 0 || n == 0) return(rep(0, n))
  if (share >= 1) return(rep(1, n))
  w[!is.finite(w) | w < 0] <- 0
  if (all(w == 0)) w[] <- 1
  if (mean(w > 0) < share) {
    stop("init.immune asks for ", signif(share, 3), " of a band immune, but only ",
         signif(mean(w > 0), 3), " of it has a positive propensity", call. = FALSE)
  }
  lo <- 0; hi <- share / mean(w[w > 0]) * 1e6
  for (i in 1:200) {
    mid <- (lo + hi) / 2
    if (mean(pmin(1, mid * w)) < share) lo <- mid else hi <- mid
  }
  pmin(1, hi * w)
}

# The run's starting state of infection and immunity: the seeds (status e from
# init.net(e.num)) are infections in the run, and init.immune makes a share of
# each severity band immune for the run
init_immunity <- function(dat) {
  status <- get_attr(dat, "status")
  active <- get_attr(dat, "active")
  seeds <- which(active == 1 & status == "e")
  if (length(seeds) > 0) {
    dat <- set_attr(dat, "infTime", 1, seeds)
    dat <- set_attr(dat, "inf.count", 1, seeds)
    dat <- set_attr(dat, "ever.inf", 1, seeds)
  }
  f <- param_or(dat, "init.immune", 0)
  if (all(f == 0)) return(dat)
  band <- sev_band(dat, get_attr(dat, "age"))
  attr <- immune_weight_attr(dat)
  w_all <- if (is.null(attr)) rep(1, length(active)) else as.numeric(get_attr(dat, attr))
  prior <- get_attr(dat, "prior.immune")
  for (b in which(f > 0)) {
    # Seeds are infected today, so they are not immune: the band's share is
    # drawn over everyone else in it, scaled to hold the share of the band
    ids_all <- which(active == 1 & band == b)
    ids <- ids_all[status[ids_all] == "s"]
    if (length(ids) == 0) next
    p <- weighted_inclusion(w_all[ids], min(1, f[b] * length(ids_all) / length(ids)))
    imm <- ids[stats::runif(length(ids)) < p]
    prior[imm] <- 1
  }
  imm <- which(prior == 1)
  dat <- set_attr(dat, "prior.immune", prior)
  dat <- set_attr(dat, "status", "r", imm)
  dat <- set_attr(dat, "statusTime", 1, imm)
  dat
}
