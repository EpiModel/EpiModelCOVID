# Relationships, daily activation and new contacts (GlobalMix India #63, #68)
#
# Out-of-home contact is two relationship networks plus new contacts. The
# networks hold who knows whom; this module decides each day which of those
# relationships are met:
#   - each person is out with probability o_i, a persistent trait drawn from
#     Beta(o kappa, (1 - o) kappa) by age band; never-active people are never out
#   - a frequent relationship (seen daily or 1-3 times a week) is met on a day
#     both people are out, with probability q
#   - an infrequent relationship (every 2 weeks or less often) is met on a day
#     both people are out and the tie is in an episode, a two-state Markov chain
#     per tie with persistence p11 and stationary probability s, carried from
#     day to day
#   - new contacts (never met before) are drawn fresh each day by age pair. In
#     the diary half of them fall on days with no relationship contact (51%
#     rural, 48% urban, against 30% and 10% of person-days), and they are not
#     associated with frequent activity, so each day a fixed share of new-contact
#     ends goes to the people with no relationship met that day, the never-active
#     included, and the rest to everyone else: a person with no relationship met
#     carries an odds multiplier, solved each day, over one with one
# The day's contacts by tie type go to the infection module in
# dat$run$contacts. The household layer is met every day.
#
# The parameters come from the GlobalMix diary
# (GlobalMix-India-NetVax R/A-networks/relationships/estimate_relationships.R);
# the methods are written up in that repository's
# plan/notes/relationship_activation_methods.md.

# Parameters ------------------------------------------------------------------

# The model runs the activation step when its parameters are set.
activation_on <- function(dat) {
  !is.null(get_param(dat, "act.o", override.null.error = TRUE))
}

# Out-of-home tie types and the names of the two relationship networks.
rel_layers <- c(freq = "freq", infreq = "infreq")
new_type <- "new"

param_or <- function(dat, name, default) {
  value <- get_param(dat, name, override.null.error = TRUE)
  if (is.null(value)) default else value
}

#' Activation parameters with the scenario multipliers applied
#'
#' @param dat Main `netsim_dat` object.
#' @return A list: o, q, pi, theta_F, new_rate (one value per age band), kappa,
#'   p11, s, theta_I, rho, new_mix (age band by age band), episodes (logical),
#'   window (days) and class_mult (the sex and occupation class multipliers of
#'   the activities, a table with one row per class, or NULL).
#' @keywords internal
activation_params <- function(dat) {
  nb <- length(get_param(dat, "age.grps"))
  p <- list(
    o = get_param(dat, "act.o"),
    kappa = get_param(dat, "act.kappa") * param_or(dat, "act.kappa.mult", 1),
    q = pmin(get_param(dat, "act.q") * param_or(dat, "act.q.mult", 1), 1),
    p11 = get_param(dat, "act.p11"),
    s = get_param(dat, "act.s"),
    pi = get_param(dat, "act.pi"),
    theta_F = get_param(dat, "act.theta.F"),
    theta_I = get_param(dat, "act.theta.I"),
    rho = get_param(dat, "act.rho"),
    new_rate = get_param(dat, "new.rate") * param_or(dat, "new.rate.mult", 1),
    new_mix = get_param(dat, "new.mix"),
    new_share = get_param(dat, "new.share.norel"),
    theta_N = param_or(dat, "act.theta.N", Inf),
    pi_N = param_or(dat, "act.pi.N", 0),
    episodes = isTRUE(as.logical(param_or(dat, "act.episodes", TRUE))),
    window = as.integer(param_or(dat, "act.window", 30)),
    class_mult = normalize_class_mult(param_or(dat, "act.class.mult", NULL), get_param(dat, "age.grps"))
  )
  bad <- function(what) stop("activation parameter ", what, call. = FALSE)
  for (v in c("o", "q", "pi", "theta_F", "new_rate")) {
    if (length(p[[v]]) != nb || anyNA(p[[v]])) bad(paste0("`", v, "` needs one value per age band (", nb, "), none NA"))
  }
  in01 <- function(x) all(x >= 0 & x <= 1)
  if (!in01(p$o) || !in01(p$q) || !in01(p$pi)) bad("o, q or pi outside [0, 1]")
  if (any(p$theta_F <= 0) || length(p$theta_I) != 1 || !isTRUE(p$theta_I > 0)) bad("theta must be positive")
  if (length(p$kappa) != 1 || !isTRUE(p$kappa > 0)) bad("kappa (act.kappa x act.kappa.mult) must be positive")
  if (length(p$rho) != 1 || !isTRUE(abs(p$rho) <= 1)) bad("rho must be in [-1, 1]")
  if (any(p$new_rate < 0)) bad("new.rate must be non-negative")
  if (length(p$new_share) != 1 || !isTRUE(p$new_share >= 0 && p$new_share <= 1)) {
    bad("`new.share.norel` must be one share in [0, 1]")
  }
  if (!is.matrix(p$new_mix) || any(dim(p$new_mix) != nb) || anyNA(p$new_mix) || any(p$new_mix < 0) ||
      !isSymmetric(unname(p$new_mix))) {
    bad(paste0("`new.mix` must be a symmetric non-negative ", nb, " x ", nb, " matrix"))
  }
  if (length(p$window) != 1 || !isTRUE(p$window >= 1)) bad("`act.window` must be at least 1")
  if (length(p$theta_N) != 1 || !isTRUE(p$theta_N > 0)) bad("`act.theta.N` must be positive")
  if (length(p$pi_N) != 1 || !isTRUE(p$pi_N >= 0 && p$pi_N < 1)) bad("`act.pi.N` must be in [0, 1)")
  # The off-to-on probability s (1 - p11) / (1 - s) must be a probability
  if (!isTRUE(p$s > 0 && p$p11 >= 0 && p$p11 < 1 && p$s <= 1 / (2 - p$p11))) {
    bad(paste0("episode parameters out of range: act.s ", p$s, ", act.p11 ", p$p11))
  }
  p
}

# The sex and occupation class multipliers of the frequent and infrequent
# activities (GlobalMix India #61), act.class.mult: a table with one row per
# class, band (the age band's index in age.grps, or its label), sex, occ, mult_F
# and mult_I. Returned with integer bands, sorted by class; NULL when unset, in
# which case every multiplier is 1.
normalize_class_mult <- function(cm, bands) {
  if (is.null(cm)) return(NULL)
  bad <- function(what) stop("act.class.mult ", what, call. = FALSE)
  if (!is.data.frame(cm) || !all(c("band", "sex", "occ", "mult_F", "mult_I") %in% names(cm))) {
    bad("must be a data frame with columns band, sex, occ, mult_F and mult_I")
  }
  band <- if (is.numeric(cm$band)) cm$band else match(as.character(cm$band), bands)
  if (anyNA(band) || any(band != round(band)) || any(band < 1 | band > length(bands))) {
    bad(paste0("needs each band as its index in age.grps (1 to ", length(bands), ")"))
  }
  out <- data.frame(band = as.integer(band), sex = as.character(cm$sex), occ = as.character(cm$occ),
                    mult_F = as.numeric(cm$mult_F), mult_I = as.numeric(cm$mult_I))
  if (!all(out$sex %in% sex_levels) || !all(out$occ %in% occ_levels)) {
    bad(paste0("has a sex other than ", toString(sex_levels), " or an occ other than ",
               toString(occ_levels)))
  }
  m <- c(out$mult_F, out$mult_I)
  if (!all(is.finite(m)) || any(m <= 0)) bad("multipliers must be positive")
  if (anyDuplicated(out[c("band", "sex", "occ")]) > 0) bad("has more than one row for a class")
  out <- out[order(out$band, out$sex, out$occ), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Each node's activity multipliers: the row of its class (age band, sex,
# occupation class) in the class multiplier table
class_multipliers <- function(band, sex, occ, cm) {
  n <- length(band)
  if (length(sex) != n || length(occ) != n) {
    stop("the class multipliers (act.class.mult) need each node's sex and occ", call. = FALSE)
  }
  if (!is.numeric(cm$band)) {
    stop("the class multipliers need each band as its index in age.grps", call. = FALSE)
  }
  key <- paste(band, sex, occ, sep = "|")
  row <- match(key, paste(cm$band, cm$sex, cm$occ, sep = "|"))
  if (anyNA(row)) {
    miss <- unique(data.frame(band = band, sex = sex, occ = occ)[is.na(row), , drop = FALSE])
    stop("act.class.mult has no row for ",
         paste0("band ", miss$band, ", sex ", miss$sex, ", occ ", miss$occ, collapse = "; "),
         call. = FALSE)
  }
  list(mult_F = cm$mult_F[row], mult_I = cm$mult_I[row])
}

# Newborns draw their traits with their own class: occ "child", in the band of
# arrival.age, either sex. Checked at initialization rather than at the first
# birth.
check_newborn_classes <- function(dat) {
  cm <- activation_params(dat)$class_mult
  if (is.null(cm)) return(invisible(TRUE))
  band <- cut(get_param(dat, "arrival.age"), get_param(dat, "age.breaks"), labels = FALSE, right = FALSE)
  class_multipliers(rep(band, 2), sex_levels, rep("child", 2), cm)
  invisible(TRUE)
}

# The out-propensity and trait parameters act only when a person's traits are
# drawn (at initialization, and for newborns), and the initial population's
# relationship traits are those the networks were fitted on. A parameter
# updater that changed them mid-run would reach only newborns, and a scenario
# that changed the trait parameters would make newborns differ from the initial
# population, so both stop the run.
check_trait_params <- function(dat) {
  at_start <- c("act.o", "act.kappa", "act.kappa.mult", "act.theta.N", "act.pi.N")
  traits <- c("act.pi", "act.theta.F", "act.theta.I", "act.rho", "act.class.mult")
  late <- unlist(lapply(dat$param$.param.updater.list, function(u) if (isTRUE(u$at > 1)) names(u$param)))
  hit <- intersect(late, c(at_start, traits))
  if (length(hit)) {
    stop("a parameter updater after the first step changes ", toString(hit), ", which act only when ",
         "traits are drawn; set them in the scenario (at = 1)", call. = FALSE)
  }
  fitted <- get_nwparam(dat, network = match(rel_layers[["freq"]], names(dat$run$el)))$trait_pars
  if (!is.null(fitted)) {
    now <- list(pi = get_param(dat, "act.pi"), theta_F = get_param(dat, "act.theta.F"),
                theta_I = get_param(dat, "act.theta.I"), rho = get_param(dat, "act.rho"))
    diff <- names(now)[!vapply(names(now), function(k) isTRUE(all.equal(unname(now[[k]]),
                                                                        unname(fitted[[k]]))), logical(1))]
    # The sex and occupation class multipliers (GlobalMix India #61): none on
    # both, or the same table
    bands <- get_param(dat, "age.grps")
    if (!isTRUE(all.equal(normalize_class_mult(param_or(dat, "act.class.mult", NULL), bands),
                          normalize_class_mult(fitted$class_mult, bands)))) {
      diff <- c(diff, "class_mult")
    }
    if (length(diff)) {
      stop("the trait parameters (", toString(diff), ") differ from those the networks were fitted ",
           "with; newborns would differ from the initial population. Refit the networks instead.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

# Age band index of every node (1..6), from age.grp.
band_index <- function(dat) {
  match(get_attr(dat, "age.grp"), get_param(dat, "age.grps"))
}

# Traits -----------------------------------------------------------------------

#' Draw the relationship traits of a population
#'
#' Never-active flag and the frequent and infrequent activities. The activities
#' are mean-1 gamma variates joined by a Gaussian copula with correlation `rho`:
#' z1 ~ N(0, 1), z2 = rho z1 + sqrt(1 - rho^2) e, a_F = G_F^-1(Phi(z1)) with
#' shape theta_F by age band, a_I = G_I^-1(Phi(z2)) with shape theta_I. The fitted
#' networks carry them as vertex attributes (`never`, `log_aF`, `log_aI`), and
#' the arrivals module draws them for newborns from the youngest band.
#'
#' With class multipliers (`pars$class_mult`, GlobalMix India #61), each
#' activity is the gamma variate times its class's multiplier, a_F = mult_F
#' G_F^-1(Phi(z1)) and a_I = mult_I G_I^-1(Phi(z2)): the class sets the mean
#' within the age band and the shape is unchanged. The never-active share stays
#' by band. Without them, the draws are those of EpiModelCOVID 1.7.2.
#'
#' @param band Age band index of each node (1..6).
#' @param pars List with `pi` and `theta_F` (one value per band), `theta_I` and
#'   `rho`, and optionally `class_mult`: a data frame with one row per class,
#'   `band` (the band index), `sex`, `occ`, `mult_F` and `mult_I`.
#' @param seed Optional seed.
#' @param sex,occ Sex ("F", "M") and occupation class ("child", "student",
#'   "worker", "home") of each node; required with `pars$class_mult`.
#' @return A list with `never` (0/1), `log_aF` and `log_aI`.
#' @export
draw_relationship_traits <- function(band, pars, seed = NULL, sex = NULL, occ = NULL) {
  mult <- if (!is.null(pars$class_mult)) class_multipliers(band, sex, occ, pars$class_mult)
  if (!is.null(seed)) set.seed(seed)
  n <- length(band)
  never <- as.integer(stats::runif(n) < pars$pi[band])
  z1 <- stats::rnorm(n)
  z2 <- pars$rho * z1 + sqrt(1 - pars$rho^2) * stats::rnorm(n)
  thF <- pars$theta_F[band]
  aF <- stats::qgamma(stats::pnorm(z1), shape = thF, rate = thF)
  aI <- stats::qgamma(stats::pnorm(z2), shape = pars$theta_I, rate = pars$theta_I)
  if (!is.null(mult)) {
    aF <- aF * mult$mult_F
    aI <- aI * mult$mult_I
  }
  list(never = never, log_aF = log(pmax(aF, 1e-8)), log_aI = log(pmax(aI, 1e-8)))
}

#' Draw each person's probability of going out on a day
#'
#' o_i ~ Beta(o kappa, (1 - o) kappa) with mean o by age band. A concentration
#' of 1e4 or more is treated as no variation.
#'
#' @param band Age band index of each node.
#' @param o Mean probability of going out, one value per band.
#' @param kappa Beta concentration.
#' @return A numeric vector.
#' @export
draw_out_propensity <- function(band, o, kappa) {
  ob <- o[band]
  if (!is.finite(kappa) || kappa >= 1e4) return(ob)
  stats::rbeta(length(band), ob * kappa, (1 - ob) * kappa)
}

#' Map relationship traits to a new age band, keeping each person's rank
#'
#' When a person crosses a band top (GlobalMix India #64), their traits move to
#' the new band's distributions with their rank kept: each activity by the
#' quantile of its gamma variate, under the new band's gamma shape and class
#' multiplier; the propensity to go out by its Beta quantile; and the
#' never-active flag re-drawn only as far as the new band's share needs, so a
#' person who is active stays active where the share falls. A person who stays
#' in their band keeps their traits exactly.
#'
#' @param from,to Age band index of each node before and after the crossing.
#' @param traits List with `never`, `log_aF`, `log_aI` and `o_i`, one value per
#'   node.
#' @param pars Activation parameters, as from the parameters of the run: `pi`,
#'   `theta_F`, `theta_I`, `o`, `kappa` and `class_mult`.
#' @param sex,occ Each node's sex and occupation class, when `pars$class_mult`
#'   is set.
#' @return `traits` with the mapped values.
#' @export
remap_relationship_traits <- function(from, to, traits, pars, sex = NULL, occ = NULL) {
  n <- length(from)
  if (n == 0) return(traits)
  if (length(to) != n || any(vapply(traits[c("never", "log_aF", "log_aI", "o_i")], length, 1L) != n)) {
    stop("remap_relationship_traits() needs one band and one value of each trait per node", call. = FALSE)
  }
  one <- list(mult_F = rep(1, n), mult_I = rep(1, n))
  m0 <- if (!is.null(pars$class_mult)) class_multipliers(from, sex, occ, pars$class_mult) else one
  m1 <- if (!is.null(pars$class_mult)) class_multipliers(to, sex, occ, pars$class_mult) else one
  # a = m G^-1(u; theta, theta), so u = G(a / m) and the new a is m' G^-1(u; theta')
  qmap <- function(la, mult0, mult1, th0, th1) {
    move <- mult0 != mult1 | th0 != th1
    if (!any(move)) return(la)
    u <- stats::pgamma(exp(la[move]) / mult0[move], shape = th0[move], rate = th0[move])
    a <- mult1[move] * stats::qgamma(pmin(u, 1 - 1e-15), shape = th1[move], rate = th1[move])
    la[move] <- log(pmax(a, 1e-8))
    la
  }
  traits$log_aF <- qmap(traits$log_aF, m0$mult_F, m1$mult_F, pars$theta_F[from], pars$theta_F[to])
  traits$log_aI <- qmap(traits$log_aI, m0$mult_I, m1$mult_I, rep(pars$theta_I, n), rep(pars$theta_I, n))

  # never: from share p0 to p1, flipping as few people as the new share needs
  p0 <- pars$pi[from]; p1 <- pars$pi[to]
  u <- stats::runif(n)
  nv <- traits$never
  on <- nv == 0 & p1 > p0 & u < (p1 - p0) / pmax(1 - p0, 1e-12)
  off <- nv == 1 & p1 < p0 & u < (p0 - p1) / pmax(p0, 1e-12)
  nv[on] <- 1L
  nv[off] <- 0L
  traits$never <- nv

  # going out: o_i ~ Beta(o kappa, (1 - o) kappa), mapped by its quantile; with
  # no variation (draw_out_propensity()'s rule) or a degenerate mean, the band mean
  o0 <- pars$o[from]; o1 <- pars$o[to]
  k <- pars$kappa
  move <- o0 != o1
  if (any(move)) {
    oi <- traits$o_i
    if (!is.finite(k) || k >= 1e4) {
      oi[move] <- o1[move]
    } else {
      deg <- move & (o0 <= 0 | o0 >= 1 | o1 <= 0 | o1 >= 1)
      q <- move & !deg
      oi[deg] <- o1[deg]
      if (any(q)) {
        uq <- stats::pbeta(oi[q], o0[q] * k, (1 - o0[q]) * k)
        oi[q] <- stats::qbeta(uq, o1[q] * k, (1 - o1[q]) * k)
      }
    }
    traits$o_i <- oi
  }
  traits
}

# Band tops (GlobalMix India #64) -------------------------------------------------

#' People who cross a band top: traits mapped to the new band
#'
#' Called by [aging()] on the relationship design with the nodes whose age band
#' changed on this step. Ages are continuous, so people cross band tops during a
#' run, but each age pair's relationships are fit as if its bands were closed:
#' a person who kept their ties and traits would carry the old band's
#' relationships for years, and over a season band degree drifted by up to 8%
#' and adjacent-band mixing by up to 76%. So a crosser becomes exchangeable with
#' the band they enter, as a newborn is with the initial infants: their traits
#' are mapped to the new band here ([remap_relationship_traits()]), and their
#' relationships are redrawn from the new band's cross-section after the day's
#' network step ([seed_relationship_ties()]). `rel.cross` records the step.
#' Off when `act.reset.cross` is FALSE, which keeps the traits and ties of the
#' engine before 1.9.0.
#'
#' @inheritParams aging
#' @param ids Positional ids of the crossers.
#' @param from Their age band index before the crossing.
#' @return The updated `dat`.
#' @keywords internal
cross_band_tops <- function(dat, at, ids, from) {
  if (!length(ids) || isFALSE(as.logical(param_or(dat, "act.reset.cross", TRUE)))) return(dat)
  pars <- activation_params(dat)
  to <- band_index(dat)[ids]
  cls <- has_pop_class(dat)
  traits <- list(never = get_attr(dat, "never")[ids], log_aF = get_attr(dat, "log_aF")[ids],
                 log_aI = get_attr(dat, "log_aI")[ids], o_i = get_attr(dat, "o_i")[ids])
  tr <- remap_relationship_traits(from, to, traits, pars,
                                  sex = if (cls) get_attr(dat, "sex")[ids],
                                  occ = if (cls) get_attr(dat, "occ")[ids])
  for (a in names(tr)) dat <- set_attr(dat, a, tr[[a]], posit_ids = ids)
  dat <- set_attr(dat, "act.eligible", 1 - tr$never, posit_ids = ids)
  dat <- set_attr(dat, "rel.cross", at, posit_ids = ids)
  dat
}

# The day's contacts --------------------------------------------------------------

# Key of an unordered pair of unique ids. Unique ids are never reused, and the
# key is exact in double precision for ids below 2^26.
pair_key <- function(u1, u2) {
  pmin(u1, u2) * 67108864 + pmax(u1, u2)
}

#' Advance the episode state of the infrequent ties by one day
#'
#' @param key Pair keys of today's infrequent ties.
#' @param prev List with `key` and `on` from the previous day, or NULL.
#' @param p11 Probability an episode continues to the next day.
#' @param s Stationary probability a tie is in an episode.
#' @param episodes If FALSE, ties are in an episode on independent days with
#'   probability s.
#' @return Logical state for each tie.
#' @keywords internal
episode_update <- function(key, prev, p11, s, episodes = TRUE) {
  u <- stats::runif(length(key))
  if (!episodes || is.null(prev)) return(u < s)
  p01 <- s * (1 - p11) / (1 - s)
  was <- prev$on[match(key, prev$key)]
  ifelse(is.na(was), u < s, ifelse(was, u < p11, u < p01))
}

# Symmetric iterative proportional fitting: scale S to diag(f) S diag(f) so the
# edge ends of each band equal `target`, where a diagonal cell holds within-band
# edges (each with two ends). The same fit as the GlobalMix target pipeline's
# sym_ipf().
sym_ipf_edges <- function(S, target, tol = 1e-10, maxit = 10000) {
  pos <- target > 0 & rowSums(S) > 0
  if (!any(pos)) return(S * 0)
  f <- as.numeric(pos)
  for (it in seq_len(maxit)) {
    E <- S * outer(f, f)
    ends <- rowSums(E) + diag(E)
    r <- ifelse(pos, target / pmax(ends, .Machine$double.xmin), 1)
    if (max(abs(r[pos] - 1)) < tol) break
    f <- f * sqrt(r)
  }
  S * outer(f, f)
}

#' Expected new-contact edges by age pair for the current population
#'
#' @param new_mix Per-capita age mixing shape of out-of-home contact (symmetric,
#'   diagonal = within-band edges).
#' @param new_rate New contacts per person per day, by band.
#' @param N_b Current number of people in each band.
#' @return Symmetric matrix of expected edges (diagonal = within-band edges).
#' @export
new_contact_edges <- function(new_mix, new_rate, N_b) {
  sym_ipf_edges(new_mix * sum(N_b), new_rate * N_b)
}

#' Draw the day's new contacts
#'
#' For each age pair, a Poisson number of edges with the expected count in `E`,
#' each end drawn from the band's people with probability proportional to
#' `weight`. Pairs of a person with themself are dropped.
#'
#' @param E Expected edges by age pair (from [new_contact_edges()]).
#' @param band Age band index of each node.
#' @param weight Draw weight of each node (0: never drawn).
#' @return Two-column integer matrix of positional ids, tail < head.
#' @export
draw_new_contacts <- function(E, band, weight) {
  nb <- nrow(E)
  idx <- split(seq_along(band), factor(band, levels = seq_len(nb)))
  ta <- vector("list", nb * nb); tb <- ta; k <- 0
  for (i in seq_len(nb)) for (j in i:nb) {
    m <- stats::rpois(1, E[i, j])
    if (m == 0) next
    pi_ <- idx[[i]]; pj <- idx[[j]]
    wi <- weight[pi_]; wj <- weight[pj]
    if (!length(pi_) || !length(pj) || sum(wi) <= 0 || sum(wj) <= 0) next
    a <- pi_[sample.int(length(pi_), m, replace = TRUE, prob = wi)]
    b <- pj[sample.int(length(pj), m, replace = TRUE, prob = wj)]
    k <- k + 1
    ta[[k]] <- a; tb[[k]] <- b
  }
  a <- unlist(ta[seq_len(k)]); b <- unlist(tb[seq_len(k)])
  keep <- a != b
  el <- cbind(pmin(a, b)[keep], pmax(a, b)[keep])
  storage.mode(el) <- "integer"
  el
}

# Odds multiplier on people with no relationship met today such that a share
# `share` of all new-contact ends goes to them, given each band's expected ends
# and its counts of people with none (n0) and with some (n1).
new_contact_odds <- function(ends, n0, n1, share) {
  if (sum(ends) <= 0 || sum(n0) == 0 || sum(n1) == 0) return(1)
  f <- function(ll) {
    l <- exp(ll)
    d <- l * n0 + n1
    sum(ifelse(d > 0, ends * l * n0 / d, 0)) / sum(ends) - share
  }
  lo <- f(-30); hi <- f(30)
  if (lo >= 0) return(exp(-30))
  if (hi <= 0) return(exp(30))
  exp(stats::uniroot(f, c(-30, 30), f.lower = lo, f.upper = hi, tol = 1e-8)$root)
}

#' The day's out-of-home contacts
#'
#' The engine's daily activation, as a function of its inputs so that it can be
#' run outside netsim (the GlobalMix India calibration harness calls it).
#'
#' @param el_freq,el_infreq Edgelists of the two relationship networks
#'   (positional ids).
#' @param uid Unique id of every node (for the episode keys).
#' @param band Age band index of every node.
#' @param out Logical: who is out today.
#' @param live Logical: who is alive (eligible for new contacts).
#' @param pars Output of [activation_params()].
#' @param aN Each node's propensity for new contacts (mean 1), or NULL for none.
#' @param episode Previous day's episode state (list with key and on), or NULL.
#' @param N_b Number of people in each band (for the new-contact rate).
#' @return A list: freq, infreq and new edgelists of the contacts met today,
#'   and the infrequent ties' episode state to carry to tomorrow.
#' @export
activate_day <- function(el_freq, el_infreq, uid, band, out, live, pars, episode, N_b, aN = NULL) {
  qa <- pars$q[band[el_freq[, 1]]]; qb <- pars$q[band[el_freq[, 2]]]
  onF <- out[el_freq[, 1]] & out[el_freq[, 2]] & stats::runif(nrow(el_freq)) < (qa + qb) / 2
  keyI <- pair_key(uid[el_infreq[, 1]], uid[el_infreq[, 2]])
  st <- episode_update(keyI, episode, pars$p11, pars$s, pars$episodes)
  onI <- out[el_infreq[, 1]] & out[el_infreq[, 2]] & st
  freq <- el_freq[onF, , drop = FALSE]; infreq <- el_infreq[onI, , drop = FALSE]
  ## New contacts: each person's draw weight is their persistent new-contact
  ## propensity aN. Each band's ends are fixed by the age mixing, so the share of
  ## all ends on people with no relationship met today is set by an odds
  ## multiplier lambda on those people, solved for the day's band composition
  none <- live & tabulate(c(freq, infreq), nbins = length(live)) == 0
  E <- new_contact_edges(pars$new_mix, pars$new_rate, N_b)
  w <- as.numeric(live) * (if (is.null(aN)) 1 else aN)
  wsum <- function(sel) as.numeric(tapply(w[sel], factor(band[sel], levels = seq_along(N_b)), sum))
  n0 <- wsum(none); n1 <- wsum(live & !none)
  n0[is.na(n0)] <- 0; n1[is.na(n1)] <- 0
  lambda <- new_contact_odds(rowSums(E) + diag(E), n0, n1, pars$new_share)
  w[none] <- w[none] * lambda
  el_new <- draw_new_contacts(E, band, w)
  list(freq = freq, infreq = infreq, new = el_new, episode = list(key = keyI, on = st))
}

# Module ------------------------------------------------------------------------

#' @title Daily Activation of Relationships and New Contacts
#'
#' @description Draws who is out, which relationships are met, and the day's new
#'   contacts, and records the day's contacts by tie type in
#'   `dat$run$contacts` for the infection module. Also keeps the per-person
#'   contact counters the targeting scores read, and the diary-comparable
#'   diagnostics. Does nothing when the activation parameters are not set.
#'
#' @inheritParams aging
#' @return The updated `dat`.
#' @rdname moduleset-gmc19
#' @export
activate_contacts <- function(dat, at) {
  if (!activation_on(dat)) return(dat)
  pars <- activation_params(dat)
  active <- get_attr(dat, "active")
  never <- get_attr(dat, "never")
  o_i <- get_attr(dat, "o_i")
  band <- band_index(dat)
  uid <- get_unique_ids(dat)
  n <- length(active)
  nb <- length(pars$o)
  live <- active == 1

  dat <- seed_relationship_ties(dat, at)

  # A daily activity schedule (#68) scales the propensity to go out: a vector of
  # multipliers recycled over the run (7 values for a weekly pattern). Absent or
  # all ones, it leaves the draw unchanged.
  sched <- get_param(dat, "act.schedule", override.null.error = TRUE)
  if (!is.null(sched) && any(sched != 1)) o_i <- pmin(o_i * sched[(at - 1) %% length(sched) + 1], 1)
  out <- live & never == 0 & stats::runif(n) < o_i
  N_b <- tabulate(band[live], nbins = nb)
  el <- dat$run$el
  day <- activate_day(el[[rel_layers[["freq"]]]], el[[rel_layers[["infreq"]]]], uid, band, out, live,
                      pars, dat$run$episode, N_b, aN = exp(get_attr(dat, "log_aN")))
  dat$run$episode <- day$episode
  hh <- hh_layer(dat)
  contacts <- list(day$freq, day$infreq, day$new)
  names(contacts) <- c(rel_layers[["freq"]], rel_layers[["infreq"]], new_type)
  if (!is.na(hh)) contacts[[names(el)[hh]]] <- el[[hh]]
  dat$run$contacts <- contacts[tie_types(dat)]
  dat$run$contacts_at <- at

  ## Per-person counts of the day's out-of-home contacts
  degF <- tabulate(day$freq, nbins = n)
  degI <- tabulate(day$infreq, nbins = n)
  degN <- tabulate(day$new, nbins = n)
  deg <- degF + degI + degN
  dat <- set_attr(dat, "act.out.days", get_attr(dat, "act.out.days") + out)
  dat <- set_attr(dat, "act.contact.days", get_attr(dat, "act.contact.days") + deg)
  dat <- set_attr(dat, "act.new.contacts", get_attr(dat, "act.new.contacts") + degN)
  dat <- update_recent_contacts(dat, uid, deg, degN, pars$window)

  if (!isFALSE(as.logical(param_or(dat, "act.diagnostics", TRUE)))) {
    dat <- activation_diagnostics(dat, at, day, uid, out, live, degF, degI, degN, band)
  }
  dat
}

# Contacts over the last `window` days, from a ring buffer of daily counts
# aligned to the nodes by unique id (dead nodes drop out, newborns enter at 0).
update_recent_contacts <- function(dat, uid, deg, degN, window) {
  ring <- dat$run$act_ring
  if (is.null(ring) || ncol(ring$all) != window) {
    ring <- list(uid = uid, all = matrix(0L, length(uid), window), new = matrix(0L, length(uid), window), pos = 0L)
  }
  if (!identical(ring$uid, uid)) {
    idx <- match(uid, ring$uid)
    fill <- function(m) {
      out <- m[idx, , drop = FALSE]
      out[is.na(idx), ] <- 0L
      out
    }
    ring$all <- fill(ring$all); ring$new <- fill(ring$new); ring$uid <- uid
  }
  ring$pos <- ring$pos %% window + 1L
  ring$all[, ring$pos] <- as.integer(deg)
  ring$new[, ring$pos] <- as.integer(degN)
  dat$run$act_ring <- ring
  dat <- set_attr(dat, "act.recent", rowSums(ring$all))
  dat <- set_attr(dat, "act.recent.new", rowSums(ring$new))
  dat
}

# Diary-comparable diagnostics, emitted as counts so that they can be summed
# over days and simulations before dividing:
#   act.out           share of living people out
#   act.deg.freq/infreq/new  mean contacts met per living person
#   act.iso.rel/all   share with no relationship contact / no out-of-home contact
#   act.new.norel, act.new.ends  new-contact ends on people with no relationship
#                     met that day, and all new-contact ends
#   act.z.rel, act.zz.rel  people with no relationship contact on a day (mean of
#                     yesterday and today) and on both days, among people living
#                     on both: the symmetric rejoin is 1 - sum(zz) / sum(z)
#   act.z.all, act.zz.all  the same counting new contacts too
#   act.n2            people living on both days; with act.z.rel and act.zz.rel
#                     it gives b = P(contact both days | either) =
#                     (n2 - 2 z + zz) / (n2 - zz)
#   act.iso.rel.ageN  daily share with no relationship contact by age band
#   act.n, act.sF, act.sI, act.sFI, act.sF2, act.sI2  sums of the day's frequent
#                     and infrequent contacts per person, for their correlation
#   act.both.F.num/den, act.both.I.num/den  relationships met on both days and
#                     on either day, counted from each end with relationship
#                     contact on both days (the diary's conditioning)
#   act.union2, act.daily2  distinct contacts over yesterday and today, and the
#                     sum of the two days' contacts: the two-day ratio is
#                     2 * sum(union2) / sum(daily2)
activation_diagnostics <- function(dat, at, day, uid, out, live, degF, degI, degN, band) {
  n_live <- sum(live)
  zero_rel <- degF + degI == 0
  zero_all <- zero_rel & degN == 0
  # the daily zero share by band, and the sums behind the person-day correlation
  # of frequent and infrequent contacts
  for (g in seq_along(get_param(dat, "age.grps"))) {
    dat <- set_epi(dat, paste0("act.iso.rel.age", g), mean(zero_rel[live & band == g]), at = at)
  }
  x <- degF[live]; y <- degI[live]
  dat <- set_epi(dat, "act.n", n_live, at = at)
  dat <- set_epi(dat, "act.sF", sum(x), at = at)
  dat <- set_epi(dat, "act.sI", sum(y), at = at)
  dat <- set_epi(dat, "act.sFI", sum(x * y), at = at)
  dat <- set_epi(dat, "act.sF2", sum(x^2), at = at)
  dat <- set_epi(dat, "act.sI2", sum(y^2), at = at)
  dat <- set_epi(dat, "act.out", sum(out) / n_live, at = at)
  dat <- set_epi(dat, "act.deg.freq", 2 * nrow(day$freq) / n_live, at = at)
  dat <- set_epi(dat, "act.deg.infreq", 2 * nrow(day$infreq) / n_live, at = at)
  dat <- set_epi(dat, "act.deg.new", 2 * nrow(day$new) / n_live, at = at)
  dat <- set_epi(dat, "act.iso.rel", sum(zero_rel & live) / n_live, at = at)
  dat <- set_epi(dat, "act.iso.all", sum(zero_all & live) / n_live, at = at)
  # New-contact ends on people with no relationship met today (the diary's 51%
  # rural and 48% urban)
  dat <- set_epi(dat, "act.new.norel", sum(zero_rel[day$new]), at = at)
  dat <- set_epi(dat, "act.new.ends", length(day$new), at = at)

  prev <- dat$run$act_prev
  now <- list(uid = uid, out = out, live = live, zero_rel = zero_rel, zero_all = zero_all,
              F = matrix(uid[day$freq], ncol = 2), I = matrix(uid[day$infreq], ncol = 2),
              n_new = nrow(day$new))
  if (!is.null(prev)) {
    j <- match(uid, prev$uid)
    both <- live & !is.na(j)
    both[both] <- prev$live[j[both]]
    zr_prev <- prev$zero_rel[j]; za_prev <- prev$zero_all[j]
    dat <- set_epi(dat, "act.z.rel", (sum(zr_prev[both]) + sum(zero_rel[both])) / 2, at = at)
    dat <- set_epi(dat, "act.zz.rel", sum(zr_prev[both] & zero_rel[both]), at = at)
    dat <- set_epi(dat, "act.z.all", (sum(za_prev[both]) + sum(zero_all[both])) / 2, at = at)
    dat <- set_epi(dat, "act.zz.all", sum(za_prev[both] & zero_all[both]), at = at)
    dat <- set_epi(dat, "act.n2", sum(both), at = at)
    # The diary's conditioning: a participant with relationship contact on both
    # days. A tie met on either day is counted from each such end; an end that
    # has died since yesterday is not counted.
    out_both <- rep(FALSE, length(uid))
    out_both[both] <- !zero_rel[both] & !prev$zero_rel[j[both]]
    both_days <- function(cur, old) {
      kc <- pair_key(cur[, 1], cur[, 2]); ko <- pair_key(old[, 1], old[, 2])
      old_only <- !ko %in% kc
      seen <- rbind(cur, old[old_only, , drop = FALSE])
      met2 <- c(kc %in% ko, rep(FALSE, sum(old_only)))
      ca <- out_both[match(seen[, 1], uid)] %in% TRUE
      cb <- out_both[match(seen[, 2], uid)] %in% TRUE
      c(num = sum(met2 & ca) + sum(met2 & cb), den = sum(ca) + sum(cb), union = nrow(seen))
    }
    bF <- both_days(now$F, prev$F)
    bI <- both_days(now$I, prev$I)
    dat <- set_epi(dat, "act.both.F.num", bF[["num"]], at = at)
    dat <- set_epi(dat, "act.both.F.den", bF[["den"]], at = at)
    dat <- set_epi(dat, "act.both.I.num", bI[["num"]], at = at)
    dat <- set_epi(dat, "act.both.I.den", bI[["den"]], at = at)
    dat <- set_epi(dat, "act.union2", bF[["union"]] + bI[["union"]] + now$n_new + prev$n_new, at = at)
    dat <- set_epi(dat, "act.daily2", nrow(now$F) + nrow(prev$F) + nrow(now$I) + nrow(prev$I) +
                     now$n_new + prev$n_new, at = at)
  }
  dat$run$act_prev <- now
  dat
}

#' Relationships of newborns and of people who crossed a band top
#'
#' The relationship networks are dyad independent, so a node's ties at
#' equilibrium are independent draws with probability
#' expit(theta_gh + log a_i + log a_j), theta_gh the age pair's cross-sectional
#' coefficient. The initial population starts at that equilibrium; newborns
#' would otherwise start with none and reach it only over a relationship
#' duration (years). Drawing them at arrival makes newborns and the initial
#' infants exchangeable (GlobalMix India #65). Off when `act.seed.newborns` is
#' FALSE.
#'
#' A person who crossed a band top on this step (`rel.cross`, set by
#' [cross_band_tops()]) loses their ties and draws new ones the same way, from
#' the band they entered, so the band's relationships stay at its
#' cross-section (GlobalMix India #64).
#'
#' @inheritParams aging
#' @return The updated `dat`.
#' @keywords internal
seed_relationship_ties <- function(dat, at) {
  active <- get_attr(dat, "active")
  newborn <- if (isFALSE(as.logical(param_or(dat, "act.seed.newborns", TRUE)))) integer(0) else
    which(active == 1 & get_attr(dat, "entrTime") == at)
  rc <- get_attr(dat, "rel.cross", override.null.error = TRUE)
  crossed <- if (is.null(rc)) integer(0) else which(active == 1 & !is.na(rc) & rc == at)
  seeded <- unique(c(newborn, crossed))
  if (length(seeded) == 0) return(dat)
  never <- get_attr(dat, "never")
  band <- band_index(dat)
  elig <- which(active == 1 & never == 0)
  for (layer in rel_layers) {
    k <- match(layer, names(dat$run$el))
    cell <- get_nwparam(dat, network = k)$cell_theta
    if (is.null(cell)) {
      stop("the `", layer, "` network has no cell_theta; refit it with the relationship pipeline", call. = FALSE)
    }
    la <- get_attr(dat, if (layer == rel_layers[["freq"]]) "log_aF" else "log_aI")
    add <- vector("list", length(seeded))
    for (b in seq_along(seeded)) {
      i <- seeded[b]
      if (never[i] == 1) next
      # a pair of nodes seeded on the same day is drawn once, from the earlier one
      j <- elig[elig != i & !elig %in% seeded[seq_len(b)]]
      p <- stats::plogis(cell[band[i], band[j]] + la[i] + la[j])
      hit <- j[stats::runif(length(j)) < p]
      if (length(hit)) add[[b]] <- cbind(pmin(i, hit), pmax(i, hit))
    }
    add <- do.call(rbind, add)
    el <- dat$run$el[[k]]
    gone <- length(crossed) > 0 && any(el[, 1] %in% crossed | el[, 2] %in% crossed)
    if (is.null(add) && !gone) next
    keep_attr <- attributes(el)[setdiff(names(attributes(el)), c("dim", "dimnames"))]
    old <- matrix(as.integer(el), ncol = 2)
    if (length(crossed)) old <- old[!(old[, 1] %in% crossed | old[, 2] %in% crossed), , drop = FALSE]
    new <- unique(rbind(old, matrix(as.integer(add), ncol = 2)))
    new <- new[order(new[, 1], new[, 2]), , drop = FALSE]
    for (a in names(keep_attr)) attr(new, a) <- keep_attr[[a]]
    dat$run$el[[k]] <- new
  }
  dat
}
