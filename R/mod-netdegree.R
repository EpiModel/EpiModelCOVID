#' @rdname moduleset-gmc19
#' @export
netdegree <- function(dat, at) {
  if (isFALSE(dat$param$compute.degree)) return(dat)

  # Degree on every network, by name (degree_school, degree_work and
  # degree_nonhome on the four-layer model; degree_freq and degree_infreq on the
  # relationship model, GlobalMix India #63), and on the household layer.
  el <- dat$run$el
  hh <- hh_layer(dat)
  nets <- setdiff(seq_along(el), hh)
  n_nodes <- length(get_attr(dat, "active"))
  deg_net <- lapply(nets, function(i) degree_of(el[[i]], n_nodes))
  names(deg_net) <- names(el)[nets]
  degree_hh <- if (is.na(hh)) numeric(n_nodes) else degree_of(el[[hh]], n_nodes)

  # cross-layer (non-household) degree: prioritization metrics that exclude the
  # dense household clique, so "degree" is not dominated by household size.
  degree_xlayer <- Reduce(`+`, deg_net, numeric(n_nodes))
  degree_total <- degree_xlayer + degree_hh

  # cross-layer bridging
  n_layers_active <- Reduce(`+`, lapply(deg_net, function(d) d > 0), numeric(n_nodes)) +
    (degree_hh > 0)
  is_bridge <- n_layers_active >= 2 # bridge indicator

  # mean degree by age group. Every network tie is a contact every day on the
  # four-layer model; on the relationship model the day's contacts are the
  # relationships met today, the new contacts and the household (GlobalMix India
  # #63), so the age profile of degree is the profile of daily contacts, as in
  # the diary.
  age.grp <- get_attr(dat, item = "age.grp")
  daily_degree <- degree_total
  if (activation_on(dat) && identical(dat$run$contacts_at, at)) {
    daily_degree <- Reduce(`+`, lapply(dat$run$contacts, degree_of, n = n_nodes))
  }
  mean_degree_age <- tapply(daily_degree, age.grp, mean)

  # proportion of nodes that are cross-layer bridges
  ## overall
  prop_bridge <- mean(is_bridge)
  ## by age group
  prop_bridge_age <- tapply(is_bridge, age.grp, mean)

  # degree variance by age group
  var_degree_age <- tapply(daily_degree, age.grp, var)

  for (nm in names(deg_net)) {
    dat <- set_attr(dat, paste0("degree_", nm), deg_net[[nm]])
  }
  dat <- set_attr(dat, "degree_hh", degree_hh)
  dat <- set_attr(dat, "degree_total", degree_total)
  dat <- set_attr(dat, "degree_xlayer", degree_xlayer)
  dat <- set_attr(dat, "is_bridge", is_bridge)
  dat <- set_attr(dat, "n_layers_active", n_layers_active)

  for (nm in names(deg_net)) {
    dat <- set_epi(dat, paste0("mean_degree_", nm), at, mean(deg_net[[nm]]))
  }
  dat <- set_epi(dat, "mean_degree_hh", at, mean(degree_hh))
  dat <- set_epi(dat, "mean_degree_total", at, mean(degree_total))

  dat <- set_epi(dat, "prop_bridge", at, prop_bridge)
  for (g in seq_along(prop_bridge_age)) {
    dat <- set_epi(dat, paste0("prop_bridge_age", g), at, prop_bridge_age[[g]])
  }
  for (g in seq_along(mean_degree_age)) {
    dat <- set_epi(dat, paste0("mean_degree_age", g), at, mean_degree_age[[g]])
  }
  for (g in seq_along(var_degree_age)) {
    dat <- set_epi(dat, paste0("var_degree_age", g), at, var_degree_age[[g]])
  }

  if (activation_on(dat)) dat <- relationship_scores(dat, at)
  score <- degree_score(dat)
  dat <- set_attr(dat, "degree_score", score)

  # degree quartile (rank-based, 4 equal-size groups) for outcome stratification,
  # on the score the degree strategies rank on
  degree_quartile <- pmin(pmax(
    ceiling(rank(score, ties.method = "first") / n_nodes * 4), 1L), 4L)
  dat <- set_attr(dat, "degree_quartile", degree_quartile)
  dat <- set_epi(dat, "mean_degree_daily", at, mean(daily_degree))

  dat
}

degree_of <- function(el, n) {
  if (is.null(el) || length(el) == 0) return(numeric(n))
  tabulate(el, nbins = n)
}

# Sum of a per-edge value over each node's edges.
edge_sum <- function(el, x, n) {
  out <- numeric(n)
  if (length(x) == 0) return(out)
  s <- rowsum(c(x, x), c(el[, 1], el[, 2]))
  out[as.integer(rownames(s))] <- s[, 1]
  out
}

# Targeting scores on the relationship model (GlobalMix India #63, #75):
#   degree_active  expected relationships met per day: o_i times the sum over
#                  frequent relationships of o_j q and over infrequent ones of
#                  o_j s (a contact survey of usual contacts plus how often each
#                  is seen)
#   degree_hazard  expected daily exposure: each relationship's probability of
#                  being met times its age pair's relative hazard and the tie
#                  type's exposure.rate, plus the household ties, met every day
# New contacts cannot be anticipated and are in neither.
relationship_scores <- function(dat, at) {
  pars <- activation_params(dat)
  el <- dat$run$el
  n <- length(get_attr(dat, "active"))
  band <- band_index(dat)
  o_eff <- get_attr(dat, "o_i") * (1 - get_attr(dat, "never"))
  eF <- el[[rel_layers[["freq"]]]]; eI <- el[[rel_layers[["infreq"]]]]
  pF <- o_eff[eF[, 1]] * o_eff[eF[, 2]] * (pars$q[band[eF[, 1]]] + pars$q[band[eF[, 2]]]) / 2
  pI <- o_eff[eI[, 1]] * o_eff[eI[, 2]] * pars$s
  degree_active <- edge_sum(eF, pF, n) + edge_sum(eI, pI, n)

  weights <- contact_weights(dat)
  exposure.rate <- layer_param(dat, "exposure.rate")
  wt <- function(type, e) {
    w <- if (is.null(weights)) rep(1, nrow(e)) else weights[[type]][cbind(band[e[, 1]], band[e[, 2]])]
    w * exposure.rate[[type]]
  }
  hh <- hh_layer(dat)
  degree_hazard <- edge_sum(eF, pF * wt(rel_layers[["freq"]], eF), n) +
    edge_sum(eI, pI * wt(rel_layers[["infreq"]], eI), n)
  if (!is.na(hh)) {
    eH <- el[[hh]]
    degree_hazard <- degree_hazard + edge_sum(eH, wt(names(el)[hh], eH), n)
  }
  dat <- set_attr(dat, "degree_active", degree_active)
  dat <- set_attr(dat, "degree_hazard", degree_hazard)
  dat <- set_epi(dat, "mean_degree_active", at, mean(degree_active))
  dat
}

#' The score the degree strategies rank on
#'
#' `vax.degree.score` names it: "total" (every network tie, the household
#' included; the default), "xlayer" (non-household ties), and on the
#' relationship model "relationship" (frequent relationships: the people seen
#' at least weekly), "all" (frequent and infrequent relationships), "active"
#' (`degree_active`), "hazard" (`degree_hazard`), "recent" (out-of-home
#' contacts met over the last `act.window` days), "new" (new contacts over the
#' last `act.window` days) and "out" (o_i, the propensity to go out).
#'
#' @param dat Main `netsim_dat` object.
#' @return Numeric score per node.
#' @keywords internal
degree_score <- function(dat) {
  score <- as.character(param_or(dat, "vax.degree.score", "total"))
  rel_only <- c("relationship", "all", "active", "hazard", "recent", "new", "out")
  if (score %in% rel_only && !activation_on(dat)) {
    stop("vax.degree.score \"", score, "\" needs the relationship model", call. = FALSE)
  }
  switch(score,
    total = get_attr(dat, "degree_total"),
    xlayer = get_attr(dat, "degree_xlayer"),
    relationship = get_attr(dat, paste0("degree_", rel_layers[["freq"]])),
    all = get_attr(dat, paste0("degree_", rel_layers[["freq"]])) +
      get_attr(dat, paste0("degree_", rel_layers[["infreq"]])),
    active = get_attr(dat, "degree_active"),
    hazard = get_attr(dat, "degree_hazard"),
    recent = get_attr(dat, "act.recent"),
    new = get_attr(dat, "act.recent.new"),
    out = get_attr(dat, "o_i") * (1 - get_attr(dat, "never")),
    stop("unknown vax.degree.score \"", score, "\"", call. = FALSE))
}
