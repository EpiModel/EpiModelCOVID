# Exact daily step of the relationship networks (GlobalMix India #63)
#
# Both relationship networks are separable-temporal ERGMs whose formation and
# dissolution models are dyad independent: a tie between i and j forms on a day
# with probability expit(phi_gh + log a_i + log a_j), phi_gh the formation
# coefficient of the age pair, and an existing tie persists with probability
# expit(psi_gh). tergm advances such a network by MCMC, and its per-step chain
# (burnin at most 1e5 proposals) almost never proposes the pairs of very active
# people again once their tie has dissolved, so those ties under-form and the
# density of the infrequent network fell 4-6% over a season (design review,
# 10/12/26). Because the models are dyad independent, the day can instead be
# drawn exactly: each tie persists by its own Bernoulli draw, and new ties form
# among the dyads that were not ties, by Poisson thinning (below).

#' Network resimulation for the relationship model
#'
#' The `resim_nets.FUN` of the relationship model. Applies EpiModel's edges
#' correction for the change in population size, then advances each relationship
#' network exactly with [step_relationships()] and any other layer with
#' EpiModel's own simulation (the household clique layer has none).
#'
#' @inheritParams aging
#' @return The updated `dat`.
#' @rdname moduleset-gmc19
#' @export
resim_relationships <- function(dat, at) {
  if (!is.null(get_control(dat, "dat.updates", override.null.error = TRUE))) {
    stop("resim_relationships() does not run dat.updates; drop it from the control", call. = FALSE)
  }
  if (any(unlist(get_control(dat, "tergmLite.track.duration", override.null.error = TRUE)))) {
    stop("resim_relationships() does not track tie durations (tergmLite.track.duration)", call. = FALSE)
  }
  active <- get_attr(dat, "active")
  if (any(active == 1) && isTRUE(get_control(dat, "resimulate.network"))) {
    dat <- EpiModel::edges_correct(dat, at)
    dat <- update_network_counts(dat)
    for (network in seq_len(dat$num.nw)) {
      if (activation_on(dat) && names(dat$run$el)[network] %in% rel_layers) {
        dat <- step_relationships(dat, at, network)
      } else {
        dat <- EpiModel::simulate_dat(dat = dat, at = at, network = network)
      }
    }
  }
  truncate <- get_control(dat, "truncate.el.cuml", override.null.error = TRUE)
  for (network in seq_len(dat$num.nw)) {
    dat <- EpiModel::update_cumulative_edgelist(dat, network, if (is.null(truncate)) 0 else truncate)
  }
  dat
}

# The population counts EpiModel's edges correction compares from one step to
# the next (as its internal update_sim_num(), for a one-group model).
update_network_counts <- function(dat) {
  active <- get_attr(dat, "active")
  dat$run$num <- sum(active == 1)
  attr.name <- get_control(dat, "edges.correct.attr", override.null.error = TRUE)
  if (!is.null(attr.name)) {
    e <- get_attr(dat, attr.name)
    dat$run$num.elig <- sum(active == 1 & !is.na(e) & e == 1)
  }
  dat
}

# Age-pair matrix of a coefficient vector in ergm's nodemix order: the first
# value is the reference cell's (carried by edges), the rest are the other cells
# in `cell_names` order, as differences from it.
cell_matrix <- function(coef, cell_names, ref, bands) {
  M <- matrix(NA_real_, length(bands), length(bands), dimnames = list(bands, bands))
  vals <- numeric(length(cell_names))
  vals[ref] <- coef[1]
  vals[-ref] <- coef[1] + coef[seq_len(length(cell_names) - 1) + 1]
  pr <- strsplit(cell_names, ".", fixed = TRUE)
  for (k in seq_along(pr)) {
    M[pr[[k]][1], pr[[k]][2]] <- vals[k]
    M[pr[[k]][2], pr[[k]][1]] <- vals[k]
  }
  M
}

#' Advance one relationship network by one day, exactly
#'
#' Every tie of yesterday persists with probability expit(psi_gh), from the
#' dissolution model's Persist coefficients (coef.adj, with the run's deaths
#' netted out). Every pair that was not a tie forms one with probability
#' expit(phi_gh + log a_i + log a_j), from the formation coefficients (coef.form,
#' after EpiModel's edges correction); never-active people form none. The
#' formation draw uses Poisson thinning: for each age pair, a Poisson number of
#' candidate pairs with mean exp(phi_gh) S_g S_h (S the sum of activities), ends
#' drawn in proportion to activity, each kept with probability
#' log(1 + e^x) / e^x, x = phi_gh + log a_i + log a_j, which makes the chance that
#' a pair forms at least one tie exactly expit(x). Candidates that were ties
#' yesterday, pairs of a person with themself and repeats are dropped.
#'
#' @inheritParams aging
#' @param network Index of the relationship network.
#' @return The updated `dat`.
#' @keywords internal
step_relationships <- function(dat, at, network) {
  nwp <- get_nwparam(dat, network = network)
  tg <- nwp$targets
  if (is.null(tg$cell_names) || is.null(nwp$coef.diss$cell_duration)) {
    stop("the `", names(dat$run$el)[network], "` network lacks its age-pair coefficients; ",
         "refit it with the relationship pipeline", call. = FALSE)
  }
  bands <- get_param(dat, "age.grps")
  k <- length(tg$cell_names)
  phi <- cell_matrix(nwp$coef.form[seq_len(k)], tg$cell_names, tg$ref, bands)
  diss_cells <- names(nwp$coef.diss$cell_duration)
  # the dissolution coefficients list the reference cell first
  psi <- cell_matrix(nwp$coef.diss$coef.adj, diss_cells, 1L, bands)

  el <- dat$run$el[[network]]
  keep_attr <- attributes(el)[setdiff(names(attributes(el)), c("dim", "dimnames"))]
  n <- length(get_attr(dat, "active"))
  band <- band_index(dat)
  attr <- if (identical(names(dat$run$el)[network], rel_layers[["freq"]])) "log_aF" else "log_aI"
  la <- get_attr(dat, attr)
  elig <- get_attr(dat, "active") == 1 & get_attr(dat, "never") == 0

  e1 <- as.integer(el[, 1]); e2 <- as.integer(el[, 2])
  keep <- stats::runif(length(e1)) < stats::plogis(psi[cbind(band[e1], band[e2])])
  old_key <- pair_key(e1, e2)
  new <- form_dyadind(band, la, elig, phi, old_key)
  a <- c(e1[keep], new[, 1]); b <- c(e2[keep], new[, 2])
  o <- order(a, b)
  out <- cbind(a[o], b[o])
  storage.mode(out) <- "integer"
  for (nm in names(keep_attr)) attr(out, nm) <- keep_attr[[nm]]
  dat$run$el[[network]] <- out
  dat
}

# New ties among pairs of eligible nodes that were not ties (keys in old_key),
# by the Poisson thinning described in step_relationships(). Returns a
# two-column matrix with tail < head.
form_dyadind <- function(band, la, elig, phi, old_key) {
  nb <- nrow(phi)
  ids <- lapply(seq_len(nb), function(g) which(elig & band == g))
  a <- exp(la)
  S <- vapply(ids, function(i) sum(a[i]), numeric(1))
  ta <- list(); tb <- list(); m <- 0
  for (g in seq_len(nb)) for (h in g:nb) {
    if (!length(ids[[g]]) || !length(ids[[h]])) next
    lam <- exp(phi[g, h]) * S[g] * S[h] * (if (g == h) 0.5 else 1)
    nc <- stats::rpois(1, lam)
    if (nc == 0) next
    i <- ids[[g]][sample.int(length(ids[[g]]), nc, replace = TRUE, prob = a[ids[[g]]])]
    j <- ids[[h]][sample.int(length(ids[[h]]), nc, replace = TRUE, prob = a[ids[[h]]])]
    x <- phi[g, h] + la[i] + la[j]
    # log(1 + e^x) / e^x, written so it neither overflows nor loses precision
    acc <- stats::runif(nc) < ifelse(x > 30, x * exp(-x), ifelse(x < -30, 1, log1p(exp(x)) / exp(x)))
    m <- m + 1
    ta[[m]] <- i[acc]; tb[[m]] <- j[acc]
  }
  i <- unlist(ta); j <- unlist(tb)
  if (!length(i)) return(matrix(integer(0), ncol = 2))
  p <- cbind(pmin(i, j), pmax(i, j))
  p <- p[p[, 1] != p[, 2], , drop = FALSE]
  key <- pair_key(p[, 1], p[, 2])
  p[!duplicated(key) & !key %in% old_key, , drop = FALSE]
}

# The relationship networks' starting cross-section, drawn exactly: every pair of
# eligible people is a tie with probability expit(theta_gh + log a_i + log a_j),
# theta_gh the age pair's cross-sectional coefficient stored with the fit.
start_relationships <- function(dat) {
  band <- band_index(dat)
  elig <- get_attr(dat, "active") == 1 & get_attr(dat, "never") == 0
  # the eligible population the cross-section is drawn for (cell_theta_shift())
  dat$run$rel.elig.start <- sum(elig)
  for (layer in rel_layers) {
    k <- match(layer, names(dat$run$el))
    cell <- get_nwparam(dat, network = k)$cell_theta
    if (is.null(cell)) {
      stop("the `", layer, "` network has no cell_theta; refit it with the relationship pipeline", call. = FALSE)
    }
    la <- get_attr(dat, if (layer == rel_layers[["freq"]]) "log_aF" else "log_aI")
    el <- dat$run$el[[k]]
    keep_attr <- attributes(el)[setdiff(names(attributes(el)), c("dim", "dimnames"))]
    new <- form_dyadind(band, la, elig, unname(cell), numeric(0))
    new <- new[order(new[, 1], new[, 2]), , drop = FALSE]
    storage.mode(new) <- "integer"
    for (a in names(keep_attr)) attr(new, a) <- keep_attr[[a]]
    dat$run$el[[k]] <- new
  }
  dat
}
