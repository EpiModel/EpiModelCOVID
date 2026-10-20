#' @rdname moduleset-gmc19
#' @export
infect_general <- function(dat, at) {

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")

  vax <- get_attr(dat, "vax")

  ## Find infected nodes ##
  idsInf <- which(active == 1 & status %in% c("a", "ic", "ip"))

  ## Common Parameters ##
  inf.prob.a.rr <- get_param(dat, "inf.prob.a.rr")
  vax.schedule <- build_vax_schedule(dat, get_pathogen(dat))

  # The vaccine's protection against infectiousness (GlobalMix India #73): a
  # vaccinated infector transmits at 1 - VE_infectious of the per-act
  # probability, on contacts and in the background hazard alike. Before this
  # the only effect of vaccination on onward transmission was the smaller
  # symptomatic share of breakthrough infections.
  inf_rr <- rep(1, length(active))
  if (length(idsInf) > 0) {
    inf_rr[idsInf] <- compute_ve(at = at, ids = idsInf, vax = vax,
                                 last.dose.time = get_attr(dat, "last.dose.time"),
                                 vax.schedule = vax.schedule, outcome = "infectious",
                                 vax.age.group = vax_age_group_for(dat))$rr
  }

  # Per-layer parameters, one value per layer and named by layer. layer_param()
  # stops if a vector no longer lines up with the layers, e.g. after an updater.
  inf.prob <- layer_param(dat, "inf.prob")
  exposure.rate <- layer_param(dat, "exposure.rate")
  # Symptomatic contact reduction (#67, #68): a per-contact hazard multiplier,
  # by tie type, while the infector is symptomatic ("ic"). Absent, no reduction.
  sympt.rr <- if (!is.null(get_param(dat, "sympt.contact.rr", override.null.error = TRUE))) {
    layer_param(dat, "sympt.contact.rr")
  }

  # Seasonal forcing of transmission (issue #26): a sinusoidal annual multiplier on
  # inf.prob. Defaults to no forcing (amplitude 0) when the params are unset.
  seasonal.amp <- get_param(dat, "seasonal.amp", override.null.error = TRUE)
  seasonal.phase <- get_param(dat, "seasonal.phase", override.null.error = TRUE)
  if (is.null(seasonal.amp)) seasonal.amp <- 0
  if (is.null(seasonal.phase)) seasonal.phase <- 0
  season_mult <- 1 + seasonal.amp * cos(2 * pi * (at - seasonal.phase) / days_per_year)

  # Tie types (GlobalMix India #63). Without the activation step these are the
  # networks, every tie of which is a contact every day. With it they are the
  # day's contacts by tie type from activate_contacts(): the relationships met
  # today, the new contacts, and the household layer.
  types <- tie_types(dat)
  contacts <- if (activation_on(dat)) dat$run$contacts else NULL
  if (activation_on(dat) && (is.null(contacts) || !identical(dat$run$contacts_at, at) ||
                             !identical(names(contacts), types))) {
    stop("no contacts for step ", at, " by tie type (", toString(types), "); is activation.FUN ",
         "in module.order, after the network modules and before infection.FUN?", call. = FALSE)
  }
  # Relative per-contact daily hazard by age pair and tie type (#68), or NULL
  weights <- contact_weights(dat)
  band <- if (!is.null(weights)) band_index(dat)

  nLayers <- length(types)
  nInf <- rep(0, nLayers)
  allNewInf <- integer(0)  # newly infected across all layers + imports, for stratified incidence

  # Daily transmission probability on each discordant pair (sus, inf) of one tie
  # type, with every per-contact multiplier. The transmission draws and the
  # discordant-contact counters below both use it, so they cannot drift apart.
  pair_prob <- function(del, layer) {
    # Per-act transmission probability on this layer (seasonally forced).
    # inf.prob and exposure.rate are the only terms that depend on the layer,
    # so when both are equal across layers every tie carries the same daily
    # hazard. The CorporateMix modifiers that used to sit in this loop are
    # removed (GlobalMix India #66): inf.add/inf.sub windows on days 40-80
    # and 120-160, a 40% cut from day 15 on every layer but the first, and
    # isolation and masking of cases found by the testing module.
    del$transProb <- inf.prob[[layer]] * season_mult

    # The susceptible's relative acquisition probability: the vaccine (dose
    # allocation is in mod-vax.R; the direct infant product carries no infection
    # VE, #41), the age band (susc.age.rr, the #45 knob and the #70 COVID-19
    # gradient), young infants (susc.infant.rr, #65) and immunity from an
    # earlier infection in the run (#69). susceptibility_rr() is in
    # mod-immunity.R; the background hazard and the imports take the same one.
    del$vaxSus <- vax[del$sus]
    del$transProb <- del$transProb * susceptibility_rr(dat, del$sus, at, vax.schedule) * inf_rr[del$inf]

    # Asymptomatic infection
    del$stat <- status[del$inf]
    del$transProb[del$stat == "a"] <- del$transProb[del$stat == "a"] *
      inf.prob.a.rr

    # Daily probability over exposure.rate acts on the tie, each act scaled
    # by the age pair's relative intensity on this tie type and, while the
    # infector is symptomatic, by the symptomatic contact multiplier
    w <- if (is.null(weights)) 1 else weights[[layer]][cbind(band[del$sus], band[del$inf])]
    rr_ic <- if (is.null(sympt.rr)) 1 else ifelse(del$stat == "ic", sympt.rr[[layer]], 1)
    del$finalProb <- 1 - (1 - del$transProb)^(exposure.rate[[layer]] * w * rr_ic)
    del
  }

  # Discordant contacts by tie type at the start of the step (GlobalMix India
  # #77, for the #68 check of infections by tie type against the weights). Per
  # tie type: the contacts between an infectious and a susceptible person
  # (disc.<type>), their summed daily hazard, -log(1 - p) (disc.haz.<type>), and
  # the infections the tie type carries if all of a person's exposures compete
  # at once (disc.exp.<type>: each susceptible is infected with probability
  # 1 - exp(-H), H their total hazard, attributed in proportion to hazard).
  # Taken before any transmission draw and with no random numbers, so the run
  # is unchanged.
  disc <- disc_counts(if (length(idsInf) > 0) lapply(seq_len(nLayers), function(layer) {
    el <- if (is.null(contacts)) get_edgelist(dat, layer) else contacts[[layer]]
    del <- discord_pairs(el, status, active, c("a", "ic", "ip"), at, shuffle = FALSE)
    if (!is.null(del)) pair_prob(del, layer)
  }), nLayers, length(active))

  if (length(idsInf) > 0) {
    for (layer in seq_len(nLayers)) {
      ## Look up discordant edgelist ##
      del <- if (is.null(contacts)) {
        discord_edgelist(dat, at, network = layer, infstat = c("a", "ic", "ip"))
      } else {
        discord_pairs(contacts[[layer]], get_attr(dat, "status"), active, c("a", "ic", "ip"), at)
      }

      ## If any discordant pairs, proceed ##
      if (!(is.null(del))
      ) {

        del <- pair_prob(del, layer)

        # rbinom() would turn these into NAs with only a warning, and the NA
        # draws would then drop out of the which() below as non-transmissions
        bad <- is.na(del$finalProb) | del$finalProb < 0 | del$finalProb > 1
        if (any(bad)) {
          stop("infect_general: ", sum(bad), " of ", nrow(del),
               " transmission probabilities on layer `", names(inf.prob)[layer],
               "` at step ", at, " are NA or outside [0, 1] (per-act range ",
               toString(signif(range(del$transProb), 3)), ", exposure.rate ",
               exposure.rate[[layer]], ")", call. = FALSE)
        }

        # Stochastic transmission process
        transmit <- rbinom(nrow(del), 1, del$finalProb)

        # Keep rows where transmission occurred
        del <- del[which(transmit == 1), , drop = FALSE]

        # Look up new ids if any transmissions occurred
        idsNewInf <- unique(del$sus)
        nInf[layer] <- length(idsNewInf)

        # Each infection starts an episode with its own severity (#69)
        if (nInf[layer] > 0) {
          dat <- start_episode(dat, idsNewInf, at)
          allNewInf <- c(allNewInf, idsNewInf)
        }
      }
    }
  }

  ## Background mass-action hazard (#68; the #79 group-contact arm), off when
  ## bg.hazard is 0 or unset. A susceptible in age band b is infected with daily
  ## probability 1 - exp(-bg.hazard x season x sum_h W_bh prev_h x rr): W is the
  ## row-normalized out-of-home age mixing (bg.mix, or the new-contact mixing
  ## new.mix), prev_h the infectious share of band h with asymptomatic infectors
  ## weighted by inf.prob.a.rr, and rr the susceptible's vaccine and
  ## age-susceptibility multipliers, as on a network contact.
  bg.hazard <- param_or(dat, "bg.hazard", 0)
  nBg <- 0L
  if (bg.hazard > 0) {
    status_now <- get_attr(dat, "status")
    bnd <- match(get_attr(dat, "age.grp"), get_param(dat, "age.grps"))
    nb <- length(get_param(dat, "age.grps"))
    live <- active == 1
    M <- get_param(dat, "bg.mix", override.null.error = TRUE)
    if (is.null(M)) M <- get_param(dat, "new.mix", override.null.error = TRUE)
    if (is.null(M)) stop("bg.hazard needs an age mixing matrix (bg.mix)", call. = FALSE)
    W <- M + diag(diag(M))
    W <- W / pmax(rowSums(W), .Machine$double.xmin)
    infw <- ((status_now %in% c("ip", "ic")) + inf.prob.a.rr * (status_now == "a")) * inf_rr
    prev <- as.numeric(tapply(infw[live], factor(bnd[live], levels = seq_len(nb)), sum)) /
      pmax(tabulate(bnd[live], nb), 1)
    prev[is.na(prev)] <- 0
    lambda_b <- bg.hazard * season_mult * as.numeric(W %*% prev)
    idsSus <- which(live & status_now == "s")
    if (length(idsSus) > 0) {
      rr <- susceptibility_rr(dat, idsSus, at, vax.schedule)
      idsBg <- idsSus[stats::runif(length(idsSus)) < 1 - exp(-lambda_b[bnd[idsSus]] * rr)]
      if (length(idsBg) > 0) {
        dat <- start_episode(dat, idsBg, at)
        nBg <- length(idsBg)
        allNewInf <- c(allNewInf, idsBg)
      }
    }
  }

  ## Infections from outside the population (#70). import.rate is a daily
  ## probability per susceptible, seasonally forced, all year round, so a
  ## seasonal epidemic can ignite whenever transmission allows it and zero
  ## prevalence is not absorbing. The introduction pulse adds import.pulse
  ## expected introductions into a fully susceptible population over
  ## import.pulse.days days from import.pulse.start: the introduction of a novel
  ## pathogen, or the seeding of a season. Both scale with the susceptible's
  ## relative acquisition probability, so the vaccine, the age gradient and
  ## immunity protect against an imported infection as against a local one
  ## (they used to ignore all three).
  import.rate <- param_or(dat, "import.rate", 0)
  pulse <- import_pulse_prob(dat, at)
  nImport <- 0L
  if (import.rate > 0 || pulse > 0) {
    status_now <- get_attr(dat, "status")
    idsSus <- which(active == 1 & status_now == "s")
    if (length(idsSus) > 0) {
      rr <- susceptibility_rr(dat, idsSus, at, vax.schedule)
      import_prob <- pmin(1, (import.rate * season_mult + pulse) * rr)
      idsImp <- idsSus[runif(length(idsSus)) < import_prob]
      if (length(idsImp) > 0) {
        dat <- start_episode(dat, idsImp, at)
        nImport <- length(idsImp)
        allNewInf <- c(allNewInf, idsImp)
      }
    }
  }
  allNewInf <- unique(allNewInf)

  ## Summary statistics for incidence (network by layer + imports)
  dat <- set_epi(dat, "se.flow", sum(nInf) + nImport + nBg, at = at)
  dat <- set_epi(dat, "se.flow.l1", nInf[1], at = at)
  dat <- set_epi(dat, "se.flow.l2", nInf[2], at = at)
  dat <- set_epi(dat, "se.flow.l3", nInf[3], at = at)
  dat <- set_epi(dat, "se.flow.l4", if (length(nInf) >= 4) nInf[4] else 0, at = at)  # household layer
  # The same counts named by tie type (GlobalMix India #63), so a reader never has
  # to know which position a layer took
  for (k in seq_along(types)) dat <- set_epi(dat, paste0("se.flow.", types[k]), nInf[k], at = at)
  for (k in seq_along(types)) {
    dat <- set_epi(dat, paste0("disc.", types[k]), disc$n[k], at = at)
    dat <- set_epi(dat, paste0("disc.haz.", types[k]), disc$haz[k], at = at)
    dat <- set_epi(dat, paste0("disc.exp.", types[k]), disc$exp[k], at = at)
  }
  dat <- set_epi(dat, "se.import.flow", nImport, at = at)
  if (bg.hazard > 0) dat <- set_epi(dat, "se.flow.bg", nBg, at = at)

  ## First infections and reinfections (#69, the #77 se.flow.first): a first
  ## infection is a person's first episode in the run. Persons ever infected
  ## are the primary infection outcome (#58 item 5); se.flow counts episodes.
  first <- get_attr(dat, "inf.count")[allNewInf] == 1
  dat <- set_epi(dat, "se.flow.first", sum(first), at = at)
  dat <- set_epi(dat, "se.flow.reinf", sum(!first), at = at)

  ## Stratified incidence: by age group, vaccination status, and network position
  age <- get_attr(dat, "age")
  age.breaks <- get_param(dat, "age.breaks")
  ag <- cut(age[allNewInf], age.breaks, labels = FALSE, right = FALSE)
  for (g in seq_len(length(age.breaks) - 1)) {
    dat <- set_epi(dat, paste0("se.flow.age", g), sum(ag == g, na.rm = TRUE), at = at)
    dat <- set_epi(dat, paste0("se.flow.first.age", g), sum(ag == g & first, na.rm = TRUE), at = at)
  }
  if (sev_diagnostics(dat)) dat <- sev_flows(dat, at, "inf", allNewInf, first)
  if (vax_diagnostics(dat)) dat <- vd_flows(dat, at, "inf", allNewInf)
  vaxed_new <- vax[allNewInf] >= 1
  dat <- set_epi(dat, "se.flow.vax", sum(vaxed_new, na.rm = TRUE), at = at)
  dat <- set_epi(dat, "se.flow.unvax", sum(!vaxed_new, na.rm = TRUE), at = at)
  # Age x vaccination-status cross-tab. The paper's indirect-protection claim
  # (vaccinating network hubs cuts infection among the UNVACCINATED vulnerable)
  # cannot be read from the age and vax marginals alone, so emit incidence among
  # the unvaccinated (and vaccinated) within each age group.
  for (g in seq_len(length(age.breaks) - 1)) {
    dat <- set_epi(dat, paste0("se.flow.unvax.age", g), sum(ag == g & !vaxed_new, na.rm = TRUE), at = at)
    dat <- set_epi(dat, paste0("se.flow.vax.age", g),   sum(ag == g &  vaxed_new, na.rm = TRUE), at = at)
  }
  degree_quartile <- get_attr(dat, "degree_quartile")
  if (!is.null(degree_quartile)) {
    dq_new <- degree_quartile[allNewInf]
    for (q in 1:4) dat <- set_epi(dat, paste0("se.flow.degQ", q), sum(dq_new == q, na.rm = TRUE), at = at)
  }
  is_bridge <- get_attr(dat, "is_bridge")
  if (!is.null(is_bridge)) {
    br_new <- is_bridge[allNewInf]
    dat <- set_epi(dat, "se.flow.bridge", sum(br_new, na.rm = TRUE), at = at)
    dat <- set_epi(dat, "se.flow.nonbridge", sum(!br_new, na.rm = TRUE), at = at)
  }
  # Infant incidence (issue #23): reported separately (and by vax status) so the
  # infant burden and the indirect protection of infants by parent-targeting are
  # not diluted in the coarse youngest age band.
  is_infant <- get_attr(dat, "is_infant")
  if (!is.null(is_infant)) {
    inf_new <- is_infant[allNewInf]
    dat <- set_epi(dat, "se.flow.infant", sum(inf_new, na.rm = TRUE), at = at)
    dat <- set_epi(dat, "se.flow.infant.unvax", sum(inf_new & !vaxed_new, na.rm = TRUE), at = at)
    dat <- set_epi(dat, "se.flow.infant.vax", sum(inf_new & vaxed_new, na.rm = TRUE), at = at)
    # Infant severity is sharply front-loaded, so the RSV infant contrast has to
    # be readable at the 0-6 / 6-12 month split rather than pooled under 1y
    # (issue #46). Deaths and hospitalisations already split this way; incidence
    # did not, which left the infant comparison undecomposable into its
    # attack-rate and severity parts.
    a_new <- age[allNewInf]
    dat <- set_epi(dat, "se.flow.infant.young", sum(inf_new & a_new <  0.5, na.rm = TRUE), at = at)
    dat <- set_epi(dat, "se.flow.infant.old",   sum(inf_new & a_new >= 0.5, na.rm = TRUE), at = at)
  }

  ## Person-level cumulative infection (issue #46) ------------------------------
  ## `se.flow` counts S->E transitions, i.e. EPISODES. With reinfection on, a
  ## cumulative per-100 built from it can and does exceed 100, which is not an
  ## attack rate and cannot be set against any seroprevalence or cohort estimate.
  ## `ever.inf` makes the person the unit: set once, never cleared, so summing it
  ## over the living population is the number of people ever infected.
  ever.inf <- get_attr(dat, "ever.inf")
  if (!is.null(ever.inf) && length(allNewInf) > 0) {
    ever.inf[allNewInf] <- 1
    dat <- set_attr(dat, "ever.inf", ever.inf)
  }

  ## Household secondary attack rate (issue #46) --------------------------------
  ## The one external anchor that survived the validation review is the PHIRST
  ## cross-pathogen household ordering, which is a per-person risk among exposed
  ## household contacts. The household-layer incidence is the numerator already,
  ## but with no denominator it is a count, not a rate.
  ##
  ## Denominator: susceptible nodes sharing a household with at least one
  ## infectious node at this step. Summing numerator and denominator over the
  ## epidemic and dividing gives a household SAR directly comparable to a cohort
  ## study's. Emitted as two series rather than a ratio so the analysis can
  ## aggregate before dividing; a per-step ratio would be undefined whenever no
  ## household is exposed.
  hh.ids <- get_attr(dat, "hh.ids")
  if (!is.null(hh.ids)) {
    status_now <- get_attr(dat, "status")
    # NA is "in no household", not a household, so it never matches
    hh_infectious <- unique(hh.ids[active == 1 & status_now %in% c("a", "ip", "ic", "h") &
                                     !is.na(hh.ids)])
    exposed_sus <- active == 1 & status_now == "s" & hh.ids %in% hh_infectious
    dat <- set_epi(dat, "hh.sar.den", sum(exposed_sus, na.rm = TRUE), at = at)
    hh <- hh_layer(dat)
    hh_type <- if (is.na(hh)) NA_integer_ else match(names(dat$run$el)[hh], types)
    dat <- set_epi(dat, "hh.sar.num", if (is.na(hh_type)) 0 else nInf[hh_type], at = at)
    dat <- set_epi(dat, "hh.exposed.n", length(hh_infectious), at = at)
  }

  return(dat)
}


# Discordant pairs of an edgelist of the day's contacts, as EpiModel's
# discord_edgelist() builds them from a network: rows in random order, the
# susceptible end first, both ends active. shuffle = FALSE keeps the edgelist
# order and draws no random numbers.
discord_pairs <- function(el, status, active, infstat, at, shuffle = TRUE) {
  if (is.null(el) || nrow(el) == 0) return(NULL)
  if (shuffle) el <- el[sample(seq_len(nrow(el))), , drop = FALSE]
  stat <- matrix(status[el], ncol = 2)
  isInf <- matrix(stat %in% infstat, ncol = 2)
  isSus <- matrix(stat %in% "s", ncol = 2)
  SIpairs <- el[isSus[, 1] * isInf[, 2] == 1, , drop = FALSE]
  ISpairs <- el[isSus[, 2] * isInf[, 1] == 1, , drop = FALSE]
  pairs <- rbind(SIpairs, ISpairs[, 2:1])
  if (nrow(pairs) == 0) return(NULL)
  del <- data.frame(at = at, sus = pairs[, 1], inf = pairs[, 2])
  del <- del[active[del$sus] == 1 & active[del$inf] == 1, , drop = FALSE]
  if (nrow(del) == 0) NULL else del
}

# Discordant-contact counters by tie type from each tie type's discordant pairs
# with their daily transmission probability (NULL for none): the number of
# pairs, the summed hazard -log(1 - p), and the competing-risk expected
# infections, sum over susceptibles of (1 - exp(-H)) h / H. A probability of 1
# is held just below it, so its hazard stays finite.
disc_counts <- function(dels, nLayers, n) {
  out <- list(n = numeric(nLayers), haz = numeric(nLayers), exp = numeric(nLayers))
  if (is.null(dels)) return(out)
  for (k in seq_len(nLayers)) {
    if (!is.null(dels[[k]])) dels[[k]]$h <- -log1p(-pmin(dels[[k]]$finalProb, 1 - 1e-12))
  }
  sus <- unlist(lapply(dels, function(d) d$sus))
  if (length(sus) == 0) return(out)
  h <- unlist(lapply(dels, function(d) d$h))
  type <- rep(seq_len(nLayers), vapply(dels, NROW, integer(1)))
  H <- numeric(n)
  agg <- rowsum(h, sus)
  H[as.integer(rownames(agg))] <- agg[, 1]
  Hi <- H[sus]
  share <- ifelse(Hi > 0, -expm1(-Hi) * h / Hi, 0)
  out$n <- tabulate(type, nLayers)
  out$haz <- vapply(seq_len(nLayers), function(k) sum(h[type == k]), numeric(1))
  out$exp <- vapply(seq_len(nLayers), function(k) sum(share[type == k]), numeric(1))
  out
}

# Relative per-contact daily hazard by age pair (#68): `contact.weight` is a
# list with one age band by age band matrix per tie type. NULL when unset or
# when `contact.weight.equal` is TRUE (the equal-weights sensitivity).
contact_weights <- function(dat) {
  w <- get_param(dat, "contact.weight", override.null.error = TRUE)
  if (is.null(w) || isTRUE(as.logical(param_or(dat, "contact.weight.equal", FALSE)))) return(NULL)
  types <- tie_types(dat)
  nb <- length(get_param(dat, "age.grps"))
  if (!is.list(w) || !all(types %in% names(w))) {
    stop("`contact.weight` needs a matrix for each tie type (", toString(types), ")", call. = FALSE)
  }
  w <- w[types]
  ok <- vapply(w, function(m) is.matrix(m) && all(dim(m) == nb) && all(is.finite(m)) && all(m >= 0), logical(1))
  if (!all(ok)) {
    stop("`contact.weight` for ", toString(types[!ok]), " is not a non-negative ", nb, " x ", nb,
         " matrix", call. = FALSE)
  }
  w
}


infect_covid_corporate <- function(dat, at) { #maria's repo
  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  dxStatus <- get_attr(dat, "dxStatus")
  isolate <- get_attr(dat, "isolate")
  
  vax <- get_attr(dat, "vax")
  vax1Time <- get_attr(dat, "vax1Time")
  vax2Time <- get_attr(dat, "vax2Time")
  vax3Time <- get_attr(dat, "vax3Time")
  vax4Time <- get_attr(dat, "vax4Time")
  
  ## Find infected nodes ##
  idsInf <- which(active == 1 & status %in% c("a", "ic", "ip"))
  
  ## Common Parameters ##
  inf.prob.a.rr <- get_param(dat, "inf.prob.a.rr")
  vax1.rr.infect <- get_param(dat, "vax1.rr.infect")
  vax2.rr.infect <- get_param(dat, "vax2.rr.infect")
  vax3.rr.infect <- get_param(dat, "vax3.rr.infect")
  vax4.rr.infect <- get_param(dat, "vax4.rr.infect")
  half.life <- get_param(dat, "half.life")
  inf.add <- get_param(dat, "inf.add")
  inf.sub <- get_param(dat, "inf.sub")
  inf.boost.start <- get_param(dat, "inf.boost.start")
  inf.boost.stop <- get_param(dat, "inf.boost.stop")
  inf.supp.start <- get_param(dat, "inf.supp.start")
  inf.supp.stop <- get_param(dat, "inf.supp.stop")
  inf.prob.mask.rr <- get_param(dat, "inf.prob.mask.rr")
  act.rate.iso.inter.time <- get_param(dat, "act.rate.iso.inter.time")
  act.rate.iso.inter.rr <- get_param(dat, "act.rate.iso.inter.rr")
  
  nLayers <- dat$num.nw
  nInf <- rep(0, nLayers)
  
  if (length(idsInf) > 0) {
    for (layer in seq_len(nLayers)) {
      ## Look up discordant edgelist ##
      del <- discord_edgelist(dat, at, network = layer,
                              infstat = c("a", "ic", "ip"))
      
      ## If any discordant pairs, proceed ##
      if (!(is.null(del))
          ) {
        
        ## Parameters ##
        inf.prob <- get_param(dat, "inf.prob")[layer]
        act.rate <- get_param(dat, "act.rate")[layer]
        inf.prob.inter.rr <- get_param(dat, "inf.prob.inter.rr")[layer]
        inf.prob.inter.time <- get_param(dat, "inf.prob.inter.time")[layer]
        act.rate.inter.rr <- get_param(dat, "act.rate.inter.rr")[layer]
        act.rate.inter.time <- get_param(dat, "act.rate.inter.time")[layer]
        
        # Update inf.prob to account for seasonal variation
        for (i in seq_along(inf.boost.start)) {
          if (at >= inf.boost.start[i] & at <= inf.boost.stop[i]){
            inf.prob <- inf.prob + inf.add[i]
          }
        }
        
        for (i in seq_along(inf.supp.start)) {
          if (at >= inf.supp.start[i] & at <= inf.supp.stop[i])
            inf.prob <- inf.prob - inf.sub[i]
        }
        
        # Set parameters on discordant edgelist data frame
        del$transProb <- inf.prob
        
        # Vaccination
        del$vaxSus <- vax[del$sus]
        del$transProb[del$vaxSus == 1] <- del$transProb[del$vaxSus == 1] * vax1.rr.infect
        del$transProb[del$vaxSus == 2] <- del$transProb[del$vaxSus == 2] * vax2.rr.infect
        del$transProb[del$vaxSus == 3] <- del$transProb[del$vaxSus == 3] * vax3.rr.infect
        del$transProb[del$vaxSus == 4] <- del$transProb[del$vaxSus == 4] * vax4.rr.infect
        
        #Waning vaccine immunity
        sinceVax1 <- at - vax1Time[del$sus] # days since 1st vax
        sinceVax2 <- at - vax2Time[del$sus] # days since 2nd vax
        sinceVax3 <- at - vax3Time[del$sus] # days since 3rd vax
        sinceVax4 <- at - vax4Time[del$sus] # days since 4th vax
        
        latest.vax <- pmin(sinceVax1, sinceVax2, sinceVax3, sinceVax4, na.rm = TRUE) # identify the most recent vax by taking the minimum time since vaccination across doses.
        latest.vax[is.na(latest.vax)] <- 0 #  individuals with no vax are assigned 0 - no waning 
        
        del$latest.vax <- latest.vax
        del$transProb <- pmin(del$transProb * (2 ^ (del$latest.vax / half.life)), inf.prob) # Vaccine protection as exponential waning: the transmission probability increases by a factor of 2 every half.life time units since the most recent dose, capped at the baseline transmission probability (inf.prob)
        
        # Asymptomatic infection
        del$stat <- status[del$inf]
        del$transProb[del$stat == "a"] <- del$transProb[del$stat == "a"] *
          inf.prob.a.rr
        
        # Generic inf.prob and act.rate interventions
        if (at >= inf.prob.inter.time) {
          del$transProb <- del$transProb * inf.prob.inter.rr
        }
        del$actRate <- act.rate
        if (at >= act.rate.inter.time) {
          del$actRate <- del$actRate * act.rate.inter.rr
        }
        
        # Case isolation for those in isolation process
        del$iso <- ifelse(!is.na(isolate[del$inf]),isolate[del$inf],0)
        if (at >= act.rate.iso.inter.time) {
          del$actRate[del$iso %in% c(1,2)] <- del$actRate[del$iso %in% c(1,2)] *
            act.rate.iso.inter.rr
        }
        
        # Masking for those in isolation process
        if (at >= act.rate.iso.inter.time) {
          del$transProb[del$iso %in% c(1,2,3,4)] <- del$transProb[del$iso %in% c(1,2,3,4)] *
            inf.prob.mask.rr
        }
        
        del$finalProb <- 1 - (1 - del$transProb)^del$actRate
        
        # Stochastic transmission process
        transmit <- rbinom(nrow(del), 1, del$finalProb)
        
        # Keep rows where transmission occurred
        del <- del[which(transmit == 1), , drop = FALSE]
        
        # Look up new ids if any transmissions occurred
        idsNewInf <- unique(del$sus)
        nInf[layer] <- length(idsNewInf)
        
        # Set new attributes for those newly infected
        if (nInf[layer] > 0) {
          dat <- set_attr(dat, "status", "e", idsNewInf)
          dat <- set_attr(dat, "infTime", at, idsNewInf)
          dat <- set_attr(dat, "statusTime", at, idsNewInf)
        }
      }
    }
  }
  
  ## Summary statistics for incidence
  dat$epi$se.flow[at] <- sum(nInf)
  dat$epi$se.flow.l1[at] <- nInf[1]
  dat$epi$se.flow.l2[at] <- nInf[2]
  dat$epi$se.flow.l3[at] <- nInf[3]
  
  return(dat)
}

# Relative per-contact acquisition probability of an infant younger than
# susc.infant.age years: a place for protection by maternal antibody (GlobalMix
# India #65). Whether the infant severity inputs are already net of it, and the
# value if not, are #71's; at the default, 1, it has no effect and returns NULL.
infant_susc_rr <- function(dat, age) {
  rr <- get_param(dat, "susc.infant.rr", override.null.error = TRUE)
  if (is.null(rr) || isTRUE(all(rr == 1))) return(NULL)
  cut_age <- get_param(dat, "susc.infant.age", override.null.error = TRUE)
  if (!(is.numeric(rr) && length(rr) == 1 && is.finite(rr) && rr >= 0) ||
      !(is.numeric(cut_age) && length(cut_age) == 1 && is.finite(cut_age) && cut_age > 0)) {
    stop("susc.infant.rr must be one non-negative relative risk and susc.infant.age one ",
         "positive age in years (GlobalMix India #65)", call. = FALSE)
  }
  ifelse(age < cut_age, rr, 1)
}
