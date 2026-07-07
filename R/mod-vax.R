#' @rdname moduleset-gmc19
#' @export
vax_general <- function(dat, at) {

  ######## extract attribute ########
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  vax <- get_attr(dat, "vax")
  last.dose.time <- get_attr(dat, "last.dose.time") # timestep of each node's most recent dose

  vax.age.group <- vax_age_group_for(dat)

  dxStatus <- get_attr(dat, "dxStatus")
  dxTime <- get_attr(dat, "dxTime")
  
  degree_total <- get_attr(dat, "degree_total") # total nodal degree
  is_bridge <- get_attr(dat, "is_bridge")
  n_layers_active <- get_attr(dat, "n_layers_active")
  hh.ids <- get_attr(dat, "hh.ids")             # household id (parent-targeting)
  is_infant <- get_attr(dat, "is_infant")       # infant flag (parent-targeting)

  ######## extract parameters ########
  vax.strategy <- get_param(dat, "vax.strategy") # determine which stategy to use
  vax.supply.rate <- get_param(dat, "vax.supply.rate")
  vax.supply.total <- get_param(dat, "vax.supply.total")
  # Infant eligibility for the (adult, actively-administered) vaccine. For RSV the
  # adult vaccine is never given to infants (they are protected only by passive
  # maternal/monoclonal products and indirectly via household contacts), so
  # vax.infant.eligible = FALSE removes is_infant nodes from the eligible pool.
  # Default TRUE preserves covid/flu behaviour. infant_hh still targets the adult
  # co-residents of infants, so gating infants here does not affect it.
  vax.infant.eligible <- get_param(dat, "vax.infant.eligible", override.null.error = TRUE)
  if (is.null(vax.infant.eligible)) vax.infant.eligible <- TRUE
  season.length <- get_param(dat, "vax.season.length", override.null.error = TRUE)
  if (is.null(season.length)) season.length <- 364  # days per annual vaccination season
  # Imperfect hub-targeting (issue #40 / degraded-uptake arm): fraction of doses
  # that reach the intended top-ranked targets under a ranked network strategy;
  # the remainder are mis-allocated to random eligibles. Default 1 = perfect
  # targeting (unchanged behaviour).
  vax.target.uptake <- get_param(dat, "vax.target.uptake", override.null.error = TRUE)
  if (is.null(vax.target.uptake)) vax.target.uptake <- 1

  # Parallel multi-criteria hybrid (issue #41 / age_plus_* arms): fraction `f` of
  # each timestep's doses reserved for the 65+ elderly, with the remaining 1-f
  # allocated by the transmission criterion (the field proxy for age_plus_proxy,
  # the idealized true degree for the age_plus_degree ceiling). Default 0 = no
  # reserve, which leaves every non-hybrid strategy unchanged.
  vax.hybrid.frac <- get_param(dat, "vax.hybrid.frac", override.null.error = TRUE)
  if (is.null(vax.hybrid.frac)) vax.hybrid.frac <- 0

  # Field-deployable proxy score (issue #40), computed once per call from cheaply
  # observable covariates only (NOT the true contact degree). NULL for every
  # non-proxy strategy. proxy_hh ranks by active household size; proxy_fit ranks
  # by the predicted degree from a regression of true degree on household size,
  # school enrollment and employment (the best linear proxy a planner could fit
  # from survey data, then apply using only the observables).
  proxy_score <- NULL
  if (vax.strategy %in% c("proxy_hh", "proxy_fit", "age_plus_proxy")) {
    active_hh <- hh.ids[active == 1]
    hh_counts <- table(active_hh)
    hh_size <- as.numeric(hh_counts[match(hh.ids, names(hh_counts))]) # active co-residents
    if (vax.strategy == "proxy_hh") {
      proxy_score <- hh_size
    } else {
      degree_school <- get_attr(dat, "degree_school")
      degree_work   <- get_attr(dat, "degree_work")
      student  <- as.integer(degree_school > 0)  # enrolled (observable)
      employed <- as.integer(degree_work > 0)    # employed (observable)
      fit_df <- data.frame(deg = degree_total, hh_size = hh_size,
                           student = student, employed = employed)
      ok <- active == 1 & is.finite(degree_total) & is.finite(hh_size)
      proxy_lm <- stats::lm(deg ~ hh_size + student + employed,
                            data = fit_df[ok, , drop = FALSE])
      proxy_score <- as.numeric(stats::predict(proxy_lm, newdata = fit_df))
    }
  }

  vax.schedule <- build_vax_schedule(dat, get_pathogen(dat)) # data.frame containing vax-related details
  # vax.schedule stores both dose-administration parameters and per-dose RR values.
  # mod-vax.R only uses the administration columns: dose, start, interval, rate, annual.
  # The rr.infect / rr.clinical / rr.hosp columns are consumed downstream in
  # mod-infection.R and mod-progress.R through compute_ve().
  
  n_pop <- sum(active == 1) # active nodes
  
  # vax.supply.rate caps per-timestep doses
  remaining_supply <- if (is.infinite(vax.supply.rate)) {
    Inf # no cap
  } else {
    floor(vax.supply.rate * n_pop) # number of eligible population who'll receive vaccine
  }
  
  #  vax.supply.total caps cumulative first-dose coverage
  n_greater_than_equal_1_dose <- sum(active == 1 & vax >= 1) # number of active people who already got at least 1 dose
  
  remaining_firstdose <- if (is.infinite(vax.supply.total)) { 
    Inf
  } else {  # how many more first doses are allowed 
    floor(vax.supply.total * n_pop) - # total number of population receive vaccines
      n_greater_than_equal_1_dose
  }
  
  remaining_firstdose <- max(0, remaining_firstdose) # number of remaining first dose
  
  ids_newly_vaxed <- integer(0) 
  ids_elig_all <- integer(0)
  

  ######## sequential multiple-dose event (initiated using start and interval) ########
  # Note: no parallel pathways (regular/boost) as in CorporateMix
  daily_vax_by_dose <- rep(0L, nrow(vax.schedule)) # length is the number of doses
  
  for (dose_i in seq_len(nrow(vax.schedule)) # for each (dose_i) dose
       ) {

    dose_row <- vax.schedule[dose_i, , drop = FALSE] 
    
    idsElig <- get_ids_eligible_for_dose(
      dose_row = dose_row, # contains column in vax.schedule
      dose_i = dose_i,
      active = active,
      status = status,
      dxStatus = dxStatus,
      dxTime = dxTime,
      vax = vax,
      last.dose.time = last.dose.time,
      vax.age.group = vax.age.group,
      at = at,
      season.length = season.length
    )

    # Infants are not eligible for the adult vaccine when vax.infant.eligible is
    # FALSE (RSV); they are protected only indirectly (household) or by separate
    # passive products not modelled in this allocation.
    if (!vax.infant.eligible && !is.null(is_infant) && length(idsElig) > 0) {
      idsElig <- idsElig[!is_infant[idsElig]]
    }

    nElig <- length(idsElig)
    
    if (nElig > 0) {
      
    rate_vec <- dose_row$rate[[1]] # Uptake rate across age group per time step
    
    effective_supply <- if (dose_i == 1) {
      min(remaining_supply, remaining_firstdose)
    } else {
      remaining_supply
    }
    
    idsVax <- allocation_strategy(
      idsElig = idsElig,
      rate = rate_vec,
      vax.age.group = vax.age.group,
      vax.strategy = vax.strategy,
      degree_total = degree_total,
      is_bridge = is_bridge,
      n_layers_active = n_layers_active,
      remaining_supply = effective_supply,
      hh.ids = hh.ids,
      is_infant = is_infant,
      proxy_score = proxy_score,
      vax.target.uptake = vax.target.uptake,
      vax.hybrid.frac = vax.hybrid.frac
    )

    ids_elig_all <- c(ids_elig_all, idsElig)
    ids_newly_vaxed <- c(ids_newly_vaxed, idsVax)
    
    nVax <- length(idsVax)
    daily_vax_by_dose[dose_i] <- nVax
    
    if (nVax > 0) {
      vax[idsVax] <- dose_i
      last.dose.time[idsVax] <- at # record the timestep of this (most recent) dose
      remaining_supply <- remaining_supply - nVax
      
      if (dose_i == 1) {
        remaining_firstdose <- remaining_firstdose - nVax
      }
    }
    }
  }
  
  # Replace attr
  dat <- set_attr(dat, "vax", vax)
  dat <- set_attr(dat, "last.dose.time", last.dose.time)
  
  # Summary statistics 
  ## Number of regular doses received at this timestep
  for (dose_i in seq_len(nrow(vax.schedule))) {
    dat <- set_epi(dat, paste0("nVax", dose_i), at, daily_vax_by_dose[dose_i])
  }
  
  
  ## Age-specific vax coverage #TODO: this needs to be turned to loop like the above
  dat <- set_epi(dat, "cov_vax1_0to4", at, length(which(vax.age.group == 1 & vax >= 1)) / length(which(vax.age.group == 1)))
  dat <- set_epi(dat, "cov_vax1_5to17", at, length(which(vax.age.group == 2 & vax >= 1)) / length(which(vax.age.group == 2)))
  dat <- set_epi(dat, "cov_vax1_18to64", at, length(which((vax.age.group == 3 | vax.age.group == 4) & vax >= 1)) / length(which(vax.age.group == 3 | vax.age.group == 4)))
  dat <- set_epi(dat, "cov_vax1_65p", at, length(which(vax.age.group == 5 & vax >= 1)) / length(which(vax.age.group == 5)))
  
  dat <- set_epi(dat, "cov_vax2_0to4", at, length(which(vax.age.group == 1 & vax >= 2)) / length(which(vax.age.group == 1)))
  dat <- set_epi(dat, "cov_vax2_5to17", at, length(which(vax.age.group == 2 & vax >= 2)) / length(which(vax.age.group == 2)))
  dat <- set_epi(dat, "cov_vax2_18to64", at, length(which((vax.age.group == 3 | vax.age.group == 4) & vax >= 2)) / length(which(vax.age.group == 3 | vax.age.group == 4)))
  dat <- set_epi(dat, "cov_vax2_65p", at, length(which(vax.age.group == 5 & vax >= 2)) / length(which(vax.age.group == 5)))
  
  dat <- set_epi(dat, "cov_vax3_0to4", at, length(which(vax.age.group == 1 & vax >= 3)) / length(which(vax.age.group == 1)))
  dat <- set_epi(dat, "cov_vax3_5to17", at, length(which(vax.age.group == 2 & vax >= 3)) / length(which(vax.age.group == 2)))
  dat <- set_epi(dat, "cov_vax3_18to49", at, length(which(vax.age.group == 3 & vax >= 3)) / length(which(vax.age.group == 3)))
  dat <- set_epi(dat, "cov_vax3_50to64", at, length(which(vax.age.group == 4 & vax >= 3)) / length(which(vax.age.group == 4)))
  dat <- set_epi(dat, "cov_vax3_65p", at, length(which(vax.age.group == 5 & vax >= 3)) / length(which(vax.age.group == 5)))

  ## Dose-1 coverage by NETWORK POSITION (fixed degree_quartile / is_bridge attrs).
  ## Verifies the network strategies actually reach high-degree / bridge nodes
  ## (the degree-strategy acceptance criterion, Q4 > Q1; closes issue #21).
  degree_quartile <- get_attr(dat, "degree_quartile")
  if (!is.null(degree_quartile)) {
    for (q in 1:4) dat <- set_epi(dat, paste0("cov_vax1_degQ", q), at,
      length(which(degree_quartile == q & vax >= 1)) / max(1, length(which(degree_quartile == q))))
  }
  if (!is.null(is_bridge)) {
    dat <- set_epi(dat, "cov_vax1_bridge", at,
      length(which(is_bridge == 1 & vax >= 1)) / max(1, length(which(is_bridge == 1))))
    dat <- set_epi(dat, "cov_vax1_nonbridge", at,
      length(which(is_bridge == 0 & vax >= 1)) / max(1, length(which(is_bridge == 0))))
  }

  # mean deg and num of bridge in those ever vaccinated and never vaccinated
  ## individual ever vaccinated and never vaccinated
  vaccinated_ids <- which(active == 1 & vax >= 1)
  unvaccinated_ids <- which(active == 1 & vax == 0)
  ##  mean degree and bridge number in ever vaccinated individuals vs. in those unvaccinated
  mean_deg_vax <- mean(degree_total[vaccinated_ids])
  mean_deg_unvax <- mean(degree_total[unvaccinated_ids])
  mean_n_layers_active_vax <- mean(n_layers_active[vaccinated_ids])
  mean_n_layers_active_unvax <- mean(n_layers_active[unvaccinated_ids])
  
  ## cumulative coverage rate
  cov_cumulative <- length(vaccinated_ids) / sum(active == 1)
  
  dat <- set_epi(dat, "mean_deg_vax", at, mean_deg_vax)
  dat <- set_epi(dat, "mean_deg_unvax", at, mean_deg_unvax)
  dat <- set_epi(dat, "mean_n_layers_active_vax", at, mean_n_layers_active_vax)
  dat <- set_epi(dat, "mean_n_layers_active_unvax", at, mean_n_layers_active_unvax)
  dat <- set_epi(dat, "cov_cumulative", at, cov_cumulative)
  
  # acceptance check, the mean degree of newly vaccinated individuals should be > that of eligible but not selected individuals.
  eligible_not_selected_ids <- setdiff(unique(ids_elig_all), unique(ids_newly_vaxed))
  
  # this tracker pools across all dose levels per timestep. If we ever need per-dose-level verification, we'd need separate trackers.
  mean_deg_newly_vaxed <- if (length(ids_newly_vaxed) > 0) { 
    mean(degree_total[unique(ids_newly_vaxed)])
  } else NA_real_
  mean_deg_elig_not_selected <- if (length(eligible_not_selected_ids) > 0) {
    mean(degree_total[eligible_not_selected_ids])
  } else NA_real_
  mean_n_layers_active_newly_vaxed <- if (length(ids_newly_vaxed) > 0) {
    mean(n_layers_active[unique(ids_newly_vaxed)])
  } else NA_real_
  mean_n_layers_active_elig_not_selected <- if (length(eligible_not_selected_ids) > 0) {
    mean(n_layers_active[eligible_not_selected_ids])
  } else NA_real_
  
  dat <- set_epi(dat, "mean_deg_newly_vaxed", at, mean_deg_newly_vaxed)
  dat <- set_epi(dat, "mean_deg_elig_not_selected", at, mean_deg_elig_not_selected)
  dat <- set_epi(dat, "mean_n_layers_active_newly_vaxed", at, mean_n_layers_active_newly_vaxed)
  dat <- set_epi(dat, "mean_n_layers_active_elig_not_selected", at, mean_n_layers_active_elig_not_selected)
  
  # acceptance check, Verified: with a cap of 0.005, daily vaccination count never exceeds 0.5% of population size
  daily_vax_total <- sum(daily_vax_by_dose)
  
  daily_vax_cap <- if (is.infinite(vax.supply.rate)) {
    Inf
  } else {
    floor(vax.supply.rate * n_pop)
  }
  
  dat <- set_epi(dat, "remaining_supply", at, remaining_supply)
  dat <- set_epi(dat, "remaining_firstdose", at, remaining_firstdose)
  dat <- set_epi(dat, "daily_vax_total", at, daily_vax_total)
  dat <- set_epi(dat, "daily_vax_cap", at, daily_vax_cap)
  
  return(dat)
}

vax_covid_corporate <- function(dat, at) {
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  age <- get_attr(dat, "age")
  vax <- get_attr(dat, "vax")
  vax1Time <- get_attr(dat, "vax1Time")
  vax2Time <- get_attr(dat, "vax2Time")
  vax3Time <- get_attr(dat, "vax3Time")
  vax4Time <- get_attr(dat, "vax4Time")
  vax.age.breaks<- get_param(dat, "vax.age.breaks")
  vax.age.group <-
    cut(
      age,
      breaks = vax.age.breaks,
      right = F,
      labels = 1:5
      
    ) |> as.character() |> as.integer()
    
  dxStatus <- get_attr(dat, "dxStatus")
  dxTime <- get_attr(dat, "dxTime")
  
  
  degree_total <- get_attr(dat, "degree_total") # total nodal degree
  is_bridge <- get_attr(dat, "is_bridge")
  n_layers_active <- get_attr(dat, "n_layers_active")
  
  vax.strategy <- get_param(dat, "vax.strategy") # determine which stategy to use
  vax.supply.rate <- get_param(dat, "vax.supply.rate")
  vax.supply.total <- get_param(dat, "vax.supply.total")
  
  n_pop <- sum(active == 1) # active nodes
  
  # vax.supply.rate caps per-timestep doses
  remaining_supply <- if (is.infinite(vax.supply.rate)) {
    Inf # no cap
  } else {
    floor(vax.supply.rate * n_pop) # number of eligible population who'll receive vaccine
  }
  
  #  vax.supply.total caps cumulative first-dose coverage
  n_greater_than_equal_1_dose <- sum(active == 1 & vax >= 1) # number of active people who already got at least 1 dose
  
  remaining_firstdose <- if (is.infinite(vax.supply.total)) { 
    Inf
  } else {  # how many more first doses are  allowed 
    floor(vax.supply.total * n_pop) - # total number of population receive vaccines
      n_greater_than_equal_1_dose
  }

  remaining_firstdose <- max(0, remaining_firstdose) # number of remaining firstdose
  
  
  vax1.start <- get_param(dat, "vax1.start")
  vax2.interval <- get_param(dat, "vax2.interval")
  vax3.start <- get_param(dat, "vax3.start")
  vax3.interval <- get_param(dat, "vax3.interval")
  vax4.start <- get_param(dat, "vax4.start")
  vax4.interval <- get_param(dat, "vax4.interval")
  
  vax1.rate <- get_param(dat, "vax1.rate")
  vax2.rate <- get_param(dat, "vax2.rate")
  vax3.rate <- get_param(dat, "vax3.rate")
  vax4.rate <- get_param(dat, "vax4.rate")
  
  vax1.boost <- get_param(dat, "vax1.boost")
  vax2.boost <- get_param(dat, "vax2.boost")
  vax3.boost <- get_param(dat, "vax3.boost")
  vax4.boost <- get_param(dat, "vax4.boost")
  vax1.boost.start <- get_param(dat, "vax1.boost.start")
  vax2.boost.start <- get_param(dat, "vax2.boost.start")
  vax3.boost.start <- get_param(dat, "vax3.boost.start")
  vax4.boost.start <- get_param(dat, "vax4.boost.start")
  
  # initialize these 4 strings for checking the number of people received boosters
  ids.vax1.boost <- integer(0)
  ids.vax2.boost <- integer(0)
  ids.vax3.boost <- integer(0)
  ids.vax4.boost <- integer(0)
  
  ids_newly_vaxed <- integer(0)
  ids_elig_all <- integer(0)
  
  if (any(at == vax1.boost.start)) {
    idsElig.vax1.boost <- which(active == 1 & !(status %in% c("ic", "h"))
                                & !(dxStatus == 2 & (at - dxTime <= 10)) & vax == 0
                                & at == vax1.boost.start[vax.age.group])
    nElig.vax1.boost <- length(idsElig.vax1.boost)
    
    if (nElig.vax1.boost > 0) {
      
      effective_supply <- min(remaining_supply, remaining_firstdose) # when remaining_firstdose = Inf, use remaining supply
      
      ids.vax1.boost <- allocation_strategy(
        idsElig = idsElig.vax1.boost,
        rate = vax1.boost,
        vax.age.group = vax.age.group,
        vax.strategy = vax.strategy,
        degree_total = degree_total,
        is_bridge = is_bridge,
        n_layers_active = n_layers_active,
        remaining_supply = effective_supply
      )
      
      ids_elig_all <- c(ids_elig_all, idsElig.vax1.boost)     
      ids_newly_vaxed <- c(ids_newly_vaxed, ids.vax1.boost) 
      
      nVax1.boost <- length(ids.vax1.boost)
      
      if (nVax1.boost > 0) {
        vax[ids.vax1.boost] <- 1
        vax1Time[ids.vax1.boost] <- at
        remaining_supply <- remaining_supply - nVax1.boost
        remaining_firstdose <- remaining_firstdose - nVax1.boost
      }
    }
  }
  
  if (any(at == vax2.boost.start)) {
    idsElig.vax2.boost <- which(active == 1 & !(status %in% c("ic", "h"))
                                & !(dxStatus == 2 & (at - dxTime <= 10)) & vax == 1
                                & at == vax2.boost.start[vax.age.group])
    nElig.vax2.boost <- length(idsElig.vax2.boost)
    
    if (nElig.vax2.boost > 0) {

      ids.vax2.boost <- allocation_strategy(
        idsElig = idsElig.vax2.boost,
        rate = vax2.boost,
        vax.age.group = vax.age.group,
        vax.strategy = vax.strategy,
        degree_total = degree_total,
        is_bridge = is_bridge,
        n_layers_active = n_layers_active,
        remaining_supply = remaining_supply
      )
      
      ids_elig_all <- c(ids_elig_all, idsElig.vax2.boost)     
      ids_newly_vaxed <- c(ids_newly_vaxed, ids.vax2.boost) 
      
      nVax2.boost <- length(ids.vax2.boost)
      
      if (nVax2.boost > 0) {
      
        vax[ids.vax2.boost] <- 2
        vax2Time[ids.vax2.boost] <- at
        remaining_supply <- remaining_supply - nVax2.boost
        
      }
    }
  }
  
  if (any(at == vax3.boost.start)) {
    idsElig.vax3.boost <- which(active == 1 & !(status %in% c("ic", "h"))
                                & !(dxStatus == 2 & (at - dxTime <= 10)) & vax == 2
                                & at == vax3.boost.start[vax.age.group])
    nElig.vax3.boost <- length(idsElig.vax3.boost)
    
    if (nElig.vax3.boost > 0) {
      
      ids.vax3.boost <- allocation_strategy(
        idsElig = idsElig.vax3.boost,
        rate = vax3.boost,
        vax.age.group = vax.age.group,
        vax.strategy = vax.strategy,
        degree_total = degree_total,
        is_bridge = is_bridge,
        n_layers_active = n_layers_active,
        remaining_supply = remaining_supply
      )
      
      ids_elig_all <- c(ids_elig_all, idsElig.vax3.boost)     
      ids_newly_vaxed <- c(ids_newly_vaxed, ids.vax3.boost) 
      
      nVax3.boost <- length(ids.vax3.boost)
      
      if (nVax3.boost > 0) {
  
        vax[ids.vax3.boost] <- 3
        vax3Time[ids.vax3.boost] <- at
        remaining_supply <- remaining_supply - nVax3.boost
        
      }
    }
  }
  
  if (any(at == vax4.boost.start)) {
    idsElig.vax4.boost <- which(active == 1 & !(status %in% c("ic", "h"))
                                & !(dxStatus == 2 & (at - dxTime <= 10)) & vax == 3
                                & at == vax4.boost.start[vax.age.group])
    nElig.vax4.boost <- length(idsElig.vax4.boost)
    
    if (nElig.vax4.boost > 0) {
      
      ids.vax4.boost <- allocation_strategy(
        idsElig = idsElig.vax4.boost,
        rate = vax4.boost,
        vax.age.group = vax.age.group,
        vax.strategy = vax.strategy,
        degree_total = degree_total,
        is_bridge = is_bridge,
        n_layers_active = n_layers_active,
        remaining_supply = remaining_supply
      )
      
      ids_elig_all <- c(ids_elig_all, idsElig.vax4.boost)     
      ids_newly_vaxed <- c(ids_newly_vaxed, ids.vax4.boost) 
      
      nVax4.boost <- length(ids.vax4.boost)
      
      if (nVax4.boost > 0) {
      
        vax[ids.vax4.boost] <- 4
        vax4Time[ids.vax4.boost] <- at
        remaining_supply <- remaining_supply - nVax4.boost
        
      }
    }
  }
  
  ## First vax
  nVax1 <- 0
  idsElig.vax1 <- which(active == 1 & !(status %in% c("ic", "h"))
                        & !(dxStatus == 2 & (at - dxTime <= 10)) & vax == 0
                        & at >= vax1.start[vax.age.group]) # acceptance: Respects existing age-eligibility time windows

  nElig.vax1 <- length(idsElig.vax1)
  if (nElig.vax1 > 0) {
    effective_supply <- min(remaining_supply, remaining_firstdose)
  
    idsVax1 <- allocation_strategy(
      idsElig = idsElig.vax1,
      rate = vax1.rate,
      vax.age.group = vax.age.group,
      vax.strategy = vax.strategy,
      degree_total = degree_total,
      is_bridge = is_bridge,
      n_layers_active = n_layers_active,
      remaining_supply = effective_supply
    )
    
    ids_elig_all <- c(ids_elig_all, idsElig.vax1)     
    ids_newly_vaxed <- c(ids_newly_vaxed, idsVax1) 
    
    nVax1 <- length(idsVax1)
    
    if (nVax1>0) {
      vax[idsVax1] <- 1
      vax1Time[idsVax1] <- at
      remaining_supply <- remaining_supply - length(idsVax1)
      remaining_firstdose <- remaining_firstdose - nVax1
    }
  }
  
  ## Second vax
  nVax2 <- 0
  idsElig.vax2 <- which(active == 1 & !(status %in% c("ic", "h"))
                        & !(dxStatus == 2 & (at - dxTime <= 10)) & vax == 1
                        & (at - vax1Time >= vax2.interval))
  nElig.vax2 <- length(idsElig.vax2)
  if (nElig.vax2 > 0) {

    idsVax2 <- allocation_strategy(
      idsElig = idsElig.vax2,
      rate = vax2.rate,
      vax.age.group = vax.age.group,
      vax.strategy = vax.strategy,
      degree_total = degree_total,
      is_bridge = is_bridge,
      n_layers_active = n_layers_active,
      remaining_supply = remaining_supply
    )
    
    ids_elig_all <- c(ids_elig_all, idsElig.vax2)     
    ids_newly_vaxed <- c(ids_newly_vaxed, idsVax2) 
    
    nVax2 <- length(idsVax2)
    
    if (nVax2>0) {
      vax[idsVax2] <- 2
      vax2Time[idsVax2] <- at
      remaining_supply <- remaining_supply - nVax2
    }
    
  }
  
  ## Third vax 
  nVax3 <- 0
  idsElig.vax3 <- which(active == 1 & !(status %in% c("ic", "h"))
                        & !(dxStatus == 2 & (at - dxTime <= 10)) & vax == 2
                        & (at - vax2Time >= vax3.interval)
                        & at >= vax3.start[vax.age.group])
  nElig.vax3 <- length(idsElig.vax3)
  if (nElig.vax3 > 0) {
    
    idsVax3 <- allocation_strategy(
      idsElig = idsElig.vax3,
      rate = vax3.rate,
      vax.age.group = vax.age.group,
      vax.strategy = vax.strategy,
      degree_total = degree_total,
      is_bridge = is_bridge,
      n_layers_active = n_layers_active,
      remaining_supply = remaining_supply
    )
    
    ids_elig_all <- c(ids_elig_all, idsElig.vax3)     
    ids_newly_vaxed <- c(ids_newly_vaxed, idsVax3) 
    
    nVax3 <- length(idsVax3)
    
    if (nVax3>0) {
      vax[idsVax3] <- 3
      vax3Time[idsVax3] <- at
      remaining_supply <- remaining_supply - nVax3
    }
    
  }
  
  
  ## Fourth vax 
  nVax4 <- 0
  idsElig.vax4 <- which(active == 1 & !(status %in% c("ic", "h"))
                        & !(dxStatus == 2 & (at - dxTime <= 10)) & vax == 3
                        & (at - vax3Time >= vax4.interval)
                        & at >= vax4.start[vax.age.group])
  nElig.vax4 <- length(idsElig.vax4)
  if (nElig.vax4 > 0) {
    
    idsVax4 <- allocation_strategy(
      idsElig = idsElig.vax4,
      rate = vax4.rate,
      vax.age.group = vax.age.group,
      vax.strategy = vax.strategy,
      degree_total = degree_total,
      is_bridge = is_bridge,
      n_layers_active = n_layers_active,
      remaining_supply = remaining_supply
    )
    
    ids_elig_all <- c(ids_elig_all, idsElig.vax4)     
    ids_newly_vaxed <- c(ids_newly_vaxed, idsVax4) 
    
    nVax4 <- length(idsVax4)
    
    if (nVax4 > 0) {
      vax[idsVax4] <- 4
      vax4Time[idsVax4] <- at
      remaining_supply <- remaining_supply - nVax4
    }
    
  }
  
  ## Replace attr
  dat <- set_attr(dat, "vax", vax)
  dat <- set_attr(dat, "vax1Time", vax1Time)
  dat <- set_attr(dat, "vax2Time", vax2Time)
  dat <- set_attr(dat, "vax3Time", vax3Time)
  dat <- set_attr(dat, "vax4Time", vax4Time)
  
  ## Summary statistics ##
  dat <- set_epi(dat, "nVax1", at, nVax1)
  dat <- set_epi(dat, "nVax2", at, nVax2)
  dat <- set_epi(dat, "nVax3", at, nVax3)
  dat <- set_epi(dat, "nVax4", at, nVax4)
  
  ## Summary stats -- vaccine coverage ##
  dat <- set_epi(dat, "cov_vax1_0to4", at, length(which(vax.age.group == 1 & vax >= 1)) / length(which(vax.age.group == 1)))
  dat <- set_epi(dat, "cov_vax1_5to17", at, length(which(vax.age.group == 2 & vax >= 1)) / length(which(vax.age.group == 2)))
  dat <- set_epi(dat, "cov_vax1_18to64", at, length(which((vax.age.group == 3 | vax.age.group == 4) & vax >= 1)) / length(which(vax.age.group == 3 | vax.age.group == 4)))
  dat <- set_epi(dat, "cov_vax1_65p", at, length(which(vax.age.group == 5 & vax >= 1)) / length(which(vax.age.group == 5)))

  dat <- set_epi(dat, "cov_vax2_0to4", at, length(which(vax.age.group == 1 & vax >= 2)) / length(which(vax.age.group == 1)))
  dat <- set_epi(dat, "cov_vax2_5to17", at, length(which(vax.age.group == 2 & vax >= 2)) / length(which(vax.age.group == 2)))
  dat <- set_epi(dat, "cov_vax2_18to64", at, length(which((vax.age.group == 3 | vax.age.group == 4) & vax >= 2)) / length(which(vax.age.group == 3 | vax.age.group == 4)))
  dat <- set_epi(dat, "cov_vax2_65p", at, length(which(vax.age.group == 5 & vax >= 2)) / length(which(vax.age.group == 5)))

  dat <- set_epi(dat, "cov_vax3_5to17", at, length(which(vax.age.group == 2 & vax >= 3)) / length(which(vax.age.group == 2)))
  dat <- set_epi(dat, "cov_vax3_18to49", at, length(which(vax.age.group == 3 & vax >= 3)) / length(which(vax.age.group == 3)))
  dat <- set_epi(dat, "cov_vax3_50to64", at, length(which(vax.age.group == 4 & vax >= 3)) / length(which(vax.age.group == 4)))
  dat <- set_epi(dat, "cov_vax3_65p", at, length(which(vax.age.group == 5 & vax >= 3)) / length(which(vax.age.group == 5)))

  dat <- set_epi(dat, "cov_vax4_50to64", at, length(which(vax.age.group == 4 & vax >= 4)) / length(which(vax.age.group == 4)))
  dat <- set_epi(dat, "cov_vax4_65p", at, length(which(vax.age.group == 5 & vax >= 4)) / length(which(vax.age.group == 5)))
  
  # mean deg and num of bridge in those ever vaccinated and never vaccinated
  ## individual ever vaccinated and never vaccinated
  vaccinated_ids <- which(active == 1 & vax >= 1)
  unvaccinated_ids <- which(active == 1 & vax == 0)
  ##  mean degree and bridge number in ever vaccinated individuals vs. in those unvaccinated
  mean_deg_vax <- mean(degree_total[vaccinated_ids])
  mean_deg_unvax <- mean(degree_total[unvaccinated_ids])
  mean_n_layers_active_vax <- mean(n_layers_active[vaccinated_ids])
  mean_n_layers_active_unvax <- mean(n_layers_active[unvaccinated_ids])
  
  ## cumulative coverage rate
  cov_cumulative <- length(vaccinated_ids)/sum(active == 1)
  
  dat <- set_epi(dat, "mean_deg_vax", at, mean_deg_vax)
  dat <- set_epi(dat, "mean_deg_unvax", at, mean_deg_unvax)
  dat <- set_epi(dat, "mean_n_layers_active_vax", at, mean_n_layers_active_vax)
  dat <- set_epi(dat, "mean_n_layers_active_unvax", at, mean_n_layers_active_unvax)
  dat <- set_epi(dat, "cov_cumulative", at, cov_cumulative)

  # acceptance check, the mean degree of newly vaccinated individuals should be > that of eligible but not selected individuals.
  eligible_not_selected_ids <- setdiff(unique(ids_elig_all), unique(ids_newly_vaxed))
  
  # this tracker pools across all dose levels per timestep. If we ever need per-dose-level verification, we'd need separate trackers.
  mean_deg_newly_vaxed <- if (length(ids_newly_vaxed) > 0) { 
    mean(degree_total[unique(ids_newly_vaxed)])
  } else NA_real_
  mean_deg_elig_not_selected <- if (length(eligible_not_selected_ids) > 0) {
    mean(degree_total[eligible_not_selected_ids])
  } else NA_real_
  mean_n_layers_active_newly_vaxed <- if (length(ids_newly_vaxed) > 0) {
      mean(n_layers_active[unique(ids_newly_vaxed)])
    } else NA_real_
  mean_n_layers_active_elig_not_selected <- if (length(eligible_not_selected_ids) > 0) {
      mean(n_layers_active[eligible_not_selected_ids])
    } else NA_real_
  
  dat <- set_epi(dat, "mean_deg_newly_vaxed", at, mean_deg_newly_vaxed)
  dat <- set_epi(dat, "mean_deg_elig_not_selected", at, mean_deg_elig_not_selected)
  dat <- set_epi(dat, "mean_n_layers_active_newly_vaxed", at, mean_n_layers_active_newly_vaxed)
  dat <- set_epi(dat, "mean_n_layers_active_elig_not_selected", at, mean_n_layers_active_elig_not_selected)

  
  # acceptance check, Verified: with a cap of 0.005, daily vaccination count never exceeds 0.5% of population size
  daily_vax_total <- nVax1 + nVax2 + nVax3 + nVax4 +
    length(ids.vax1.boost) + length(ids.vax2.boost) +
    length(ids.vax3.boost) + length(ids.vax4.boost)
  
  
  daily_vax_cap <- if (is.infinite(vax.supply.rate)) {
    Inf
  } else {
    floor(vax.supply.rate * n_pop)
  }
  
  dat <- set_epi(dat, "remaining_supply", at, remaining_supply)
  dat <- set_epi(dat, "remaining_firstdose", at, remaining_firstdose)
  dat <- set_epi(dat, "daily_vax_total", at, daily_vax_total)
  dat <- set_epi(dat, "daily_vax_cap", at, daily_vax_cap)
  
  return(dat)
}

# vaccine allocation strategy
allocation_strategy <- function(idsElig, rate, vax.age.group, vax.strategy,
                                degree_total, # argument for degree strategy
                                is_bridge, n_layers_active,
                                remaining_supply,
                                hh.ids = NULL, is_infant = NULL,
                                proxy_score = NULL, vax.target.uptake = 1,
                                vax.hybrid.frac = 0) {
                                # hh.ids + is_infant default NULL so existing call
                                # sites (vax_covid_corporate) are unaffected; only
                                # the "infant_hh" parent-targeting branch needs them.
                                # proxy_score (issue #40) is the field-observable
                                # ranking score for the proxy_* strategies, precomputed
                                # by the caller; vax.target.uptake (<1) degrades
                                # hub-targeting precision. Both default to the
                                # unchanged behaviour.
  nElig <- length(idsElig) # number of eligible people 
  
  if (nElig == 0) { # if nobody is eligible, return an empty vector, meaning nobody is eligible and exit the function
    return(integer(0)) 
  }
  
  if (remaining_supply <= 0) { # supply mechanism: if no vaccine doses are left for this timestep-nobody would be vaccinated and exit the function
    return(integer(0))
  }

  # First step: priority mechanism-rank eligible people by prioritization strategy
  if (vax.strategy == "degree") { 
    
    # First, degree-based allocation - rank all eligible individuals by degree score from high to low
    degree_jitter <- degree_total[idsElig] + runif(nElig, 0, 1e-8) # add small random noise breaks ties randomly
    
    names(degree_jitter) <- idsElig  # attach eligible ids as names so they stay linked to their degree_jitters
    
    ids_ranked <- as.integer(names(sort(degree_jitter, decreasing = TRUE))) # sort jittered degree from high to low, then recover the ranked IDs
     
  } else if (vax.strategy == "bridge") {
    
    max_deg <- max(degree_total) + 1 # define a multiplier so n_layers_active has higher priority than degree_total
    
    # First, bridge-based allocation - rank all eligible individuals by bridge score from high to low
    score <- as.integer(is_bridge[idsElig]) * ((max(n_layers_active) +1) * max_deg) + # bridge gets the biggest weight, (max(n_layers_active) +1)=5 in the 4-layer model
      n_layers_active[idsElig] * max_deg + # then n_layers_active (i.e., 5> max(n_layers_active))
      degree_total[idsElig] +  # then degree_total (i.e., max_deg> degree_total)
      runif(nElig, 0, 1e-8) # tiny random noise to break ties
    
    names(score) <- idsElig  # attach IDs of who are eligible as names
    
    ids_ranked <- as.integer(names(sort(score, decreasing = TRUE)))  # sort from high to low and get ranked IDs
    
    
  }   else if (vax.strategy == "age") {
    # First, rank all eligible individuals by age from high to low
    age_priority <- rate[vax.age.group[idsElig]] # among those who eligible, use the age-specific vaccination rate as the priority ordering
    
    age_priority_jitter <- age_priority + runif(nElig, 0, 1e-8) # higher rate = higher priority; add random noise to break ties randomly
    
    names(age_priority_jitter) <- idsElig # attach IDs as names so they stay linked to their ordering
    
    ids_ranked <- as.integer(names(sort(age_priority_jitter, decreasing = TRUE)))# sort from high to low and get the ranked IDs
    
  } else if (vax.strategy == "random") {
    # First uniformly rank eligible people at random
    ids_ranked <- sample(idsElig, nElig)

  } else if (vax.strategy == "age_then_degree") {
    # Two-stage hybrid: vaccinate the oldest first (vax age groups 4-5, ~50+) to
    # protect the vulnerable, then give any remaining supply to the highest-degree
    # of everyone else to cut transmission. Aims to capture BOTH the mortality
    # benefit of age-targeting and the infection benefit of network-targeting,
    # which no single pure strategy can (the corners of the supply/objective frontier).
    vag <- vax.age.group[idsElig]
    is_old <- vag >= 4
    old <- idsElig[is_old]; old_age <- vag[is_old]
    rest <- idsElig[!is_old]; rest_deg <- degree_total[rest]
    old_ord  <- old[order(-old_age, runif(length(old)))]                          # 65+ before 50-64
    rest_ord <- rest[order(-(rest_deg + runif(length(rest), 0, 1e-8)))]           # then by degree
    ids_ranked <- c(old_ord, rest_ord)

  } else if (vax.strategy == "infant_hh") {
    # Parent-targeting (issue #23, RSV): infants are network-peripheral (household
    # only) and cannot be reached by degree/bridge targeting, but the ADULTS who
    # share their household can. Rank eligible co-residents of an infant household
    # first (the "parents"), by their age-specific rate, then everyone else by the
    # same rate. Household membership comes from hh.ids; infant households are
    # those containing ANY infant (infants themselves are not vax-eligible).
    if (is.null(hh.ids) || is.null(is_infant)) {
      stop("vax.strategy 'infant_hh' requires hh.ids and is_infant")
    }
    hh_with_infant <- unique(hh.ids[which(is_infant)])
    in_infant_hh   <- hh.ids[idsElig] %in% hh_with_infant
    age_priority   <- rate[vax.age.group[idsElig]] + runif(nElig, 0, 1e-8)
    parents <- idsElig[in_infant_hh]
    others  <- idsElig[!in_infant_hh]
    parents_ord <- parents[order(-age_priority[in_infant_hh])]
    others_ord  <- others[order(-age_priority[!in_infant_hh])]
    ids_ranked  <- c(parents_ord, others_ord)

  } else if (vax.strategy %in% c("proxy_hh", "proxy_fit")) {
    # Field-deployable proxy (issue #40): rank by a score built ONLY from cheaply
    # observable covariates (household size for proxy_hh; a fitted predicted-degree
    # composite of household size, school enrollment and employment for proxy_fit),
    # NOT by the true, field-unobservable contact degree. proxy_score is precomputed
    # by the caller from those observables. Measures how much of degree-targeting's
    # gain a deployable proxy recovers.
    if (is.null(proxy_score)) {
      stop("vax.strategy '", vax.strategy, "' requires proxy_score")
    }
    score <- proxy_score[idsElig] + runif(nElig, 0, 1e-8) # tiny noise breaks ties
    names(score) <- idsElig
    ids_ranked <- as.integer(names(sort(score, decreasing = TRUE)))

  } else if (vax.strategy %in% c("age_plus_proxy", "age_plus_degree")) {
    # Parallel multi-criteria hybrid (issue #41): reserve a fraction `f`
    # (vax.hybrid.frac) of each timestep's doses for the 65+ elderly and allocate
    # the remaining 1-f by a transmission criterion, the field-deployable proxy
    # (age_plus_proxy) or the idealized true degree (age_plus_degree ceiling).
    # Unlike the sequential age_then_degree, which fills all 50+ first and, because
    # that band is large relative to 5-20% supply, exhausts the budget before the
    # network stage, this splits EVERY step's supply so both criteria get doses.
    # The split is at allocation, so this branch returns its doses directly.
    if (vax.strategy == "age_plus_proxy" && is.null(proxy_score)) {
      stop("vax.strategy 'age_plus_proxy' requires proxy_score")
    }
    # Stochastic daily demand throttled by supply (same model as the network arms).
    rate_person <- rate[vax.age.group[idsElig]]
    n_willing <- sum(rbinom(nElig, size = 1, prob = rate_person))
    if (n_willing == 0) return(integer(0))
    n_take <- min(remaining_supply, n_willing)
    n_age  <- floor(vax.hybrid.frac * n_take)

    # Age reserve: 65+ (vax age band 5), random order within the band (rate + jitter).
    age_pool   <- idsElig[vax.age.group[idsElig] == 5]
    age_ranked <- age_pool[order(-(rate[vax.age.group[age_pool]] +
                                     runif(length(age_pool), 0, 1e-8)))]
    ids_age    <- utils::head(age_ranked, n_age)

    # Transmission remainder: everyone NOT taken for the reserve, ranked by the
    # field proxy (age_plus_proxy) or true degree (age_plus_degree). n_proxy takes
    # ALL leftover doses, so a short 65+ pool spills into the remainder and the
    # total handed out is always min(n_take, nElig).
    n_proxy    <- n_take - length(ids_age)
    trans_all  <- if (vax.strategy == "age_plus_proxy") proxy_score else degree_total
    proxy_pool <- setdiff(idsElig, ids_age)
    proxy_rank <- proxy_pool[order(-(trans_all[proxy_pool] +
                                      runif(length(proxy_pool), 0, 1e-8)))]
    ids_proxy  <- utils::head(proxy_rank, n_proxy)

    return(c(ids_age, ids_proxy))

  } else if (vax.strategy == "infant_direct") {
    # Direct infant product (issue #41, RSV): a nirsevimab-like monoclonal or
    # maternal-derived protection applied to the infant itself, so it requires
    # vax.infant.eligible = TRUE on these rows (set per-scenario in the grid) to
    # keep infants in the eligible pool. Only infants receive it, so allocation is
    # CAPPED at the infant pool: no spillover to adults even when the nominal
    # supply cap exceeds the number of infants (infants are ~1% of the population,
    # supply is 5-20%). Infants are network-peripheral, so this protects the infant
    # directly with little herd effect; the modelled value is a clean direct-vs-
    # cocooning comparison, not a transmission change.
    if (is.null(is_infant)) {
      stop("vax.strategy 'infant_direct' requires is_infant")
    }
    inf_ids <- idsElig[is_infant[idsElig]]
    if (length(inf_ids) == 0) return(integer(0))
    rate_inf <- rate[vax.age.group[inf_ids]]
    n_willing <- sum(rbinom(length(inf_ids), size = 1, prob = rate_inf))
    if (n_willing == 0) return(integer(0))
    n_take <- min(remaining_supply, n_willing, length(inf_ids))
    return(inf_ids[sample.int(length(inf_ids), n_take)]) # infants are equivalent targets

  } else if (vax.strategy == "infant_direct_cocoon") {
    # Direct product PLUS household cocooning (issue #41): dose infants first
    # (the direct product), then their eligible adult co-residents (the infant_hh
    # cocooning logic), then everyone else. Tests cocooning as a public complement
    # to a private monoclonal. Falls through to the network-style cap, so once
    # infants and infant-household adults are covered any remaining supply spills
    # to the general population exactly as infant_hh does, keeping adult coverage
    # comparable at matched supply.
    if (is.null(hh.ids) || is.null(is_infant)) {
      stop("vax.strategy 'infant_direct_cocoon' requires hh.ids and is_infant")
    }
    hh_with_infant <- unique(hh.ids[which(is_infant)])
    is_inf_elig    <- is_infant[idsElig]
    in_infant_hh   <- hh.ids[idsElig] %in% hh_with_infant
    age_priority   <- rate[vax.age.group[idsElig]] + runif(nElig, 0, 1e-8)
    # tier 2 = infants themselves; tier 1 = adult co-residents of infants; tier 0 = rest
    tier <- ifelse(is_inf_elig, 2L, ifelse(in_infant_hh, 1L, 0L))
    ids_ranked <- idsElig[order(-tier, -age_priority)]

  } else {
    stop("Unknown vax.strategy: ", vax.strategy)
  }
  
  # Second step: determine vaccination allocation
  ## Get each ranked eligible person's age-specific vax probability.
  rate_person <- rate[vax.age.group[ids_ranked]]

  if (vax.strategy %in% c("degree", "bridge", "infant_hh", "proxy_hh", "proxy_fit",
                          "infant_direct_cocoon")) {
    ## For network-based priority strategies: stochastic demand (driven by
    ## age-specific rates) determines how many get vaccinated each timestep;
    ## the priority ranking determines who.  This ensures high-priority
    ## individuals are selected even when demand < supply.
    n_willing <- sum(rbinom(nElig, size = 1, prob = rate_person))
    if (n_willing == 0) return(integer(0))
    n_take <- min(remaining_supply, n_willing)
    if (vax.target.uptake < 1 &&
        vax.strategy %in% c("degree", "bridge", "proxy_hh", "proxy_fit")) {
      ## Imperfect hub-targeting (degraded-uptake arm): only a fraction
      ## (vax.target.uptake) of the doses reach the intended top-ranked targets;
      ## the rest are mis-allocated to random eligibles. Total doses (coverage)
      ## are unchanged, so this isolates targeting PRECISION, interpolating from
      ## the pure network strategy (uptake = 1) toward random (uptake = 0). age,
      ## random and infant_hh are excluded: their targets are trivially
      ## identifiable in the field.
      n_targeted <- round(vax.target.uptake * n_take)
      ids_top <- ids_ranked[seq_len(n_targeted)]
      pool <- setdiff(idsElig, ids_top)
      n_rand <- n_take - n_targeted
      ids_rand <- if (n_rand > 0 && length(pool) > 0) {
        pool[sample.int(length(pool), min(n_rand, length(pool)))]
      } else integer(0)
      ids.vax <- c(ids_top, ids_rand)
    } else {
      ids.vax <- ids_ranked[seq_len(n_take)]
    }
  } else {
    ## For random/age strategies: individual show-up filter then cap by supply
    show_up <- rbinom(nElig, size = 1, prob = rate_person) == 1
    idsShow <- ids_ranked[show_up]
    nShow <- length(idsShow)
    if (nShow == 0) return(integer(0))
    n_take <- min(remaining_supply, nShow)
    ids.vax <- idsShow[1:n_take]
  }

  return(ids.vax)
}

# helper getting generally eligible for a dose
get_ids_eligible_for_dose <- function(
    dose_row, dose_i,
    active, status, dxStatus, dxTime,
    vax,
    last.dose.time = NULL,
    vax.age.group, at,
    season.length = 364
) {
  # Generally eligible: active, not severe (ic/h), not recently diagnosed.
  base_eligible <- active == 1 &
    !(status %in% c("ic", "h")) &
    !(dxStatus == 2 & (at - dxTime <= 10))

  # Age-specific start time for this dose (length 5; indexed by vax.age.group).
  start_vec <- dose_row$start[[1]]
  start_age <- start_vec[vax.age.group]

  # Minimum interval since the previous dose (NA for dose 1).
  interval <- dose_row$interval

  annual <- isTRUE(dose_row$annual)

  if (annual) {
    # Annual (seasonal) revaccination: eligible once per season, reopening each
    # year. Seasons recur every season.length days from the age-specific start;
    # a node is eligible if the current season's campaign has opened and it has
    # not been dosed since it opened (so prior-season vaccinees re-up).
    seasons_elapsed <- pmax(0, floor((at - start_age) / season.length))
    season_start <- start_age + seasons_elapsed * season.length
    idsElig <- which(
      base_eligible &
        at >= season_start &
        (is.na(last.dose.time) | last.dose.time < season_start)
    )

  } else if (dose_i == 1) {
    idsElig <- which(
      base_eligible &
        vax == 0 &
        at >= start_age
    )

  } else {
    # Dose 2+: needs exactly the previous dose, the interval elapsed since it
    # (last.dose.time is the previous dose for a node with vax == dose_i - 1),
    # and the age-specific start time reached.
    idsElig <- which(
      base_eligible &
        vax == dose_i - 1 &
        !is.na(last.dose.time) &
        (at - last.dose.time >= interval) &
        at >= start_age
    )
  }

  return(idsElig)
}

# Compute current vaccine relative risk for each node.
# This helper is shared by mod-infection.R and mod-progress.R despite saved in mod-vax.R
#
# rr_col can be:
#   "rr.infect"   for susceptibility / transmission risk
#   "rr.clinical" for symptomatic disease risk
#   "rr.hosp"     for hospitalization risk
#' @rdname moduleset-gmc19
#' @export
compute_ve <- function(at, ids, vax,
                       vax.age.group,
                       last.dose.time = NULL,
                       vax.schedule,
                       outcome = c("infect", "clinical", "hosp", "death")) {

  outcome <- match.arg(outcome)

  ve_peak_col <- paste0("ve.peak.", outcome)
  ve_halflife_col <- paste0("ve.halflife.", outcome)
  ve_floor_col <- paste0("ve.floor.", outcome)

  # Start everyone at no vaccine effect (rr = 1).
  rr <- rep(1, length(ids))
  ve <- rep(0, length(ids))
  time.since.last.dose <- rep(0, length(ids))

  vax_ids <- vax[ids]

  # A node's protection is set by its most recent dose: last.dose.time is the
  # timestep of that dose and vax (the dose count) selects which dose row's VE
  # profile applies. Supports any number of doses and annual revaccination.
  for (dose_i in seq_len(nrow(vax.schedule))) {

    ids_this_dose_position <- which(vax_ids == dose_i)

    if (length(ids_this_dose_position) > 0) {

      ids_this_dose <- ids[ids_this_dose_position]

      time_since_dose <- at - last.dose.time[ids_this_dose]
      if (any(is.na(time_since_dose))) {
        stop("Missing dose time for vaccinated individuals in compute_ve().")
      }

      age_group_this <- vax.age.group[ids_this_dose]

      ve_peak <- get_schedule_value(vax.schedule, ve_peak_col, dose_i, age_group_this)
      ve_halflife <- get_schedule_value(vax.schedule, ve_halflife_col, dose_i, age_group_this)
      ve_floor <- get_schedule_value(vax.schedule, ve_floor_col, dose_i, age_group_this)
      ve_delay <- get_schedule_value(vax.schedule, "ve.delay", dose_i, age_group_this)

      # VE holds near peak until ve_delay days post-dose, then decays toward the
      # floor by the half-life.
      t_eff <- pmax(0, time_since_dose - ve_delay)
      current_ve <- ve_floor + (ve_peak - ve_floor) * 0.5^(t_eff / ve_halflife)

      if (any(current_ve < 0 | current_ve > 1)) {
        warning("VE outside [0, 1] in compute_ve(). Check ve.peak, ve.floor, and ve.halflife in vax.schedule.")
      }

      ve[ids_this_dose_position] <- current_ve
      rr[ids_this_dose_position] <- 1 - current_ve
      time.since.last.dose[ids_this_dose_position] <- time_since_dose
    }
  }

  return(list(
    rr = rr,
    time.since.last.dose = time.since.last.dose
  ))
}

# Helper to pull either scalar or age-specific schedule values
get_schedule_value <- function(schedule, col, dose_i, age_group) {
  
  # If the schedule column is a list-column, each row may contain
  # an age-specific vector, e.g., c(0.50, 0.50, 0.50, 0.50, 0.60).
  if (is.list(schedule[[col]])) {
    x <- schedule[[col]][[dose_i]]
    
  } else {
    # Otherwise, the column is a regular vector with one scalar value per dose.
    x <- schedule[[col]][dose_i]
  }
  
  # If x is scalar, use the same value for everyone in ids.
  if (length(x) == 1L) {
    rep(x, length(age_group))
    
  } else {
    # If x is age-specific, select the value matching each person's age group.
    x[age_group]
  }
}

# Helper to assign vaccine age groups using vax.age.breaks
#' @rdname moduleset-gmc19
#' @export
vax_age_group_for <- function(dat) {
  age <- get_attr(dat, "age")
  vax.age.breaks <- get_param(dat, "vax.age.breaks")
  
  vax.age.group <- 
    cut(
      age,
      breaks = vax.age.breaks,
      right = FALSE,
      labels = 1:5
    ) |> as.character() |> as.integer()
  
  return(vax.age.group)
}

# Resolve the active pathogen for disease-specific vaccine-schedule selection.
# Defaults to "covid" when the `pathogen` parameter is absent, preserving
# backward compatibility for models that predate multi-pathogen support.
get_pathogen <- function(dat) {
  pathogen <- get_param(dat, "pathogen", override.null.error = TRUE)
  if (is.null(pathogen) || length(pathogen) == 0) "covid" else pathogen
}

# Build dose-indexed vaccine schedule from model_parameters.csv
build_vax_schedule <- function(dat, disease = "covid") {
  
  prefix <- paste0("vax.", disease, ".")
  n_doses <- get_param(dat, paste0(prefix, "n.doses"))
  
  get_dose_param <- function(attr, dose_i) {
    get_param(dat, paste0(prefix, attr, ".dose", dose_i))
  }
  
  vax.schedule <- data.frame(
    dose = seq_len(n_doses),
    interval = sapply(seq_len(n_doses), function(i) get_dose_param("interval", i)),
    annual = sapply(seq_len(n_doses), function(i) get_dose_param("annual", i)),
    
    rr.infect = sapply(seq_len(n_doses), function(i) get_dose_param("rr.infect", i)),
    rr.clinical = sapply(seq_len(n_doses), function(i) get_dose_param("rr.clinical", i)),
    rr.hosp = sapply(seq_len(n_doses), function(i) get_dose_param("rr.hosp", i)),
    
    ve.peak.infect = sapply(seq_len(n_doses), function(i) get_dose_param("ve.peak.infect", i)),
    ve.halflife.infect = sapply(seq_len(n_doses), function(i) get_dose_param("ve.halflife.infect", i)),
    ve.floor.infect = sapply(seq_len(n_doses), function(i) get_dose_param("ve.floor.infect", i)),
    
    ve.peak.clinical = sapply(seq_len(n_doses), function(i) get_dose_param("ve.peak.clinical", i)),
    ve.halflife.clinical = sapply(seq_len(n_doses), function(i) get_dose_param("ve.halflife.clinical", i)),
    ve.floor.clinical = sapply(seq_len(n_doses), function(i) get_dose_param("ve.floor.clinical", i)),
    
    ve.peak.hosp = sapply(seq_len(n_doses), function(i) get_dose_param("ve.peak.hosp", i)),
    ve.halflife.hosp = sapply(seq_len(n_doses), function(i) get_dose_param("ve.halflife.hosp", i)),
    ve.floor.hosp = sapply(seq_len(n_doses), function(i) get_dose_param("ve.floor.hosp", i)),

    # Direct severity/mortality protection (death-VE): reduces the disease-death
    # hazard for a vaccinated case who still develops severe disease. This is the
    # dominant, durable arm of the real vaccines (covid death-VE ~0.85, rsv ~0.80,
    # flu ~0.40) and is what lets age-targeting directly protect high-IFR groups.
    ve.peak.death = sapply(seq_len(n_doses), function(i) get_dose_param("ve.peak.death", i)),
    ve.halflife.death = sapply(seq_len(n_doses), function(i) get_dose_param("ve.halflife.death", i)),
    ve.floor.death = sapply(seq_len(n_doses), function(i) get_dose_param("ve.floor.death", i)),

    ve.delay = sapply(seq_len(n_doses), function(i) get_dose_param("ve.delay", i))
  )
  
  vax.schedule$start <- I(lapply(seq_len(n_doses), function(i) {
    get_dose_param("start", i)
  }))
  
  vax.schedule$rate <- I(lapply(seq_len(n_doses), function(i) {
    get_dose_param("rate", i)
  }))
  
  return(vax.schedule)
}

