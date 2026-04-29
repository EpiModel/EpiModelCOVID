#' @rdname moduleset-common
#' @export
vax_general <- function(dat, at) {
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  age <- get_attr(dat, "age")
  vax <- get_attr(dat, "vax")
  vax1Time <- get_attr(dat, "vax1Time")
  vax2Time <- get_attr(dat, "vax2Time")
  vax3Time <- get_attr(dat, "vax3Time")
  
  vax.age.breaks <- get_param(dat, "vax.age.breaks")
  vax.age.group <- cut(
    age,
    breaks = vax.age.breaks,
    right = FALSE,
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
  
  vax.schedule <- get_param(dat, "vax.schedule") # data.frame containing vax-related details, replacing individual get_param for each par
  
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
  
  # the following keeps the original logic in CorporateMix where booster campaign doses are attempted first.
  
  ######## 1. booster campaign doses initiated using boost.start ########
  
  # Track booster vaccinations separately from regular vaccinations
  daily_boost_by_dose <- rep(0L, nrow(vax.schedule)) # Track the number of people who receive each booster dose at this timestep.
  
  for (dose_i in seq_len(nrow(vax.schedule))) {
    
    # Pull the schedule row for this dose.
    # dose_row contains boost.start, boost.rate, start, interval...
    dose_row <- vax.schedule[dose_i, , drop = FALSE]
    
    # Pull age-specific booster campaign start times for this dose
    boost_start_vec <- dose_row$boost.start[[1]] # this is dose-specific because dose_row is the current row of vax.schedule
    
    # Identify people eligible for the booster campaign
    # This uses boost.start and dose_type = "boost" inside get_ids_eligible_for_dose().
    # vax1Time, vax2Time not required here
    if (any(at == boost_start_vec)) {
      
      idsEligBoost <- get_ids_eligible_for_dose(
        dose_row = dose_row,
        dose_i = dose_i,
        dose_type = "boost",
        active = active,
        status = status,
        dxStatus = dxStatus,
        dxTime = dxTime,
        vax = vax,
        vax.age.group = vax.age.group,
        at = at
      )
      
  
      nEligBoost <- length(idsEligBoost) # number of people eligible 
      
      if (nEligBoost > 0) { #  only run allocation if at least one person is eligible
      
      # Pull age-specific booster uptake probabilities for this dose
      boost_rate_vec <- dose_row$boost.rate[[1]]
      
      # For dose 1, both the daily supply cap and the cumulative first-dose cap apply.
      # For later doses, only the daily supply cap applies.
      effective_supply <- if (dose_i == 1) {
        min(remaining_supply, remaining_firstdose)
      } else {
        remaining_supply
      }
      
      # Select who receives the booster campaign dose using the allocation strategy
      idsBoost <- allocation_strategy(
        idsElig = idsEligBoost,
        rate = boost_rate_vec,
        vax.age.group = vax.age.group,
        vax.strategy = vax.strategy,
        degree_total = degree_total,
        is_bridge = is_bridge,
        n_layers_active = n_layers_active,
        remaining_supply = effective_supply
      )
      
      # Add eligible and vaccinated IDs to trackers
      ids_elig_all <- c(ids_elig_all, idsEligBoost)
      ids_newly_vaxed <- c(ids_newly_vaxed, idsBoost)
      
      # Count how many people received this booster dose.
      nBoost <- length(idsBoost)
      daily_boost_by_dose[dose_i] <- nBoost
      
      if (nBoost > 0) {
        
        # Update vaccine dose status.
        # Booster campaign pathway moves people to the same dose level as the regular pathway
        vax[idsBoost] <- dose_i
        
        # Record the time this dose was received.
        if (dose_i == 1) {
          vax1Time[idsBoost] <- at
        } else if (dose_i == 2) {
          vax2Time[idsBoost] <- at
        } else if (dose_i == 3) {
          vax3Time[idsBoost] <- at
        } else {
          stop("This version only supports dose_i = 1, 2, or 3.")
        }
        
        # Subtract used doses from the remaining daily supply.
        remaining_supply <- remaining_supply - nBoost
        
        # If these were first doses, also subtract from the first-dose cap.
        if (dose_i == 1) {
          remaining_firstdose <- remaining_firstdose - nBoost
        }
      }
      }
    }
  }
  ######## 2. regular doses initiated using start, interval ########
  daily_vax_by_dose <- rep(0L, nrow(vax.schedule)) # length is the number of doses
  
  for (dose_i in seq_len(nrow(vax.schedule)) # for each (dose_i) dose
       ) {
    
    dose_row <- vax.schedule[dose_i, , drop = FALSE] 
    
    idsElig <- get_ids_eligible_for_dose(
      dose_row = dose_row,
      dose_i = dose_i,
      dose_type = "regular",
      active = active,
      status = status,
      dxStatus = dxStatus,
      dxTime = dxTime,
      vax = vax,
      vax1Time = vax1Time,
      vax2Time = vax2Time,
      vax.age.group = vax.age.group,
      at = at
    )
    
    nElig <- length(idsElig)
    
    if (nElig > 0) {
      
    rate_vec <- dose_row$rate[[1]]
    
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
      remaining_supply = effective_supply
    )
    
    ids_elig_all <- c(ids_elig_all, idsElig)
    ids_newly_vaxed <- c(ids_newly_vaxed, idsVax)
    
    nVax <- length(idsVax)
    daily_vax_by_dose[dose_i] <- nVax
    
    if (nVax > 0) {
      vax[idsVax] <- dose_i
      if (dose_i == 1) {
        vax1Time[idsVax] <- at
      } else if (dose_i == 2) {
        vax2Time[idsVax] <- at
      } else if (dose_i == 3) {
        vax3Time[idsVax] <- at
      } else {
        stop("This version only supports dose_i = 1, 2, or 3.")
      }
      remaining_supply <- remaining_supply - nVax
      
      if (dose_i == 1) {
        remaining_firstdose <- remaining_firstdose - nVax
      }
    }
    }
  }
  
  # Replace attr
  dat <- set_attr(dat, "vax", vax)
  dat <- set_attr(dat, "vax1Time", vax1Time)
  dat <- set_attr(dat, "vax2Time", vax2Time)
  dat <- set_attr(dat, "vax3Time", vax3Time)
  
  # Summary statistics 
  ## Number of regular doses received at this timestep
  for (dose_i in seq_len(nrow(vax.schedule))) {
    dat <- set_epi(dat, paste0("nVax", dose_i), at, daily_vax_by_dose[dose_i])
  }
  
  ## Number of booster/campaign doses received at this timestep
  ### This is not in the old script
  for (dose_i in seq_len(nrow(vax.schedule))) {
    dat <- set_epi(dat, paste0("nBoost", dose_i), at, daily_boost_by_dose[dose_i])
  }
  
  ## Age-specific vax coverage
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
  daily_vax_total <- sum(daily_vax_by_dose) + sum(daily_boost_by_dose)
  # daily_vax_total <- nVax1 + nVax2 + nVax3 + nVax4 + # old syntax
  #   length(ids.vax1.boost) + length(ids.vax2.boost) +
  #   length(ids.vax3.boost) + length(ids.vax4.boost)
  
  
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
                                remaining_supply) {
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
    
  } else { 
    stop("Unknown vax.strategy: ", vax.strategy) 
  }
  
  # Second step: determine vaccination allocation
  ## Get each ranked eligible person's age-specific vax probability.
  rate_person <- rate[vax.age.group[ids_ranked]]

  if (vax.strategy %in% c("degree", "bridge")) {
    ## For network-based priority strategies: stochastic demand (driven by
    ## age-specific rates) determines how many get vaccinated each timestep;
    ## the priority ranking determines who.  This ensures high-priority
    ## individuals are selected even when demand < supply.
    n_willing <- sum(rbinom(nElig, size = 1, prob = rate_person))
    if (n_willing == 0) return(integer(0))
    n_take <- min(remaining_supply, n_willing)
    ids.vax <- ids_ranked[1:n_take]
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

# function to calculation each person's current (at) vaccine-related relative risk
update_vax_waning_attrs <- function(dat, at, vax.schedule, vax, vax_times,
                                    vax.age.group, active) {
  
  # dat = EpiModel simulation object inside netsim
  # at = current simulation time step
  # vax.schedule = data frame defining dose schedule and VE parameters
  # vax = each person's current highest dose received
  # vax_times = list of vaccination time vectors, e.g., vax1Time, vax2Time
  # vax.age.group = each person's vaccine age group
  # active = whether each person is active in the simulation
  
  # Total number of people in the simulation
  n <- length(vax)
  
  # Number of vaccine age groups
  n_age <- max(vax.age.group, na.rm = TRUE)
  
  # Initialize everyone as having no vaccine protection
  # rr_infect = 1 means no reduction in infection risk
  rr_infect <- rep(1, n)
  
  # Loop over each dose row in vax.schedule
  for (i in seq_len(nrow(vax.schedule))) {
    
    # Get the dose number for this schedule row
    dose <- vax.schedule$dose[i]
    
    # Identify active people who received at least this dose
    ids <- which(active == 1 & vax >= dose & vax_times[[dose]] > -Inf)
    
    # If nobody received this dose, skip to the next schedule row
    if (length(ids) == 0) {
      next
    }
    
    # Extract age-specific peak VE for this dose
    ve_peak <- sched_vec("ve.infect.peak", i, default = 0)
    
    # Extract age-specific VE half-life for this dose
    ve_halflife <- sched_vec("ve.infect.halflife", i, default = Inf)
    
    # Extract age-specific VE floor for this dose
    ve_floor <- sched_vec("ve.infect.floor", i, default = 0)
    
    # Calculate days since each selected person received this dose
    t_since <- at - vax_times[[dose]][ids]
    
    # Calculate current VE for each selected person
    ve_now <- ve_decay(
      t_since = t_since,
      ve.peak = ve_peak[vax.age.group[ids]],
      ve.halflife = ve_halflife[vax.age.group[ids]],
      ve.floor = ve_floor[vax.age.group[ids]]
    )
    
    # Convert VE to relative risk:
    # VE = 0.60 means relative risk = 0.40
    # If multiple dose rows apply, keep the strongest protection
    rr_infect[ids] <- pmin(rr_infect[ids], 1 - ve_now)
  }
  
  # Store the current vaccine-related infection relative risk as a node attribute
  dat <- set_attr(dat, "rr_infect_vax", rr_infect)
  
  # Return the updated dat object
  return(dat)
}


get_ids_eligible_for_dose <- function(
    dose_row, dose_i, 
    dose_type, # "regular" or "boost"
    active, status, dxStatus, dxTime,
    # For 3-dose COVID, only need vax1Time and vax2Time for eligibility checking, 
    # because dose 1 has no previous dose, dose 2 checks vax1Time, and dose 3 checks vax2Time. 
    vax, 
    vax1Time = NULL, 
    vax2Time = NULL, 
    vax.age.group, at
) {
  # Define people who are generally eligible for vaccination, regardless of dose number.
  # active == 1: person is active in the simulation.
  # status not in c("ic", "h"): exclude severe individuals.
  # ! dxStatus == 2 and recent dxTime: exclude recently diagnosed people.
  base_eligible <- active == 1 &
    !(status %in% c("ic", "h")) &
    !(dxStatus == 2 & (at - dxTime <= 10))
  
  if (dose_type == "boost") {
    
    # Pull the age-specific booster campaign start times for this dose.
    boost_start_vec <- dose_row$boost.start[[1]]
    
    # Identify people eligible for the booster/campaign pathway.
    # base_eligible: active, not hospitalized/ICU, and not recently diagnosed.
    # vax == dose_i - 1: person must have completed the previous dose level.
    #   For dose_i = 1, this means vax == 0.
    #   For dose_i = 2, this means vax == 1.
    #   For dose_i = 3, this means vax == 2.
    # at == boost_start_vec[vax.age.group]: booster is only offered on the
    # exact age-specific campaign start day, matching the original script.
    idsElig <- which(
      base_eligible & # person must be generally eligible
        vax == dose_i - 1 & # person must have completed the previous dose level (regular or boost)
        at == boost_start_vec[vax.age.group] #  booster is only offered on the exact age-specific campaign start day
    )
    
    return(idsElig)
  }
  
  if (dose_type == "regular") {
  # Pull the age-specific start times for this dose from vax.schedule
  # start_vec has length of 5
  # Example: start_vec[1] is the start time for age group 0–4.
  start_vec <- dose_row$start[[1]]
  
  # Pull the minimum required interval since the previous dose.
  # NA for dose 1 because there is no previous dose.
  interval <- dose_row$interval

  # Dose 1 eligibility
  if (dose_i == 1) {
    
    idsElig <- which(
      base_eligible & # person must be generally eligible
        vax == 0 & # have received no prior vaccine dose
        at >= start_vec[vax.age.group]  # and have reached their age-specific dose 1 start time
    )
    
  } else { # Dose 2 or 3 eligibility determination
    
    # For dose 2 or dose 3, identify the previous dose time.
    # If dose_i == 2, previous dose is dose 1.
    # If dose_i == 3, previous dose is dose 2.
    prev_time <- if (dose_i == 2) {
      vax1Time
    } else if (dose_i == 3) {
      vax2Time
    } else {
      stop("This version only supports dose_i = 1, 2, or 3.")
    }
    
    # Dose 2 or 3 eligibility
    idsElig <- which(
      base_eligible &  # person must be generally eligible
        vax == dose_i - 1 & # have received exactly the previous dose (regular or boost)
        !is.na(prev_time) & # have a recorded previous dose time
        (at - prev_time >= interval) & # have waited long enough since the previous dose
        at >= start_vec[vax.age.group]  # and have reached their age-specific start time for this dose
    )
  }
  
  # Return IDs of eligible individuals for this dose at this time step.
  return(idsElig)
  }
}


# Helper function to extract a schedule column as an age-specific vector
sched_vec <- function(col, i, default = NA_real_) {
  
  # If the column does not exist, return the default value for all age groups
  if (!col %in% names(vax.schedule)) {
    return(rep(default, n_age))
  }
  
  # If this column is a list-column, extract the vector from row i
  x <- if (is.list(vax.schedule[[col]])) {
    vax.schedule[[col]][[i]]
    
    # Otherwise, extract the single value from row i
  } else {
    vax.schedule[[col]][i]
  }
  
  # If only one value is provided, recycle it to all age groups
  if (length(x) == 1) {
    rep(x, n_age)
    
    # If a vector is provided, use it directly
  } else {
    x
  }
}

# Calculate vaccine effectiveness at a given time since vaccination
ve_decay <- function(t_since, ve.peak, ve.halflife, ve.floor, delay = 14) {
  
  # t_since = number of days since the person received the dose
  # ve.peak = maximum vaccine effectiveness after the dose
  # ve.halflife = number of days for VE to decline by 50%
  # ve.floor = minimum remaining VE after waning
  # delay = days before full vaccine effect begins, default is 14 days
  
  # Before the delay period ends, treat effective waning time as 0
  t_eff <- pmax(0, t_since - delay)
  
  # Exponential waning:
  # VE starts near ve.peak and gradually declines toward ve.floor
  ve.floor + (ve.peak - ve.floor) * 0.5^(t_eff / ve.halflife)
}

