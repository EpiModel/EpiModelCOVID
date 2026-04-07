
#' @rdname moduleset-common
#' @export
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
  n_firstdose_done <- sum(active == 1 & vax >= 1) # number of active people who already got at least 1 dose
  
  remaining_firstdose <- if (is.infinite(vax.supply.total)) { 
    Inf
  } else {  # how many more first doses are  allowed 
    floor(vax.supply.total * n_pop) - # total number of population receive vaccines
      n_firstdose_done
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
  
  # initialize these 4 strings for checking the number of people recieved boosters
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
      
      if (length(ids.vax1.boost) > 0) {
        vax[ids.vax1.boost] <- 1
        vax1Time[ids.vax1.boost] <- at
        remaining_supply <- remaining_supply - length(ids.vax1.boost)
        remaining_firstdose <- remaining_firstdose - length(ids.vax1.boost)
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
      
      if (length(ids.vax2.boost) > 0) {
        vax[ids.vax2.boost] <- 2
        vax2Time[ids.vax2.boost] <- at
        remaining_supply <- remaining_supply - length(ids.vax2.boost)
        
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
      
      if (length(ids.vax3.boost) > 0) {
        vax[ids.vax3.boost] <- 3
        vax3Time[ids.vax3.boost] <- at
        remaining_supply <- remaining_supply - length(ids.vax3.boost)
        
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
      
      if (length(ids.vax4.boost) > 0) {
        vax[ids.vax4.boost] <- 4
        vax4Time[ids.vax4.boost] <- at
        remaining_supply <- remaining_supply - length(ids.vax4.boost)
        
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
    # if (vax.strategy == "degree" && length(idsElig.vax1) > 0) browser() # check acceptance: When vax.strategy = "degree", highest-degree nodes are vaccinated first within each dose level
    # if (vax.strategy == "bridge" && length(idsElig.vax1) > 0) browser() # check acceptance: When vax.strategy = "bridge", cross-layer bridge nodes are vaccinated first
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
    
    if (nVax1) {
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
    
    if (nVax2) {
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
    if (length(idsVax4)) {
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
  
  # mnean deg and num of bridge in those ever vaccinated and never vaccinated
  ## individual ever vaccinated and never vaccinated
  vaccinated_ids <- which(active == 1 & vax >= 1)
  unvaccinated_ids <- which(active == 1 & vax == 0)
  ##  mean degree and bridge number in ever vaccinated individuals vs. in those unvaccinated
  mean_deg_vax <- mean(degree_total[vaccinated_ids])
  mean_deg_unvax <- mean(degree_total[unvaccinated_ids])
  mean_n_layers_vax <- mean(n_layers_active[vaccinated_ids])
  mean_n_layers_unvax <- mean(n_layers_active[unvaccinated_ids])
  
  dat <- set_epi(dat, "mean_deg_vax", at, mean_deg_vax)
  dat <- set_epi(dat, "mean_deg_unvax", at, mean_deg_unvax)
  dat <- set_epi(dat, "mean_n_layers_vax", at, mean_n_layers_vax)
  dat <- set_epi(dat, "mean_n_layers_unvax", at, mean_n_layers_unvax)

  # acceptance check, the mean degree of newly vaccinated individuals should be > that of eligible but not selected individuals.
  eligible_not_selected_ids <- setdiff(unique(ids_elig_all), unique(ids_newly_vaxed))
  
  mean_deg_newly_vaxed <- if (length(ids_newly_vaxed) > 0) {
    mean(degree_total[unique(ids_newly_vaxed)])
  } else NA_real_
  mean_deg_elig_not_selected <- if (length(eligible_not_selected_ids) > 0) {
    mean(degree_total[eligible_not_selected_ids])
  } else NA_real_
  mean_n_layers_newly_vaxed <- if (length(ids_newly_vaxed) > 0) {
      mean(n_layers_active[unique(ids_newly_vaxed)])
    } else NA_real_
  mean_n_layers_elig_not_selected <- if (length(eligible_not_selected_ids) > 0) {
      mean(n_layers_active[eligible_not_selected_ids])
    } else NA_real_
  
  dat <- set_epi(dat, "mean_deg_newly_vaxed", at, mean_deg_newly_vaxed)
  dat <- set_epi(dat, "mean_deg_elig_not_selected", at, mean_deg_elig_not_selected)
  dat <- set_epi(dat, "mean_n_layers_newly_vaxed", at, mean_n_layers_newly_vaxed)
  dat <- set_epi(dat, "mean_n_layers_elig_not_selected", at, mean_n_layers_elig_not_selected)

  
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
  
  if (vax.strategy == "degree") {
    if (remaining_supply <= 0) { # if no vaccine doses are left for this timestep-nobody is eligible and exit the function
      return(integer(0))
    }
    
    n_take <- min(remaining_supply, nElig) # per-timestep supply-cap mechanism - vaccineate no more than doses available & number of eligible people
    # degree-based allocatiom - rank eligible people by degree_total from high to low
    degree_jitter <- degree_total[idsElig] + runif(nElig, 0, 1e-8) # add small random noise breaks ties randomly
    
    names(degree_jitter) <- idsElig  # attach eligible ids as names so they stay linked to their degree_jitters
    
    idsElig_ranked <- as.integer(names(sort(degree_jitter, decreasing = TRUE))) # sort jittered degree from high to low, then recover the ranked IDs
    ids.vax <- idsElig_ranked[1:n_take] # take the top n_take people from the ranked list, i.e., nodes w/ highest degree are vaccinated first
    return(ids.vax)  
  } else if (vax.strategy == "bridge") {
    if (remaining_supply <= 0) {
      return(integer(0)) # same supply-cap mechanism as the degree strategy
    } 
    
    n_take <- min(remaining_supply, nElig) # same supply-cap mechanism as the degree strategy
    
    max_deg <- max(degree_total) + 1 # define a multiplier so n_layers_active has higher priority than degree_total
    
    score <- as.integer(is_bridge[idsElig]) * (5 * max_deg) + # bridge gets the biggest weight, 5 is (max(n_layers_active) +1)
      n_layers_active[idsElig] * max_deg + # then n_layers_active (i.e., 5> max(n_layers_active))
      degree_total[idsElig] +  # then degree_total (i.e., max_deg> degree_total)
      runif(nElig, 0, 1e-8) # tiny random noise to break ties
    
    names(score) <- idsElig  # attach eligible IDs as names
    
    idsElig_ranked <- as.integer(names(sort(score, decreasing = TRUE)))  # sort from high to low and get ranked IDs
    
    
    ids.vax <- idsElig_ranked[1:n_take] # same per-timestep supply cap as the degree strategy
    
    return(ids.vax)
    
  }   else if (vax.strategy == "age") {
    
    if (remaining_supply <= 0) {
      return(integer(0))
    } # same supply-cap mechanism as the degree strategy
    
    n_take <- min(remaining_supply, nElig)  # same supply-cap mechanism as the degree strategy
    
    age_priority <- rate[vax.age.group[idsElig]] # use the age-specific vaccination rate as the priority ordering
    
    age_priority_jitter <- age_priority + runif(nElig, 0, 1e-8) # higher rate = higher priority; add random noise to break ties randomly
    
    names(age_priority_jitter) <- idsElig # attach eligible IDs as names so they stay linked to their ordering
    
    idsElig_ranked <- as.integer(names(sort(age_priority_jitter, decreasing = TRUE)))# sort from high to low and get the ranked IDs
    
    ids.vax <- idsElig_ranked[1:n_take] # take the top n_take eligible people from the ranked 
    
    return(ids.vax)
    
  } else if (vax.strategy == "random") {
    
    if (remaining_supply <= 0) {
      return(integer(0))
    }# same supply-cap mechanism as the degree strategy
    
    n_take <- min(remaining_supply, nElig) # same supply-cap mechanism as the degree strategy
  
    
    ids.vax <- sample(idsElig, n_take) # choose n_take eligible people uniformly at random
    
    return(ids.vax) 
}
  else { stop("Unknown vax.strategy: ", vax.strategy) }
  
}



