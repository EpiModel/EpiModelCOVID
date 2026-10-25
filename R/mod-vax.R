#' @rdname moduleset-gmc19
#' @export
vax_general <- function(dat, at) {
  ## The allocation model (GlobalMix India #74): one campaign for every arm.
  ## Doses come from a cumulative budget, vax.budget x the living population on
  ## the day the campaign opens (vax.start), counted in doses administered and
  ## never refilled when a recipient dies, which an exogenous throughput
  ## schedule releases over vax.delivery.days calendar days (vax_allowed()):
  ## every arm in a cell delivers the same doses on the same days, and a
  ## strategy can only reorder who gets them. Every arm allocates only among
  ## the willing eligibles: alive, susceptible or recovered (the exposed, the
  ## asymptomatic and the presymptomatic are not dosed; they used to be, and
  ## got severity protection for the episode under way), inside the product's
  ## age window (vax.<pathogen>.age.min and .age.max, licensure by #58 item
  ## 10), and willing (vax.willing, a flag drawn once per person by vaccine
  ## age group at initialization or birth). A strategy is an ordering of that
  ## set and nothing else (allocation_strategy()). The supply caps
  ## (vax.supply.rate, vax.supply.total), the per-age uptake rates and start
  ## days (vax.<pathogen>.rate.dose<i>, .start.dose<i>) and vax.infant.eligible
  ## are gone: the caps counted the living vaccinated, so a recipient's death
  ## freed a dose and the age arm delivered 35% more doses than the degree arm
  ## (#54); the uptake rates made the age arm rank by rate (#55) and gave the
  ## network arms perfect compliance while the age and random arms waited for
  ## show-ups (#56). The strategy set is #75's: what an arm ranks on is fixed
  ## when the campaign opens where the arm is frozen (vax_campaign_rank()). The
  ## infant products are a dose line of their own beside this budget (#76,
  ## infant_program_catchup() here and infant_program_births() at birth), and
  ## the cocoon is an ordering of this budget like any other.

  ######## attributes ########
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  vax <- get_attr(dat, "vax")
  last.dose.time <- get_attr(dat, "last.dose.time")
  age <- get_attr(dat, "age")
  willing <- get_attr(dat, "vax.willing")
  ever.inf <- get_attr(dat, "ever.inf")
  vax.age.group <- vax_age_group_for(dat)

  degree_total <- get_attr(dat, "degree_total")
  # The score the degree strategies rank on (vax.degree.score; netdegree() sets
  # it, and it is degree_total unless the parameter names another score)
  degree_rank <- get_attr(dat, "degree_score", override.null.error = TRUE)
  if (is.null(degree_rank)) degree_rank <- degree_total
  is_bridge <- get_attr(dat, "is_bridge")
  n_layers_active <- get_attr(dat, "n_layers_active")
  is_infant <- get_attr(dat, "is_infant")
  uid <- get_attr(dat, "unique_id")
  # the cocoon tier (#76): the adults of the households with a young infant today
  tier <- cocoon_tier(dat)

  ######## parameters ########
  vax.strategy <- get_param(dat, "vax.strategy")
  season.length <- param_or(dat, "vax.season.length", days_per_year)
  pa <- get_pathogen(dat)

  # The dose schedule: one row per dose with interval, annual and the VE
  # profile (read downstream through compute_ve())
  vax.schedule <- build_vax_schedule(dat, pa)

  ######## the campaign ########
  dat <- vax_campaign_open(dat, at, active)
  # the infant program's catch-up, on the day the campaign opens (#76)
  cu <- infant_program_catchup(dat, at)
  dat <- cu$dat
  camp <- dat$run$vax.campaign
  start <- get_param(dat, "vax.start")
  n_today <- if (is.null(camp)) 0L else max(0L, vax_allowed(camp, at) - camp$given)

  # Eligibility shared by every strategy, before any ranking
  elig_any <- vax_eligible(dat, pa, active, status, age)
  eligible <- elig_any & willing == 1
  if (anyNA(eligible)) {
    stop("vaccine eligibility is NA for ", sum(is.na(eligible)), " node(s) (GlobalMix India #74)", call. = FALSE)
  }

  # What the strategy ranks on (#75), fixed on the day the campaign opens where
  # the arm is frozen, and the hybrid's reserve doses for today
  if (!is.null(camp) && is.null(camp$ranked)) {
    camp <- vax_campaign_rank(dat, camp, vax.strategy, eligible, degree_rank, age)
  }
  score <- if (n_today > 0) {
    if (vax.strategy == "cocoon") tier else strategy_score(dat, camp, vax.strategy, degree_rank, age, uid)
  }
  n_cocoon_elig <- sum(tier == 1 & eligible & vax == 0)
  n_reserve <- if (vax.strategy == "hybrid" && !is.null(camp)) hybrid_quota(camp, at) else 0L
  reserve.age <- param_or(dat, "vax.hybrid.age", 60)

  ids_newly_vaxed <- integer(0)
  ids_elig_all <- integer(0)
  n_elig_any <- 0L
  nVax_reserve <- 0L
  daily_vax_by_dose <- rep(0L, nrow(vax.schedule))

  for (dose_i in seq_len(nrow(vax.schedule))) {
    dose_row <- vax.schedule[dose_i, , drop = FALSE]
    n_elig_any <- n_elig_any + length(get_ids_eligible_for_dose(
      dose_row, dose_i, elig_any, vax, last.dose.time, at, start, season.length))
    idsElig <- get_ids_eligible_for_dose(dose_row, dose_i, eligible, vax, last.dose.time, at, start, season.length)
    ids_elig_all <- c(ids_elig_all, idsElig)
    if (length(idsElig) == 0 || n_today == 0) next
    idsVax <- allocation_strategy(
      idsElig = idsElig, n_take = n_today, vax.strategy = vax.strategy, age = age, score = score,
      n_reserve = n_reserve, reserve.age = reserve.age)
    n_res <- if (is.null(attr(idsVax, "n_reserve"))) 0L else attr(idsVax, "n_reserve")
    idsVax <- as.vector(idsVax)
    n_reserve <- n_reserve - n_res
    nVax_reserve <- nVax_reserve + n_res
    nVax <- length(idsVax)
    if (nVax > n_today || any(!idsVax %in% idsElig) || anyDuplicated(idsVax)) {
      stop("vax.strategy ", vax.strategy, " handed out ", nVax, " doses for ", n_today,
           " allowed, or outside the eligible set (GlobalMix India #74)", call. = FALSE)
    }
    daily_vax_by_dose[dose_i] <- nVax
    if (nVax > 0) {
      vax[idsVax] <- dose_i
      last.dose.time[idsVax] <- at
      n_today <- n_today - nVax
      camp$given <- camp$given + nVax
      ids_newly_vaxed <- c(ids_newly_vaxed, idsVax)
    }
  }
  if (!is.null(camp$reserve)) camp$reserve_given <- camp$reserve_given + nVax_reserve
  if (!is.null(camp)) dat$run$vax.campaign <- camp

  dat <- set_attr(dat, "vax", vax)
  dat <- set_attr(dat, "last.dose.time", last.dose.time)

  ######## summary statistics ########
  ## doses by dose number, and the dose accounting (#74; #77 reads it)
  for (dose_i in seq_len(nrow(vax.schedule))) {
    dat <- set_epi(dat, paste0("nVax", dose_i), daily_vax_by_dose[dose_i], at = at)
  }
  daily_vax_total <- sum(daily_vax_by_dose)
  dat <- set_epi(dat, "daily_vax_total", daily_vax_total, at = at)
  dat <- set_epi(dat, "vax.doses", if (is.null(camp)) 0L else camp$given, at = at)
  dat <- set_epi(dat, "vax.budget.doses", if (is.null(camp)) NA_integer_ else camp$budget, at = at)
  dat <- set_epi(dat, "vax.allowed", if (is.null(camp)) 0L else vax_allowed(camp, at), at = at)
  dat <- set_epi(dat, "n.elig", n_elig_any, at = at)
  dat <- set_epi(dat, "n.elig.willing", length(ids_elig_all), at = at)
  ## recipients by their state at the dose: susceptible and never infected,
  ## susceptible after an infection, recovered
  st_new <- status[ids_newly_vaxed]
  dat <- set_epi(dat, "nVax.sus", sum(st_new == "s" & ever.inf[ids_newly_vaxed] == 0), at = at)
  dat <- set_epi(dat, "nVax.prior", sum(st_new == "s" & ever.inf[ids_newly_vaxed] == 1), at = at)
  dat <- set_epi(dat, "nVax.rec", sum(st_new == "r"), at = at)
  ## and two the eligibility rule should keep empty or small (#77): recipients
  ## in an infected state at the dose (none under #74, which doses s and r only)
  ## and recipients immune from before the run (#70; in r for the whole run, so
  ## a subset of nVax.rec). The four states plus nVax.inf sum to daily_vax_total.
  dat <- set_epi(dat, "nVax.inf", sum(st_new %in% c("e", "a", "ip", "ic", "h")), at = at)
  prior.immune <- get_attr(dat, "prior.immune", override.null.error = TRUE)
  dat <- set_epi(dat, "nVax.immune",
                 if (is.null(prior.immune)) 0L else sum(prior.immune[ids_newly_vaxed] %in% 1), at = at)
  ## Doses by model age band, infants and vaccine age group, per step (#77).
  ## Summed over the run they are cumulative doses by band whatever happens to
  ## the recipients afterwards; the cov_vax* shares below count the living
  ## vaccinated over the living, so they fall when a vaccinated person dies.
  age.breaks <- get_param(dat, "age.breaks")
  ag_new <- cut(age[ids_newly_vaxed], age.breaks, labels = FALSE, right = FALSE)
  for (g in seq_len(length(age.breaks) - 1)) {
    dat <- set_epi(dat, paste0("nVax.age", g), sum(ag_new == g, na.rm = TRUE), at = at)
  }
  dat <- set_epi(dat, "nVax.infant",
                 if (is.null(is_infant)) 0L else sum(is_infant[ids_newly_vaxed] %in% TRUE), at = at)
  vg_new <- vax.age.group[ids_newly_vaxed]
  for (k in seq_len(length(get_param(dat, "vax.age.breaks")) - 1)) {
    dat <- set_epi(dat, paste0("nVax.vaxage", k), sum(vg_new == k, na.rm = TRUE), at = at)
  }
  ## The day the campaign opens (the first step at or after vax.start), with the
  ## eligible and the willing eligible pools that day, which size the budget
  if (!is.null(camp) && at == max(camp$start, 2)) {
    dat <- set_epi(dat, "alloc.day", at, at = at)
    dat <- set_epi(dat, "n.elig.alloc", n_elig_any, at = at)
    dat <- set_epi(dat, "n.elig.willing.alloc", length(ids_elig_all), at = at)
    ## the precision arm's realized rank correlation with the score over that
    ## day's willing eligibles, and the hybrid's reserve (#75)
    dat <- set_epi(dat, "vax.rank.spearman.alloc", if (is.null(camp$spearman)) NA_real_ else camp$spearman, at = at)
    dat <- set_epi(dat, "vax.reserve.doses", if (is.null(camp$reserve)) NA_integer_ else camp$reserve, at = at)
  }
  ## doses from the hybrid's 60+ reserve (#75; 0 on every other arm)
  dat <- set_epi(dat, "nVax.reserve", nVax_reserve, at = at)
  ## The cocoon tier (#76), on every arm: its living members, those of them
  ## willing, eligible and undosed before today's doses, and the doses they got
  dat <- set_epi(dat, "n.cocoon.tier", sum(tier == 1), at = at)
  dat <- set_epi(dat, "n.cocoon.elig", n_cocoon_elig, at = at)
  dat <- set_epi(dat, "nVax.cocoon", sum(tier[ids_newly_vaxed] == 1), at = at)
  ## The infant program's dose line (#76)
  dat <- infant_program_outputs(dat, at, cu$n)

  ## Age-specific vax coverage #TODO: this needs to be turned to loop like the above
  dat <- set_epi(dat, "cov_vax1_0to4", length(which(vax.age.group == 1 & vax >= 1)) / length(which(vax.age.group == 1)), at = at)
  dat <- set_epi(dat, "cov_vax1_5to17", length(which(vax.age.group == 2 & vax >= 1)) / length(which(vax.age.group == 2)), at = at)
  dat <- set_epi(dat, "cov_vax1_18to64", length(which((vax.age.group == 3 | vax.age.group == 4) & vax >= 1)) / length(which(vax.age.group == 3 | vax.age.group == 4)), at = at)
  dat <- set_epi(dat, "cov_vax1_65p", length(which(vax.age.group == 5 & vax >= 1)) / length(which(vax.age.group == 5)), at = at)

  dat <- set_epi(dat, "cov_vax2_0to4", length(which(vax.age.group == 1 & vax >= 2)) / length(which(vax.age.group == 1)), at = at)
  dat <- set_epi(dat, "cov_vax2_5to17", length(which(vax.age.group == 2 & vax >= 2)) / length(which(vax.age.group == 2)), at = at)
  dat <- set_epi(dat, "cov_vax2_18to64", length(which((vax.age.group == 3 | vax.age.group == 4) & vax >= 2)) / length(which(vax.age.group == 3 | vax.age.group == 4)), at = at)
  dat <- set_epi(dat, "cov_vax2_65p", length(which(vax.age.group == 5 & vax >= 2)) / length(which(vax.age.group == 5)), at = at)

  dat <- set_epi(dat, "cov_vax3_0to4", length(which(vax.age.group == 1 & vax >= 3)) / length(which(vax.age.group == 1)), at = at)
  dat <- set_epi(dat, "cov_vax3_5to17", length(which(vax.age.group == 2 & vax >= 3)) / length(which(vax.age.group == 2)), at = at)
  dat <- set_epi(dat, "cov_vax3_18to49", length(which(vax.age.group == 3 & vax >= 3)) / length(which(vax.age.group == 3)), at = at)
  dat <- set_epi(dat, "cov_vax3_50to64", length(which(vax.age.group == 4 & vax >= 3)) / length(which(vax.age.group == 4)), at = at)
  dat <- set_epi(dat, "cov_vax3_65p", length(which(vax.age.group == 5 & vax >= 3)) / length(which(vax.age.group == 5)), at = at)

  ## Dose-1 coverage by NETWORK POSITION (fixed degree_quartile / is_bridge attrs).
  ## Verifies the network strategies actually reach high-degree / bridge nodes
  ## (the degree-strategy acceptance criterion, Q4 > Q1; closes issue #21).
  degree_quartile <- get_attr(dat, "degree_quartile")
  if (!is.null(degree_quartile)) {
    for (q in 1:4) dat <- set_epi(dat, paste0("cov_vax1_degQ", q),
      length(which(degree_quartile == q & vax >= 1)) / max(1, length(which(degree_quartile == q))), at = at)
  }
  # and by the baseline quartile fixed on day 2 (#77; netdegree()): the daily
  # quartile follows the current score, so a hub's coverage read on it moves
  # with the hub's ties
  degree_quartile_base <- get_attr(dat, "degree_quartile_base", override.null.error = TRUE)
  if (!is.null(degree_quartile_base)) {
    for (q in 1:4) dat <- set_epi(dat, paste0("cov_vax1_degQb", q),
      length(which(degree_quartile_base == q & vax >= 1)) / max(1, length(which(degree_quartile_base == q))), at = at)
  }
  if (!is.null(is_bridge)) {
    dat <- set_epi(dat, "cov_vax1_bridge",
      length(which(is_bridge == 1 & vax >= 1)) / max(1, length(which(is_bridge == 1))), at = at)
    dat <- set_epi(dat, "cov_vax1_nonbridge",
      length(which(is_bridge == 0 & vax >= 1)) / max(1, length(which(is_bridge == 0))), at = at)
  }

  # mean deg and num of bridge in those ever vaccinated and never vaccinated
  vaccinated_ids <- which(active == 1 & vax >= 1)
  unvaccinated_ids <- which(active == 1 & vax == 0)
  mean_deg_vax <- mean(degree_rank[vaccinated_ids])
  mean_deg_unvax <- mean(degree_rank[unvaccinated_ids])
  mean_n_layers_active_vax <- mean(n_layers_active[vaccinated_ids])
  mean_n_layers_active_unvax <- mean(n_layers_active[unvaccinated_ids])

  ## cumulative coverage rate
  cov_cumulative <- length(vaccinated_ids) / sum(active == 1)

  dat <- set_epi(dat, "mean_deg_vax", mean_deg_vax, at = at)
  dat <- set_epi(dat, "mean_deg_unvax", mean_deg_unvax, at = at)
  dat <- set_epi(dat, "mean_n_layers_active_vax", mean_n_layers_active_vax, at = at)
  dat <- set_epi(dat, "mean_n_layers_active_unvax", mean_n_layers_active_unvax, at = at)
  dat <- set_epi(dat, "cov_cumulative", cov_cumulative, at = at)

  # acceptance check, the mean degree of newly vaccinated individuals should be > that of eligible but not selected individuals.
  eligible_not_selected_ids <- setdiff(unique(ids_elig_all), unique(ids_newly_vaxed))

  # this tracker pools across all dose levels per timestep. If we ever need per-dose-level verification, we'd need separate trackers.
  mean_deg_newly_vaxed <- if (length(ids_newly_vaxed) > 0) {
    mean(degree_rank[unique(ids_newly_vaxed)])
  } else NA_real_
  mean_deg_elig_not_selected <- if (length(eligible_not_selected_ids) > 0) {
    mean(degree_rank[eligible_not_selected_ids])
  } else NA_real_
  mean_n_layers_active_newly_vaxed <- if (length(ids_newly_vaxed) > 0) {
    mean(n_layers_active[unique(ids_newly_vaxed)])
  } else NA_real_
  mean_n_layers_active_elig_not_selected <- if (length(eligible_not_selected_ids) > 0) {
    mean(n_layers_active[eligible_not_selected_ids])
  } else NA_real_

  dat <- set_epi(dat, "mean_deg_newly_vaxed", mean_deg_newly_vaxed, at = at)
  dat <- set_epi(dat, "mean_deg_elig_not_selected", mean_deg_elig_not_selected, at = at)
  dat <- set_epi(dat, "mean_n_layers_active_newly_vaxed", mean_n_layers_active_newly_vaxed, at = at)
  dat <- set_epi(dat, "mean_n_layers_active_elig_not_selected", mean_n_layers_active_elig_not_selected, at = at)

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
  dat <- set_epi(dat, "nVax1", nVax1, at = at)
  dat <- set_epi(dat, "nVax2", nVax2, at = at)
  dat <- set_epi(dat, "nVax3", nVax3, at = at)
  dat <- set_epi(dat, "nVax4", nVax4, at = at)
  
  ## Summary stats -- vaccine coverage ##
  dat <- set_epi(dat, "cov_vax1_0to4", length(which(vax.age.group == 1 & vax >= 1)) / length(which(vax.age.group == 1)), at = at)
  dat <- set_epi(dat, "cov_vax1_5to17", length(which(vax.age.group == 2 & vax >= 1)) / length(which(vax.age.group == 2)), at = at)
  dat <- set_epi(dat, "cov_vax1_18to64", length(which((vax.age.group == 3 | vax.age.group == 4) & vax >= 1)) / length(which(vax.age.group == 3 | vax.age.group == 4)), at = at)
  dat <- set_epi(dat, "cov_vax1_65p", length(which(vax.age.group == 5 & vax >= 1)) / length(which(vax.age.group == 5)), at = at)

  dat <- set_epi(dat, "cov_vax2_0to4", length(which(vax.age.group == 1 & vax >= 2)) / length(which(vax.age.group == 1)), at = at)
  dat <- set_epi(dat, "cov_vax2_5to17", length(which(vax.age.group == 2 & vax >= 2)) / length(which(vax.age.group == 2)), at = at)
  dat <- set_epi(dat, "cov_vax2_18to64", length(which((vax.age.group == 3 | vax.age.group == 4) & vax >= 2)) / length(which(vax.age.group == 3 | vax.age.group == 4)), at = at)
  dat <- set_epi(dat, "cov_vax2_65p", length(which(vax.age.group == 5 & vax >= 2)) / length(which(vax.age.group == 5)), at = at)

  dat <- set_epi(dat, "cov_vax3_5to17", length(which(vax.age.group == 2 & vax >= 3)) / length(which(vax.age.group == 2)), at = at)
  dat <- set_epi(dat, "cov_vax3_18to49", length(which(vax.age.group == 3 & vax >= 3)) / length(which(vax.age.group == 3)), at = at)
  dat <- set_epi(dat, "cov_vax3_50to64", length(which(vax.age.group == 4 & vax >= 3)) / length(which(vax.age.group == 4)), at = at)
  dat <- set_epi(dat, "cov_vax3_65p", length(which(vax.age.group == 5 & vax >= 3)) / length(which(vax.age.group == 5)), at = at)

  dat <- set_epi(dat, "cov_vax4_50to64", length(which(vax.age.group == 4 & vax >= 4)) / length(which(vax.age.group == 4)), at = at)
  dat <- set_epi(dat, "cov_vax4_65p", length(which(vax.age.group == 5 & vax >= 4)) / length(which(vax.age.group == 5)), at = at)
  
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
  
  dat <- set_epi(dat, "mean_deg_vax", mean_deg_vax, at = at)
  dat <- set_epi(dat, "mean_deg_unvax", mean_deg_unvax, at = at)
  dat <- set_epi(dat, "mean_n_layers_active_vax", mean_n_layers_active_vax, at = at)
  dat <- set_epi(dat, "mean_n_layers_active_unvax", mean_n_layers_active_unvax, at = at)
  dat <- set_epi(dat, "cov_cumulative", cov_cumulative, at = at)

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
  
  dat <- set_epi(dat, "mean_deg_newly_vaxed", mean_deg_newly_vaxed, at = at)
  dat <- set_epi(dat, "mean_deg_elig_not_selected", mean_deg_elig_not_selected, at = at)
  dat <- set_epi(dat, "mean_n_layers_active_newly_vaxed", mean_n_layers_active_newly_vaxed, at = at)
  dat <- set_epi(dat, "mean_n_layers_active_elig_not_selected", mean_n_layers_active_elig_not_selected, at = at)

  
  # acceptance check, Verified: with a cap of 0.005, daily vaccination count never exceeds 0.5% of population size
  daily_vax_total <- nVax1 + nVax2 + nVax3 + nVax4 +
    length(ids.vax1.boost) + length(ids.vax2.boost) +
    length(ids.vax3.boost) + length(ids.vax4.boost)
  
  
  daily_vax_cap <- if (is.infinite(vax.supply.rate)) {
    Inf
  } else {
    floor(vax.supply.rate * n_pop)
  }
  
  dat <- set_epi(dat, "remaining_supply", remaining_supply, at = at)
  dat <- set_epi(dat, "remaining_firstdose", remaining_firstdose, at = at)
  dat <- set_epi(dat, "daily_vax_total", daily_vax_total, at = at)
  dat <- set_epi(dat, "daily_vax_cap", daily_vax_cap, at = at)
  
  return(dat)
}

# The campaign (GlobalMix India #74). It opens on vax.start with a budget of
# round(vax.budget x the living) doses, fixed then and kept in dat$run, with the
# doses administered so far; a recipient's death returns nothing. The
# throughput schedule releases the budget evenly over vax.delivery.days
# calendar days: by day t of the campaign, floor(budget x min(1, t / days))
# doses may have been given, so a day whose willing eligibles fall short of its
# release carries the shortfall forward.
vax_campaign_open <- function(dat, at, active) {
  if (!is.null(dat$run$vax.campaign)) return(dat)
  start <- get_param(dat, "vax.start")
  if (at < start) return(dat)
  dat$run$vax.campaign <- list(
    start = start, days = get_param(dat, "vax.delivery.days"),
    budget = as.integer(round(get_param(dat, "vax.budget") * sum(active == 1))), given = 0L)
  dat
}

vax_allowed <- function(camp, at) {
  as.integer(floor(camp$budget * min(1, (at - camp$start + 1) / camp$days)))
}

# Who a product may be given to, before willingness and the dose rules: alive,
# susceptible or recovered, and inside the product's age window (age.min
# inclusive, age.max exclusive, in years; licensure by #58 item 10). The
# exposed, the asymptomatic and the presymptomatic are not dosed, nor are the
# symptomatic and the admitted. NA stops the run (a waned, never-tested
# susceptible used to turn NA and drop out for good).
vax_eligible <- function(dat, pa, active, status, age) {
  amin <- get_param(dat, paste0("vax.", pa, ".age.min"))
  amax <- get_param(dat, paste0("vax.", pa, ".age.max"))
  e <- active == 1 & status %in% c("s", "r") & age >= amin & age < amax
  if (anyNA(e)) {
    stop("vaccine eligibility is NA for ", sum(is.na(e)), " node(s) (GlobalMix India #74)", call. = FALSE)
  }
  e
}

# The eligibles for one dose: the shared eligibility (vax_eligible(), and
# willingness) plus the dose rules. Dose 1 goes to the undosed once the campaign
# has started; a later dose needs exactly the previous dose and the interval
# since it; an annual dose reopens every season.length days from the campaign
# start to whoever has not been dosed since the season opened.
get_ids_eligible_for_dose <- function(dose_row, dose_i, eligible, vax, last.dose.time, at, start,
                                      season.length = days_per_year) {
  interval <- dose_row$interval[[1]]
  annual <- isTRUE(dose_row$annual[[1]])
  if (annual) {
    seasons_elapsed <- max(0, floor((at - start) / season.length))
    season_start <- start + seasons_elapsed * season.length
    which(eligible & at >= season_start & (is.na(last.dose.time) | last.dose.time < season_start))
  } else if (dose_i == 1) {
    which(eligible & vax == 0 & at >= start)
  } else {
    which(eligible & vax == dose_i - 1 & !is.na(last.dose.time) & (at - last.dose.time >= interval) & at >= start)
  }
}

# The strategy set (GlobalMix India #75, decided October 25, 2026). A strategy
# is an ordering of the willing eligibles of the product's age window (#74):
# the first n_take of its order get today's doses, and ties are broken at
# random everywhere.
#   age_oldest        continuous age, oldest first
#   age_transmission  five-year age bands, the band with the highest mean
#                     degree score among the willing eligibles on the day the
#                     campaign opens first (fixed then), at random within a band:
#                     the best order an age-only rule could pick (#58 item 9)
#   random            a uniform permutation
#   degree            the degree score (vax.degree.score), ranked on the day the
#                     campaign opens when vax.degree.frozen, and with fixed
#                     per-person rank noise calibrated to the Spearman
#                     correlation vax.rank.spearman (the precision arm)
#   acquaintance      everyone alive names one of their frequent relationships
#                     at random, and the nominees are taken in the order named:
#                     successive sampling of the eligibles in proportion to
#                     w_j = sum over j's frequent partners i of 1 / degree_i,
#                     redrawn each day; the eligibles nobody can name come last,
#                     at random
#   proxy_obs         a score from traits a program could record (age band, sex,
#                     occupation class, roster household size), with
#                     coefficients fit once to the GlobalMix diary
#                     (vax.proxy.coef), fixed on the day the campaign opens
#   hybrid            two queues from the same willing pool: a reserve of
#                     vax.hybrid.f x the willing eligibles aged vax.hybrid.age
#                     and over (60) on the day the campaign opens, released in
#                     step with the budget, oldest first, and the rest by the
#                     degree score over everyone not taken by the reserve
#   cocoon            the adults vax.cocoon.age.min to .age.max (18-49) of the
#                     households with a living infant under
#                     vax.cocoon.infant.months (6) that day, at random, then
#                     everyone else oldest first (#76); the tier is read each
#                     day, so a newborn brings its household in
# Retired by #75 (retired_strategies): age (renamed age_oldest), bridge,
# age_then_degree, proxy_hh, proxy_fit, and the hybrids age_plus_proxy and
# age_plus_degree; vax.target.uptake gave way to the rank noise and
# vax.hybrid.frac to vax.hybrid.f. Retired by #76: the infant arms, which drew
# the infant products from the adult budget; the infant products are now a dose
# line of their own (vax.infant.program) beside any arm.
vax_strategies <- c("age_oldest", "age_transmission", "random", "degree", "acquaintance", "proxy_obs",
                    "hybrid", "cocoon")
retired_strategies <- c(
  age = "age_oldest, the same ordering under its #75 name",
  bridge = "nothing: on the relationship networks every bridge measure is non-household degree under another name (#75)",
  age_then_degree = "nothing: under gradual rollout it was the age arm (#58 item 13)",
  proxy_hh = "proxy_obs, fit once on traits a program could record (#75)",
  proxy_fit = "proxy_obs, fit once on traits a program could record (#75)",
  age_plus_proxy = "hybrid, with the reserve defined on the 60+ (#75)",
  age_plus_degree = "hybrid, with the reserve defined on the 60+ (#75)",
  infant_hh = "cocoon, the adults 18-49 of the households with an infant under 6 months, then oldest first (#76)",
  infant_direct = paste("vax.infant.program 'nirsevimab' or 'maternal' beside an adult arm such as age_oldest:",
                        "the infant products are a dose line of their own (#76)"),
  infant_direct_cocoon = "vax.infant.program 'nirsevimab' or 'maternal' with vax.strategy 'cocoon' (#76)"
)

allocation_strategy <- function(idsElig, n_take, vax.strategy, age, score = NULL, n_reserve = 0L,
                                reserve.age = 60) {
  nElig <- length(idsElig)
  if (nElig == 0 || n_take <= 0) return(integer(0))
  by_desc <- function(ids, s) ids[order(-s, stats::runif(length(ids)))]
  need_score <- c("degree", "age_transmission", "proxy_obs", "acquaintance", "hybrid", "cocoon")
  if (vax.strategy %in% need_score && is.null(score)) {
    stop("vax.strategy '", vax.strategy, "' requires a score", call. = FALSE)
  }

  if (vax.strategy == "hybrid") {
    old <- idsElig[age[idsElig] >= reserve.age]
    ids_age <- utils::head(by_desc(old, age[old]), max(0L, min(n_reserve, n_take)))
    rest <- setdiff(idsElig, ids_age)
    out <- c(ids_age, utils::head(by_desc(rest, score[rest]), n_take - length(ids_age)))
    attr(out, "n_reserve") <- length(ids_age)
    return(out)
  }

  ids_ranked <- switch(vax.strategy,
    age_oldest = by_desc(idsElig, age[idsElig]),
    random = idsElig[sample.int(nElig)],
    degree = , age_transmission = , proxy_obs = by_desc(idsElig, score[idsElig]),
    acquaintance = {
      # Exponential keys: the order of a weighted draw without replacement
      w <- score[idsElig]
      key <- ifelse(w > 0, stats::rexp(nElig) / w, Inf)
      idsElig[order(key, stats::runif(nElig))]
    },
    cocoon = {
      # the tier (score 1) at random, then the rest oldest first
      in_tier <- score[idsElig] == 1
      idsElig[order(!in_tier, ifelse(in_tier, 0, -age[idsElig]), stats::runif(nElig))]
    },
    stop("Unknown vax.strategy: ", vax.strategy))
  ids_ranked[seq_len(min(n_take, nElig))]
}

# What each person is ranked on when the campaign opens (#75), kept in the
# campaign by unique id, so a person keeps their value as nodes are removed and
# added, and a person born afterwards ranks last: the degree score of a frozen
# degree or hybrid arm, the five-year band means of age_transmission, the
# proxy_obs score, the hybrid's reserve, and the precision arm's noise.
vax_campaign_rank <- function(dat, camp, vax.strategy, eligible, degree_rank, age) {
  ids <- which(get_attr(dat, "active") == 1)
  pool <- which(eligible)
  camp$uid <- get_attr(dat, "unique_id")[ids]
  camp$frozen <- isTRUE(as.logical(param_or(dat, "vax.degree.frozen", FALSE)))
  if (vax.strategy %in% c("degree", "hybrid") && camp$frozen) camp$score <- degree_rank[ids]
  if (vax.strategy == "degree") {
    rho <- param_or(dat, "vax.rank.spearman", 1)
    camp$sigma <- 0
    camp$spearman <- 1
    if (rho < 1 && length(pool) > 1) {
      # The noise is drawn on a side stream, so the arms of a cell keep a common
      # random stream through the opening day and every precision level of a
      # simulation shares the same noise
      camp$ref <- sort(degree_rank[pool])
      z <- with_side_rng(stats::rnorm(length(ids)))
      camp$z <- z
      ns <- normal_score(degree_rank[pool], camp$ref)
      zp <- z[match(pool, ids)]
      camp$sigma <- calibrate_rank_noise(ns, zp, rho)
      camp$spearman <- spearman(ns, ns + camp$sigma * zp)
    }
  }
  if (vax.strategy == "age_transmission") {
    camp$band_mean <- tapply(degree_rank[pool], age_band5(age[pool]), mean)
  }
  if (vax.strategy == "proxy_obs") camp$proxy <- proxy_obs_score(dat)[ids]
  if (vax.strategy == "hybrid") {
    f <- param_or(dat, "vax.hybrid.f", 0)
    a0 <- param_or(dat, "vax.hybrid.age", 60)
    camp$reserve <- min(camp$budget, as.integer(round(f * sum(eligible & age >= a0))))
    camp$reserve_given <- 0L
  }
  camp$ranked <- TRUE
  camp
}

# The score today's ordering uses (#75); NULL for the arms that rank on age or
# at random. Frozen values are looked up by unique id; NA (born since the
# campaign opened) ranks last.
strategy_score <- function(dat, camp, vax.strategy, degree_rank, age, uid) {
  if (is.null(camp)) return(NULL)
  m <- match(uid, camp$uid)
  frozen_value <- function(x) {
    v <- x[m]
    v[is.na(v)] <- -Inf
    v
  }
  switch(vax.strategy,
    degree = {
      s <- if (camp$frozen) frozen_value(camp$score) else degree_rank
      if (camp$sigma > 0) {
        z <- camp$z[m]
        z[is.na(z)] <- 0
        s <- normal_score(s, camp$ref) + camp$sigma * z
      }
      s
    },
    hybrid = if (camp$frozen) frozen_value(camp$score) else degree_rank,
    age_transmission = {
      v <- unname(camp$band_mean[as.character(age_band5(age))])
      v[is.na(v)] <- -Inf
      v
    },
    proxy_obs = frozen_value(camp$proxy),
    acquaintance = acquaintance_weight(dat),
    NULL)
}

# The hybrid's reserve doses for today: the reserve is released in step with the
# budget, so by day t of the campaign floor(reserve x allowed / budget) of them
# may have gone to the 60+; a day when the 60+ queue runs short gives the rest to
# the degree queue, and the reserve catches up later if it can.
hybrid_quota <- function(camp, at) {
  if (is.null(camp$reserve) || camp$budget == 0) return(0L)
  target <- floor(camp$reserve * vax_allowed(camp, at) / camp$budget)
  as.integer(max(0, target - camp$reserve_given))
}

# Five-year age bands (the lower limit, 75 and over in one band)
age_band5 <- function(age) pmin(floor(age / 5), 15) * 5

# Normal score of x on the reference distribution ref (sorted): the mid-rank
# of x among ref, as a share, through the normal quantile. A block of tied
# values takes the score of its mid-rank.
normal_score <- function(x, ref) {
  n <- length(ref)
  p <- (findInterval(x, ref, left.open = TRUE) + findInterval(x, ref)) / (2 * n)
  stats::qnorm(pmin(pmax(p, 0.5 / n), 1 - 0.5 / n))
}

spearman <- function(x, y) stats::cor(rank(x), rank(y))

# The noise scale that gives the target Spearman correlation between a score
# and its noisy version, ns + sigma z, over one population, by bisection with
# the noise fixed (the precision arm, #75). The correlation falls as sigma grows.
calibrate_rank_noise <- function(ns, z, rho, tol = 1e-4) {
  if (rho >= 1) return(0)
  g <- function(s) spearman(ns, ns + s * z) - rho
  hi <- 1
  while (g(hi) > 0 && hi < 1e4) hi <- hi * 2
  lo <- 0
  for (i in 1:60) {
    mid <- (lo + hi) / 2
    if (g(mid) > 0) lo <- mid else hi <- mid
    if (hi - lo < tol) break
  }
  (lo + hi) / 2
}

# The acquaintance weight of every node (#75): the chance, up to a constant,
# that a person drawn at random among the living names this one when naming a
# frequent relationship at random, sum over frequent partners i of 1 / degree_i.
# The frequent relationships are the contacts a person can name: the infrequent
# ones are rarely met, and new contacts are strangers.
acquaintance_weight <- function(dat) {
  if (!activation_on(dat)) {
    stop("the acquaintance strategy names a frequent relationship, which needs the relationship model (#75)",
         call. = FALSE)
  }
  n <- length(get_attr(dat, "active"))
  w <- numeric(n)
  e <- dat$run$el[[rel_layers[["freq"]]]]
  if (is.null(e) || nrow(e) == 0) return(w)
  deg <- tabulate(e, nbins = n)
  s <- rowsum(c(1 / deg[e[, 2]], 1 / deg[e[, 1]]), c(e[, 1], e[, 2]))
  w[as.integer(rownames(s))] <- s[, 1]
  w
}

# The proxy_obs score (#75): the linear predictor of a log-linear fit of daily
# non-household contacts on traits a program could record, from the GlobalMix
# diary (vax.proxy.coef, named as stats::glm() names them: "(Intercept)",
# "band<age.grp>", "sex<F|M>", "occ<class>" and "hh_size"; a level with no
# coefficient is the reference). Household size counts the living members.
proxy_obs_score <- function(dat) {
  coef <- get_param(dat, "vax.proxy.coef", override.null.error = TRUE)
  check_proxy_coef(dat, coef)
  active <- get_attr(dat, "active")
  hh <- get_attr(dat, "hh.ids")
  sex <- get_attr(dat, "sex", override.null.error = TRUE)
  occ <- get_attr(dat, "occ", override.null.error = TRUE)
  if (is.null(sex) || is.null(occ)) {
    stop("proxy_obs reads sex and occupation class, which this population does not carry (#61)", call. = FALSE)
  }
  counts <- table(hh[active == 1])
  hh_size <- as.numeric(counts[match(hh, names(counts))])
  hh_size[is.na(hh_size)] <- 1
  term <- function(nm) {
    v <- coef[nm]
    v[is.na(v)] <- 0
    unname(v)
  }
  unname(coef[["(Intercept)"]] + term(paste0("band", get_attr(dat, "age.grp"))) + term(paste0("sex", sex)) +
           term(paste0("occ", occ)) + coef[["hh_size"]] * hh_size)
}

check_proxy_coef <- function(dat, coef) {
  known <- c("(Intercept)", "hh_size", paste0("band", param_or(dat, "age.grps", character(0))),
             paste0("sex", sex_levels), paste0("occ", occ_levels))
  if (is.null(coef) || !is.numeric(coef) || is.null(names(coef)) || anyNA(coef) ||
      !all(c("(Intercept)", "hh_size") %in% names(coef)) || !all(names(coef) %in% known)) {
    stop("vax.proxy.coef must be named coefficients among ", toString(known),
         ", with an intercept and hh_size (#75)", call. = FALSE)
  }
  invisible(TRUE)
}

# The strategy and its knobs, checked at initialization (#75)
check_vax_strategy <- function(dat) {
  s <- get_param(dat, "vax.strategy", override.null.error = TRUE)
  if (is.null(s)) return(invisible(TRUE))
  if (s %in% names(retired_strategies)) {
    stop("vax.strategy '", s, "' was retired; use ", retired_strategies[[s]], call. = FALSE)
  }
  if (!s %in% vax_strategies) {
    stop("unknown vax.strategy '", s, "'; the strategies are ", toString(vax_strategies), call. = FALSE)
  }
  rho <- param_or(dat, "vax.rank.spearman", 1)
  if (length(rho) != 1 || is.na(rho) || rho <= 0 || rho > 1) {
    stop("vax.rank.spearman must be one number in (0, 1] (#75)", call. = FALSE)
  }
  f <- param_or(dat, "vax.hybrid.f", 0)
  if (length(f) != 1 || is.na(f) || f < 0 || f > 1) stop("vax.hybrid.f must be one number in [0, 1] (#75)", call. = FALSE)
  a0 <- param_or(dat, "vax.hybrid.age", 60)
  if (length(a0) != 1 || is.na(a0) || a0 < 0) stop("vax.hybrid.age must be one age in years (#75)", call. = FALSE)
  fr <- param_or(dat, "vax.degree.frozen", FALSE)
  if (length(fr) != 1 || is.na(as.logical(fr))) stop("vax.degree.frozen must be TRUE or FALSE (#75)", call. = FALSE)
  if (s == "proxy_obs") check_proxy_coef(dat, get_param(dat, "vax.proxy.coef", override.null.error = TRUE))
  if (s == "acquaintance" && !activation_on(dat)) {
    stop("the acquaintance strategy names a frequent relationship, which needs the relationship model (#75)",
         call. = FALSE)
  }
  cmin <- param_or(dat, "vax.cocoon.age.min", 18)
  cmax <- param_or(dat, "vax.cocoon.age.max", 50)
  cm <- param_or(dat, "vax.cocoon.infant.months", 6)
  if (length(cmin) != 1 || length(cmax) != 1 || anyNA(c(cmin, cmax)) || cmin < 0 || cmin >= cmax) {
    stop("vax.cocoon.age.min and .age.max must be two ages in years, the first below the second (#76)", call. = FALSE)
  }
  if (length(cm) != 1 || is.na(cm) || cm <= 0 || cm > 12) {
    stop("vax.cocoon.infant.months must be one number of months in (0, 12] (#76)", call. = FALSE)
  }
  if (s == "cocoon") {
    # under licensure no RSV product covers the tier (the older-adult vaccine is
    # 60 and over), so the cocoon is run on the all-ages product
    pa <- get_pathogen(dat)
    amin <- get_param(dat, paste0("vax.", pa, ".age.min"), override.null.error = TRUE)
    amax <- get_param(dat, paste0("vax.", pa, ".age.max"), override.null.error = TRUE)
    if (!is.null(amin) && !is.null(amax) && (amin >= cmax || amax <= cmin)) {
      stop("the cocoon strategy doses adults ", cmin, " to ", cmax, ", whom the ", pa, " product (ages ", amin,
           " to ", amax, ") does not cover; run it on a product whose age window includes them (#76)",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

# Vaccine protection (GlobalMix India #73). A node's protection comes from its
# most recent dose: vax is the dose count and last.dose.time its step, and the
# dose's row of vax.schedule holds, per outcome, a peak marginal VE, a half-life
# and a floor, each one value or one per vaccine age group. Protection is 0
# until ve.delay days after the dose (the onset delay; 0 for a passive
# product), the peak at ve.delay, and from there
#   VE(t) = floor + (peak - floor) x 0.5^((t - ve.delay) / half-life).
# ve.delay used to postpone only the waning, so a dose protected at its peak
# from the next day.
#
# vax_protection() gives the marginal VE against one outcome. compute_ve()
# gives the relative risk the model draws with at one stage: per contact for a
# susceptible ("infect", 1 - VE_infect), for the symptomatic draw given
# infection ("clinical", (1 - VE_clinical) / (1 - VE_infect)), for the
# admission and death draws given symptomatic disease ("hosp" and "death",
# (1 - VE_hosp) / (1 - VE_clinical) and (1 - VE_death) / (1 - VE_clinical)),
# and on the onward transmission probability of an infector ("infectious",
# 1 - VE_infectious). The VE inputs are therefore marginal, as studies report
# them; the stage relative risks used to be multiplied along the chain, so a
# 0.40 "death VE" on top of the earlier stages was a marginal 0.71.
#
# An infant product (#76; iprog, from infant_prog_state()) protects from its
# own profile (passive_ve()), and a person holding both takes the larger
# protection outcome by outcome, which keeps the stages ordered.
vax_protection <- function(at, ids, vax, last.dose.time, vax.schedule, vax.age.group, outcome, iprog = NULL) {
  ve <- rep(0, length(ids))
  if (length(ids) == 0) return(ve)
  if (!is.null(iprog)) {
    t <- at - iprog$time[ids]
    has <- !is.na(t)
    if (any(has)) ve[has] <- passive_ve(t[has], iprog$pars, outcome)
  }
  if (is.null(vax.schedule)) return(ve)
  vax_ids <- vax[ids]
  for (dose_i in seq_len(nrow(vax.schedule))) {
    pos <- which(vax_ids == dose_i)
    if (length(pos) == 0) next
    ids_this <- ids[pos]
    t <- at - last.dose.time[ids_this]
    if (any(is.na(t))) stop("Missing dose time for vaccinated individuals in vax_protection().")
    ag <- vax.age.group[ids_this]
    peak <- get_schedule_value(vax.schedule, paste0("ve.peak.", outcome), dose_i, ag)
    hl <- get_schedule_value(vax.schedule, paste0("ve.halflife.", outcome), dose_i, ag)
    fl <- get_schedule_value(vax.schedule, paste0("ve.floor.", outcome), dose_i, ag)
    delay <- get_schedule_value(vax.schedule, "ve.delay", dose_i, ag)
    v <- fl + (peak - fl) * 0.5^(pmax(0, t - delay) / hl)
    v[t < delay] <- 0
    ve[pos] <- pmax(ve[pos], v)
  }
  ve
}

#' @rdname moduleset-gmc19
#' @export
compute_ve <- function(at, ids, vax,
                       vax.age.group,
                       last.dose.time = NULL,
                       vax.schedule,
                       outcome = c("infect", "clinical", "hosp", "death", "infectious"),
                       iprog = NULL) {
  outcome <- match.arg(outcome)
  P <- function(o) vax_protection(at, ids, vax, last.dose.time, vax.schedule, vax.age.group, o, iprog)
  ve <- P(outcome)
  rr <- switch(outcome,
               infect = 1 - ve,
               infectious = 1 - ve,
               clinical = cond_rr(ve, P("infect")),
               hosp = cond_rr(ve, P("clinical")),
               death = cond_rr(ve, P("clinical")))
  time.since.last.dose <- rep(0, length(ids))
  if (length(ids) > 0) {
    dosed <- vax[ids] >= 1
    time.since.last.dose[dosed] <- at - last.dose.time[ids][dosed]
  }
  list(rr = rr, ve = ve, time.since.last.dose = time.since.last.dose)
}

# The schedule's VE profiles must be probabilities, wane over a positive
# half-life, and not fall along the stages at the peak or at the floor (or a
# conditional relative risk would exceed 1 and be capped): VE_infect <=
# VE_clinical <= VE_hosp and VE_clinical <= VE_death, per dose and vaccine age
# group. Checked at initialization for the run's pathogen.
check_vax_schedule <- function(dat) {
  pa <- get_pathogen(dat)
  if (is.null(get_param(dat, paste0("vax.", pa, ".n.doses"), override.null.error = TRUE))) {
    return(invisible(TRUE))
  }
  vs <- build_vax_schedule(dat, pa)
  n_age <- length(get_param(dat, "vax.age.breaks")) - 1
  for (dose_i in seq_len(nrow(vs))) {
    for (o in vax_ve_outcomes) {
      for (k in c("peak", "floor")) {
        v <- vs[[paste0("ve.", k, ".", o)]][[dose_i]]
        if (anyNA(v) || any(v < 0 | v > 1)) {
          stop("vax.", pa, ".ve.", k, ".", o, ".dose", dose_i, " must be in [0, 1]", call. = FALSE)
        }
      }
      hl <- vs[[paste0("ve.halflife.", o)]][[dose_i]]
      if (anyNA(hl) || any(hl <= 0)) {
        stop("vax.", pa, ".ve.halflife.", o, ".dose", dose_i, " must be positive", call. = FALSE)
      }
    }
    delay <- vs$ve.delay[[dose_i]]
    if (length(delay) != 1 || is.na(delay) || delay < 0) {
      stop("vax.", pa, ".ve.delay.dose", dose_i, " must be one non-negative number of days", call. = FALSE)
    }
    for (k in c("peak", "floor")) {
      v <- function(o) rep_len(vs[[paste0("ve.", k, ".", o)]][[dose_i]], n_age)
      bad <- v("infect") > v("clinical") | v("clinical") > v("hosp") | v("clinical") > v("death")
      if (any(bad)) {
        stop("vaccine efficacy must not fall along the stages (infect <= clinical <= hosp, and ",
             "clinical <= death): vax.", pa, ".ve.", k, ".*.dose", dose_i, " do in vaccine age group(s) ",
             toString(which(bad)), call. = FALSE)
      }
    }
  }
  invisible(TRUE)
}

# The campaign parameters (GlobalMix India #74), checked at initialization with
# the schedule: vax.budget a non-negative share of the living (above 1 for a
# multi-dose product), vax.delivery.days a
# number of days of at least 1, vax.start a day of at least 1, vax.willing one
# probability or one per vaccine age group, and the product's age window
# vax.<pathogen>.age.min < .age.max in years.
check_vax_campaign <- function(dat) {
  check_vax_strategy(dat)
  pa <- get_pathogen(dat)
  if (is.null(get_param(dat, paste0("vax.", pa, ".n.doses"), override.null.error = TRUE))) {
    return(invisible(TRUE))
  }
  one <- function(name, lo, hi = Inf) {
    v <- get_param(dat, name, override.null.error = TRUE)
    if (is.null(v) || length(v) != 1 || is.na(v) || v < lo || v > hi) {
      stop(name, " must be one number in [", lo, ", ", hi, "] (GlobalMix India #74)", call. = FALSE)
    }
    v
  }
  one("vax.budget", 0)
  one("vax.delivery.days", 1)
  one("vax.start", 1)
  n_age <- length(get_param(dat, "vax.age.breaks")) - 1
  w <- get_param(dat, "vax.willing", override.null.error = TRUE)
  if (is.null(w) || !length(w) %in% c(1, n_age) || anyNA(w) || any(w < 0 | w > 1)) {
    stop("vax.willing must be one probability, or one per vaccine age group (", n_age,
         "), in [0, 1] (GlobalMix India #74)", call. = FALSE)
  }
  amin <- one(paste0("vax.", pa, ".age.min"), 0)
  amax <- one(paste0("vax.", pa, ".age.max"), 0)
  if (amin >= amax) stop("vax.", pa, ".age.min must be below vax.", pa, ".age.max (GlobalMix India #74)", call. = FALSE)
  invisible(TRUE)
}

# Outcomes by vaccine age group and vaccination status (GlobalMix India #73
# acceptance), with vax.diagnostics = TRUE: per step, vd.<what>.<g>.<v> for
# what in n (the living), sus (the susceptible), inf, sympt, hosp and death,
# g the vaccine age group and v 1 for the dosed and 0 for the rest. Realized
# marginal VE within an age group is read off them: against infection from the
# infections per susceptible person-day, against symptomatic disease,
# admission and death from the events per living person.
vax_diagnostics <- function(dat) isTRUE(as.logical(param_or(dat, "vax.diagnostics", FALSE)))

vd_flows <- function(dat, at, what, ids) {
  n_age <- length(get_param(dat, "vax.age.breaks")) - 1
  g <- vax_age_group_for(dat)[ids]
  v <- get_attr(dat, "vax")[ids] >= 1
  for (k in seq_len(n_age)) {
    dat <- set_epi(dat, paste0("vd.", what, ".", k, ".1"), sum(g == k & v, na.rm = TRUE), at = at)
    dat <- set_epi(dat, paste0("vd.", what, ".", k, ".0"), sum(g == k & !v, na.rm = TRUE), at = at)
  }
  dat
}

# The day's exposures (discordant contacts) of susceptibles, by vaccine age
# group and vaccination status: vd.exp.<g>.<v> their number, vd.einf.<g>.<v>
# the infections they imply (each susceptible infected with probability
# 1 - prod(1 - p) over its contacts, before any draw) and vd.einf0.<g>.<v> the
# same with the susceptible's own vaccine protection against infection taken
# out. The ratio einf / einf0 among the dosed is the reduction in infection the
# model's own exposures imply, which a per-person-day hazard ratio measured in
# a saturating epidemic falls short of (several exposures on one day infect at
# most once), so the cohort check judges against it rather than against the
# per-contact VE.
vd_exposure <- function(dat, at, dels, vax, last.dose.time, vax.schedule) {
  n_age <- length(get_param(dat, "vax.age.breaks")) - 1
  sus <- unlist(lapply(dels, function(d) d$sus))
  zero <- function(what) {
    for (k in seq_len(n_age)) for (v in c("1", "0")) dat <<- set_epi(dat, paste0("vd.", what, ".", k, ".", v), 0, at = at)
  }
  if (length(sus) == 0) {
    for (what in c("exp", "einf", "einf0")) zero(what)
    return(dat)
  }
  p <- unlist(lapply(dels, function(d) d$finalProb))
  tp <- unlist(lapply(dels, function(d) d$transProb))
  acts <- unlist(lapply(dels, function(d) d$acts))
  rr <- compute_ve(at = at, ids = sus, vax = vax, last.dose.time = last.dose.time, vax.schedule = vax.schedule,
                   outcome = "infect", vax.age.group = vax_age_group_for(dat), iprog = infant_prog_state(dat))$rr
  p0 <- 1 - (1 - pmin(1, tp / pmax(rr, 1e-12)))^acts
  lq <- rowsum(log1p(-pmin(p, 1 - 1e-12)), sus)
  lq0 <- rowsum(log1p(-pmin(p0, 1 - 1e-12)), sus)
  ids <- as.integer(rownames(lq))
  cnt <- rowsum(rep(1, length(sus)), sus)[, 1]
  g <- vax_age_group_for(dat)[ids]
  v <- vax[ids] >= 1
  for (k in seq_len(n_age)) {
    for (vv in c(TRUE, FALSE)) {
      sel <- g == k & v == vv
      tag <- paste0(".", k, ".", as.integer(vv))
      dat <- set_epi(dat, paste0("vd.exp", tag), sum(cnt[sel]), at = at)
      dat <- set_epi(dat, paste0("vd.einf", tag), sum(1 - exp(lq[sel, 1])), at = at)
      dat <- set_epi(dat, paste0("vd.einf0", tag), sum(1 - exp(lq0[sel, 1])), at = at)
    }
  }
  dat
}

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
      labels = seq_len(length(vax.age.breaks) - 1)
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

# Build the dose-indexed vaccine schedule from the parameters
# vax.<disease>.<field>.dose<i>, one row per dose, with the fields in
# vax_dose_fields (R/params.R). compute_ve() reads the VE columns, and
# vax_general() the administration columns (start, rate, interval, annual).
build_vax_schedule <- function(dat, disease = "covid") {
  prefix <- paste0("vax.", disease, ".")
  n_doses <- get_param(dat, paste0(prefix, "n.doses"))
  if (length(n_doses) != 1 || is.na(n_doses) || n_doses < 1) {
    stop(prefix, "n.doses must be one dose count of at least 1", call. = FALSE)
  }
  br <- get_param(dat, "vax.age.breaks", override.null.error = TRUE)
  # the vaccine age groups; without vax.age.breaks (a bare parameter list) a
  # VE field may hold any length
  n_age <- if (is.null(br)) NA_integer_ else length(br) - 1
  vax.schedule <- data.frame(dose = seq_len(n_doses))
  for (field in vax_dose_fields) {
    values <- lapply(seq_len(n_doses), function(i) {
      v <- get_param(dat, paste0(prefix, field, ".dose", i))
      n <- length(v)
      ok <- if (field %in% vax_dose_scalar) n == 1 else if (is.na(n_age)) n >= 1 else n %in% c(1, n_age)
      if (!ok) {
        stop(prefix, field, ".dose", i, " has ", n, " value(s); it needs ",
             if (field %in% vax_dose_scalar) "one" else if (field %in% vax_dose_by_age) paste("one per vaccine age group,", n_age)
             else paste("one, or one per vaccine age group,", n_age), " (GlobalMix India #73)", call. = FALSE)
      }
      v
    })
    # Every field is a list column: one value per dose, or one per vaccine age
    # group (start and rate always; the VE fields optionally, #73). sapply()
    # used to turn a one-dose schedule with an age vector into five rows.
    vax.schedule[[field]] <- I(values)
  }
  if (nrow(vax.schedule) != n_doses) stop("the vaccine schedule has ", nrow(vax.schedule), " rows for ", n_doses, " doses")
  vax.schedule
}

