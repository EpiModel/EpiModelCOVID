
#' @title Aging Module
#'
#' @description Module for aging over time for active nodes in the population.
#'
#' @param dat Main data object of class `netsim_dat` containing networks,
#'            individual-level attributes, and summary statistics.
#' @param at Current time step.
#'
#' @details
#' Ages advance by `time.unit / 365` years each step, one day on the daily
#' model (GlobalMix India #61).
#'
#' On the relationship design the infant flag `is_infant` is recomputed from age
#' every step (`age < infant.age.max`; GlobalMix India #64).
#'
#' On the relationship design, a person whose age band or occupation class
#' changes on the step has their relationship traits and ties moved to the new
#' band's cross-section ([occ_transitions()], [cross_band_tops()], GlobalMix
#' India #64).
#'
#' @return
#' This function returns `dat` after updating the nodal attributes `age`,
#' `age.grp` and `is_infant`.
#'
#' @export
#'
aging <- function(dat, at) {
  # if (at>10) browser()
  
  ## Input
  # Attributes
  age     <- get_attr(dat, "age")
  active  <- get_attr(dat, "active")
  age.grp <- get_attr(dat, "age.grp") 

  # Parameters
  time.unit <- get_param(dat, "time.unit")

  ## Process
  grp_before <- age.grp
  age_before <- age
  age[active == 1] <- age[active == 1] + time.unit / days_per_year
  
  
  age.breaks <- get_param(dat, "age.breaks") 
  age.grps <- get_param(dat, "age.grps") 

  
  age.grp[active == 1] <- cut(
    age[active == 1],
    age.breaks,
    labels = age.grps,
    right = FALSE
  ) |> as.character()

  # The infant flag follows age (GlobalMix India #64, item 1 of #65). It was
  # set only at initialization and birth, so an initial infant kept it up to
  # about age 2, with the infant severity, vaccine eligibility and counts. The
  # legacy four-layer model keeps the old flag, so its runs still reproduce the
  # Redo 08 runs.
  is_infant <- get_attr(dat, "is_infant", override.null.error = TRUE)
  if (!is.null(is_infant) && activation_on(dat)) {
    infant.age.max <- get_param(dat, "infant.age.max", override.null.error = TRUE)
    if (is.null(infant.age.max)) infant.age.max <- 1
    is_infant[active == 1] <- age[active == 1] < infant.age.max
    dat <- set_attr(dat, "is_infant", is_infant)
  }

  ## Output
  # Set Attributes
  dat <- set_attr(dat, "age.grp", age.grp)
  dat <- set_attr(dat, "age", age)

  # People who crossed a band top or changed occupation class on this step move
  # to the new band's and class's relationship traits and cross-section
  # (GlobalMix India #64; occ_transitions(), cross_band_tops())
  if (activation_on(dat)) {
    ct <- occ_transitions(dat, active, age_before, age)
    dat <- ct$dat
    crossed <- which(active == 1 & !is.na(grp_before) & age.grp != grp_before)
    ids <- sort(unique(c(crossed, ct$changed)))
    dat <- cross_band_tops(dat, at, ids, match(grp_before[ids], age.grps),
                           occ_from = if (!is.null(ct$occ_before)) ct$occ_before[ids])
    dat <- set_epi(dat, "n.band.cross", length(crossed), at = at)
    dat <- set_epi(dat, "n.class.change", length(ct$changed), at = at)
  }

  return(dat)
}

# Days in a year (GlobalMix India #61): ages advance by time.unit / days_per_year
# each step, and the seasonal forcing and the vaccination season recur every
# days_per_year days. The CorporateMix code used 364, 52 weeks, which made a
# simulated year a day short of a calendar year.
days_per_year <- 365
