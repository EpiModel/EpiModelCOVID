
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
#' @return
#' This function returns `dat` after updating the nodal attributes `age` and
#' `age.grp`.
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
  age[active == 1] <- age[active == 1] + time.unit / days_per_year
  
  
  age.breaks <- get_param(dat, "age.breaks") 
  age.grps <- get_param(dat, "age.grps") 

  
  age.grp[active == 1] <- cut(
    age[active == 1],
    age.breaks,
    labels = age.grps,
    right = FALSE
  ) |> as.character()

  ## Output
  # Set Attributes
  dat <- set_attr(dat, "age.grp", age.grp)
  dat <- set_attr(dat, "age", age)

  return(dat)
}

# Days in a year (GlobalMix India #61): ages advance by time.unit / days_per_year
# each step, and the seasonal forcing and the vaccination season recur every
# days_per_year days. The CorporateMix code used 364, 52 weeks, which made a
# simulated year a day short of a calendar year.
days_per_year <- 365
