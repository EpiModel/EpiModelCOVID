
#' @title Aging Module
#'
#' @description Module for aging over time for active nodes in the population.
#'
#' @param dat Main data object of class `netsim_dat` containing networks,
#'            individual-level attributes, and summary statistics.
#' @param at Current time step.
#'
#' @return
#' This function returns `dat` after updating the nodal attributes `age` and
#' `sqrt.age`.
#'
#' @export
#'
aging <- function(dat, at) {

  
  ## Input
  # Attributes
  age     <- get_attr(dat, "age")
  active  <- get_attr(dat, "active")
  age.grp <- get_attr(dat, "age.grp") 

  # Parameters
  time.unit <- 1
    #get_param(dat, "time.unit")

  ## Process
  age[active == 1] <- age[active == 1] + time.unit / 364
  
  #netstats <- get_param(dat, "netstats")
  
  age.breaks <-c(0, 10, 20, 30, 40, 60, Inf)
  age.labels <- c("0-9y", "10-19y", "20-29y", "30-39y", "40-59y", "60+y")
  

    #netstats[["demog"]][["age.breaks"]]

  
  age.grp[active == 1] <- cut(
    age[active == 1],
    age.breaks,
    labels = age.labels,
    right = FALSE
  )

  ## Output
  # Set Attributes
  dat <- set_attr(dat, "age.grp", age.grp)
  dat <- set_attr(dat, "age", age)

  return(dat)
}
