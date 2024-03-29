
#' Model Extensions to EpiModel for the Network-Based Modeling of COVID-19
#'
#' \tabular{ll}{
#'    Package: \tab EpiModelCOVID\cr
#'    Type: \tab Package\cr
#'    Version: \tab 1.2.0\cr
#'    Date: \tab 2022-07-29\cr
#'    License: \tab GPL-3\cr
#'    LazyLoad: \tab yes\cr
#' }
#'
#' @name EpiModelCOVID-package
#' @aliases EpiModelCOVID
#'
#' @import EpiModel ergm network tergm
#' @importFrom networkLite networkLite
#' @importFrom statnet.common NVL
#' @importFrom stats rbinom simulate rpois
#' @importFrom utils tail
#'
#' @docType package
#' @keywords package
#'
NULL


#' @title EpiModel Common Module Set Across COVID Models
#'
#' @description This set of functions is associated with all EpiModel studies of
#'              COVID.
#'
#' @param dat Main `netsim_dat` class data object of network models.
#' @param at Current time step.
#'
#' @name moduleset-common
#'
NULL

#' @title EpiModel Module Set for COVID Cruise Ship Model
#'
#' @description This set of functions is associated with the EpiModel study of
#'              COVID on cruise ship environments.
#'
#' @param x An \code{EpiModel} object of class `netest`.
#' @param param An \code{EpiModel} object of class `param.net`.
#' @param init An \code{EpiModel} object of class `init.net`.
#' @param control An \code{EpiModel} object of class `control.net`.
#' @param s Simulation number, used for restarting dependent simulations.
#' @param dat Main `netsim_dat` class data object of network models.
#' @param at Current time step.
#'
#' @name moduleset-ship
#'
NULL

#' @title EpiModel Module Set for Corporate Office Model
#'
#' @description This set of functions is associated with the EpiModel study of
#'              COVID in corporate office environments.
#'
#' @param x An \code{EpiModel} object of class `netest`.
#' @param param An \code{EpiModel} object of class `param.net`.
#' @param init An \code{EpiModel} object of class `init.net`.
#' @param control An \code{EpiModel} object of class `control.net`.
#' @param s Simulation number, used for restarting dependent simulations.
#' @param dat Main `netsim_dat` class data object of network models.
#' @param at Current time step.
#'
#' @name moduleset-corporate
#'
NULL
