
#' Model Extensions to EpiModel for the Network-Based Modeling of COVID-19
#'
#' \tabular{ll}{
#'    Package: \tab EpiModelCOVID\cr
#'    Type: \tab Package\cr
#'    Version: \tab 1.0.0\cr
#'    Date: \tab 2020-05-17\cr
#'    License: \tab GPL-3\cr
#'    LazyLoad: \tab yes\cr
#' }
#'
#' @name EpiModelCOVID-package
#' @aliases EpiModelCOVID
#'
#' @import EpiModel tergmLite ergm network
#' @importFrom stats rbinom rgeom rmultinom rpois runif simulate rnbinom plogis predict
#'
#' @docType package
#' @keywords package
#'
NULL


#' @title EpiModel Module Set for COVID Cruise Ship Model
#'
#' @description This set of functions is associated with the EpiModel study of
#'              COVID on cruise ship environments.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netest}}.
#' @param param An \code{EpiModel} object of class \code{\link{param_msm}}.
#' @param init An \code{EpiModel} object of class \code{\link{init_msm}}.
#' @param control An \code{EpiModel} object of class \code{\link{control_msm}}.
#' @param s Simulation number, used for restarting dependent simulations.
#' @param dat Master list object of network models.
#' @param at Current time step.
#'
#' @name moduleset-ship
#'
NULL
