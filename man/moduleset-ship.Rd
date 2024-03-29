% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EpiModelCOVID-package.R, R/mod-demog.R,
%   R/mod-infection.R, R/mod-init.R, R/mod-prevalence.R, R/mod-progress.R,
%   R/mod-simnet.R
\name{moduleset-ship}
\alias{moduleset-ship}
\alias{deaths_covid_ship}
\alias{offload_covid_ship}
\alias{infect_covid_ship}
\alias{init_covid_ship}
\alias{prevalence_covid_ship}
\alias{progress_covid_ship}
\alias{resim_nets_covid_ship}
\title{EpiModel Module Set for COVID Cruise Ship Model}
\usage{
deaths_covid_ship(dat, at)

offload_covid_ship(dat, at)

infect_covid_ship(dat, at)

init_covid_ship(x, param, init, control, s)

prevalence_covid_ship(dat, at)

progress_covid_ship(dat, at)

resim_nets_covid_ship(dat, at)
}
\arguments{
\item{dat}{Main `netsim_dat` class data object of network models.}

\item{at}{Current time step.}

\item{x}{An \code{EpiModel} object of class `netest`.}

\item{param}{An \code{EpiModel} object of class `param.net`.}

\item{init}{An \code{EpiModel} object of class `init.net`.}

\item{control}{An \code{EpiModel} object of class `control.net`.}

\item{s}{Simulation number, used for restarting dependent simulations.}
}
\description{
This set of functions is associated with the EpiModel study of
             COVID on cruise ship environments.
}
