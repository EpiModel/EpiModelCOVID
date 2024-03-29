% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EpiModelCOVID-package.R, R/mod-demog.R,
%   R/mod-infection.R, R/mod-init.R, R/mod-prevalence.R, R/mod-simnet.R
\name{moduleset-corporate}
\alias{moduleset-corporate}
\alias{deaths_covid_corporate}
\alias{arrival_covid_corporate}
\alias{infect_covid_corporate}
\alias{init_covid_corporate}
\alias{prevalence_covid_corporate}
\alias{resim_nets_covid_corporate}
\title{EpiModel Module Set for Corporate Office Model}
\usage{
deaths_covid_corporate(dat, at)

arrival_covid_corporate(dat, at)

infect_covid_corporate(dat, at)

init_covid_corporate(x, param, init, control, s)

prevalence_covid_corporate(dat, at)

resim_nets_covid_corporate(dat, at)
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
             COVID in corporate office environments.
}
