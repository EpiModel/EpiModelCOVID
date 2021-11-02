# EpiModelCOVID

<!-- badges: start -->

[![R-CMD-check](https://github.com/EpiModel/EpiModelCOVID/workflows/R-CMD-check/badge.svg)](https://github.com/EpiModel/EpiModelCOVID/actions)

<!-- badges: end -->

Modules for simulating SARS-CoV-2 transmission dynamics in different epidemiological settings, developed as an extension to our general network-based epidemic modeling platform, [EpiModel](http://epimodel.org).

`EpiModel` and `EpiModelCOVID` use the statistical framework of temporal exponential-family random graph models to fit and simulate models of dynamic contact networks. These [statistical methods](http://onlinelibrary.wiley.com/doi/10.1111/rssb.12014/abstract) have been developed and implemented as open-source software, building on the extensive efforts of the [Statnet](https://statnet.org/) research group to build software tools for the representation, analysis, and visualization of complex network data.

These packages combine these Statnet methods with an individual-based epidemic modeling engine to simulate SARS-Cov-2 transmission over networks, allowing for complex dependencies between the network, epidemiological, and demographic changes in the simulated populations. Readers new to these methods are recommended to consult our [EpiModel](http://epimodel.org) resources, including our main methods paper [Vignette](http://doi.org/10.18637/jss.v084.i08) describing the theory and implementation.

## Installation

You can install `EpiModelCOVID` in R using `remotes`:

    install.packages("EpiModel", dependencies = TRUE)
    remotes::install_github("EpiModel/EpiModelCOVID")

Documentation on using this software package is forthcoming, although limited function documentation is provided within the package and available with the `help(package = "EpiModelCOVID")` command.
