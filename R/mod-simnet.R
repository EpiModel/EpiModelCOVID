
#' @rdname moduleset-ship
#' @export
resim_nets_covid_ship <- function(dat, at) {

  ## Edges correction
  dat <- edges_correct_covid(dat, at)

  if (at < dat$param$network.lockdown.time) {
    nets <- 1:3
  } else {
    nets <- 4:6
  }

  # Network Resimulation
  for (i in nets) {
    nwparam <- EpiModel::get_nwparam(dat, network = i)
    isTERGM <- ifelse(nwparam$coef.diss$duration > 1, TRUE, FALSE)
    dat <- tergmLite::updateModelTermInputs(dat, network = i)
    if (isTERGM == TRUE) {
      dat$el[[i]] <- tergmLite::simulate_network(p = dat$p[[i]],
                                                 el = dat$el[[i]],
                                                 coef.form = nwparam$coef.form,
                                                 coef.diss = nwparam$coef.diss$coef.adj,
                                                 save.changes = FALSE)
    } else {
      dat$el[[i]] <- tergmLite::simulate_ergm(p = dat$p[[i]],
                                              el = dat$el[[i]],
                                              coef = nwparam$coef.form)
    }
  }

  if (dat$control$save.nwstats == TRUE) {
    dat <- calc_nwstats_covid(dat, at)
  }

  return(dat)
}


#' @rdname moduleset-corporate
#' @export
resim_nets_covid_corporate <- function(dat, at) {

  ## Edges correction
  dat <- edges_correct_covid(dat, at)

  # Network Resimulation
  for (i in 1:length(dat$el)) {
    nwparam <- EpiModel::get_nwparam(dat, network = i)
    isTERGM <- ifelse(nwparam$coef.diss$duration > 1, TRUE, FALSE)
    dat <- tergmLite::updateModelTermInputs(dat, network = i)
    if (isTERGM == TRUE) {
      dat$el[[i]] <- tergmLite::simulate_network(p = dat$p[[i]],
                                                 el = dat$el[[i]],
                                                 coef.form = nwparam$coef.form,
                                                 coef.diss = nwparam$coef.diss$coef.adj,
                                                 save.changes = FALSE)
    } else {
      dat$el[[i]] <- tergmLite::simulate_ergm(p = dat$p[[i]],
                                              el = dat$el[[i]],
                                              coef = nwparam$coef.form)
    }
  }

  if (dat$control$save.nwstats == TRUE) {
    dat <- calc_nwstats_covid(dat, at)
  }

  return(dat)
}


edges_correct_covid <- function(dat, at) {

  old.num <- dat$epi$num[at - 1]
  new.num <- sum(dat$attr$active == 1, na.rm = TRUE)
  adjust <- log(old.num) - log(new.num)

  for (i in 1:length(dat$nwparam)) {
    coef.form1 <- get_nwparam(dat, network = i)$coef.form
    coef.form1[1] <- coef.form1[1] + adjust
    dat$nwparam[[i]]$coef.form <- coef.form1
  }

  return(dat)
}


calc_nwstats_covid <- function(dat, at) {

  for (nw in 1:length(dat$el)) {
    n <- attr(dat$el[[nw]], "n")
    edges <- nrow(dat$el[[nw]])
    meandeg <- round(edges * (2/n), 3)
    concurrent <- round(mean(get_degree(dat$el[[nw]]) > 1), 3)
    mat <- matrix(c(edges, meandeg, concurrent), ncol = 3, nrow = 1)
    if (at == 1) {
      dat$stats$nwstats[[nw]] <- mat
      colnames(dat$stats$nwstats[[nw]]) <- c("edges", "mdeg", "conc")
    }
    if (at > 1) {
      dat$stats$nwstats[[nw]] <- rbind(dat$stats$nwstats[[nw]], mat)
    }
  }

  return(dat)
}



#' @rdname moduleset-contacttrace
#' @export
resim_nets_covid_contacttrace <- function(dat, at) {
  
  ## Edges correction
  dat <- edges_correct_covid(dat, at)
  
  # Network Resimulation
  for (i in 1:length(dat$el)) {
    nwparam <- EpiModel::get_nwparam(dat, network = i)
    isTERGM <- ifelse(nwparam$coef.diss$duration > 1, TRUE, FALSE)
    dat <- tergmLite::updateModelTermInputs(dat, network = i)
    if (isTERGM == TRUE) {
      dat$el[[i]] <- tergmLite::simulate_network(p = dat$p[[i]],
                                                 el = dat$el[[i]],
                                                 coef.form = nwparam$coef.form,
                                                 coef.diss = nwparam$coef.diss$coef.adj,
                                                 save.changes = FALSE)
    } else {
      dat$el[[i]] <- tergmLite::simulate_ergm(p = dat$p[[i]],
                                              el = dat$el[[i]],
                                              coef = nwparam$coef.form)
    }
  }
  
  if (dat$control$save.nwstats == TRUE) {
    dat <- calc_nwstats_covid(dat, at)
  }
  
  for (n_network in seq_len(3)) {
    dat <- update_cumulative_edgelist(dat, n_network)
  }
  
  return(dat)
}


edges_correct_covid <- function(dat, at) {
  
  old.num <- dat$epi$num[at - 1]
  new.num <- sum(dat$attr$active == 1, na.rm = TRUE)
  adjust <- log(old.num) - log(new.num)
  
  for (i in 1:length(dat$nwparam)) {
    coef.form1 <- get_nwparam(dat, network = i)$coef.form
    coef.form1[1] <- coef.form1[1] + adjust
    dat$nwparam[[i]]$coef.form <- coef.form1
  }
  
  return(dat)
}


calc_nwstats_covid <- function(dat, at) {
  
  for (nw in 1:length(dat$el)) {
    n <- attr(dat$el[[nw]], "n")
    edges <- nrow(dat$el[[nw]])
    meandeg <- round(edges * (2/n), 3)
    concurrent <- round(mean(get_degree(dat$el[[nw]]) > 1), 3)
    mat <- matrix(c(edges, meandeg, concurrent), ncol = 3, nrow = 1)
    if (at == 1) {
      dat$stats$nwstats[[nw]] <- mat
      colnames(dat$stats$nwstats[[nw]]) <- c("edges", "mdeg", "conc")
    }
    if (at > 1) {
      dat$stats$nwstats[[nw]] <- rbind(dat$stats$nwstats[[nw]], mat)
    }
  }
  
  return(dat)
}
