
#' @rdname moduleset-corporate
#' @export
resim_nets_covid_corporate <- function(dat, at) {
  
  dat$num.nw <- dat$num.nw - 1 # disregard household nw for this module

  ## Edges correction
  dat <- edges_correct(dat, at)
  
  ## network resimulation
  dat.updates <- NVL(get_control(dat, "dat.updates"), function(dat, ...) dat)
  dat <- dat.updates(dat = dat, at = at, network = 0L)
  for (network in seq_len(dat$num.nw)) {
    dat <- simulate_dat(dat = dat, at = at, network = network)
    dat <- dat.updates(dat = dat, at = at, network = network)
  }

  if (get_control(dat, "cumulative.edgelist")) {
    for (n_network in seq_len(dat$num.nw)) {
      browser()
      dat <- update_cumulative_edgelist(dat, n_network,
                                        get_control(dat, "truncate.el.cuml"))
    }
  }

  dat <- summary_nets(dat, at)
  
  dat$num.nw <- dat$num.nw + 1

  return(dat)
}

#' @rdname moduleset-corporate
#' @export
update_cumulative_edgelist <- function(dat, network, truncate = 0) {
  if (!get_control(dat, "cumulative.edgelist")) {
    return(dat)
  }
  
  el <- get_edgelist(dat, network)
  el_cuml <- get_cumulative_edgelist(dat, network)
  
  el <- tibble::tibble(
    head = get_unique_ids(dat, el[, 1]),
    tail = get_unique_ids(dat, el[, 2]),
    current = TRUE
  )
  
  el_cuml <- dplyr::full_join(el_cuml, el, by = c("head", "tail"))
  
  at <- get_current_timestep(dat)
  
  new_edges <- is.na(el_cuml[["start"]])
  if (any(new_edges)) {
    el_cuml[new_edges, ][["start"]] <- at
  }
  
  terminated_edges <- is.na(el_cuml[["current"]]) & is.na(el_cuml[["stop"]])
  if (any(terminated_edges)) {
    el_cuml[terminated_edges, ][["stop"]] <- at - 1
  }
  
  if (truncate != Inf) {
    rel.age <- at - el_cuml[["stop"]]
    rel.age <- ifelse(is.na(rel.age), 0, rel.age)
    el_cuml <- el_cuml[rel.age <= truncate, ]
  }
  
  dat[["el.cuml"]][[network]] <- el_cuml[, c("head", "tail", "start", "stop")]
  
  return(dat)
}