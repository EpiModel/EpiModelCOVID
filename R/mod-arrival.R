

#' @rdname moduleset-gmc19
#' @export
arrival <- function(dat, at) {

  ## Births (GlobalMix India #64). birth.rate is the crude birth rate per
  ## person-day: each person alive after this step's deaths has a birth with
  ## that probability, so births do not depend on deaths, and the population
  ## size drifts by the difference of the two. 0 closes the population (the
  ## closed-cohort design, with background mortality also at 0). The relationship
  ## design must set it (init_gmc19() stops otherwise). The legacy four-layer
  ## model sets none and keeps its pre-redo design, in which births replace the
  ## step's deaths, so that a legacy run still reproduces the Redo 08 runs.
  birth.rate <- get_param(dat, "birth.rate", override.null.error = TRUE)
  if (is.null(birth.rate)) {
    nDep <- get_epi(dat, "d.flow", at = at)
    nNew <- if (is.null(nDep) || is.na(nDep)) 0L else as.integer(nDep)
  } else {
    nAlive <- sum(get_attr(dat, "active") == 1)
    nNew <- stats::rbinom(1, nAlive, birth.rate)
  }

  # Update Attr. The new nodes join the networks in nwupdate.FUN (EpiModel's
  # nwupdate.net), which runs after this module, so this module only appends
  # their attributes. On the relationship design each newborn's household is
  # its mother's (place_births(), GlobalMix India #65); on the household clique
  # layer, nwupdate.FUN then keeps that household with the layer's arrival rule
  # (hh_join_birth) and connects the newborn to its co-residents.
  if (nNew > 0) {
    dat <- init_new_nodes_attrs(dat, nNew)
    if (births_placed(dat)) {
      n <- length(get_attr(dat, "active"))
      dat <- place_births(dat, at, n - nNew + seq_len(nNew))
    }
  }

  ## Output
  dat <- set_epi(dat, "nNew", nNew, at = at)
  
  return(dat)
}



init_new_nodes_attrs <- function(dat, n_new) {
  current_timestep <- get_current_timestep(dat)
  # Core attributes (necessary for EpiModel)
  dat <- append_core_attr(dat, current_timestep, n_new)
  def_attrs <- get_default_attrs(dat)
  for (attr_name in names(def_attrs)) {
    dat <- append_attr(dat, attr_name, def_attrs[[attr_name]], n_new)
  }
  dat <- make_computed_attrs(dat, n_new, post_init = TRUE)

  # Every attribute must get a value here. nwupdate.net() would otherwise fill
  # a missing one by drawing from its current distribution, with no error.
  attr_lengths <- lengths(get_attr_list(dat))
  short <- names(attr_lengths)[attr_lengths < attr_lengths[["active"]]]
  if (length(short) > 0) {
    stop("no value set for new nodes for attribute(s) ",
         paste(short, collapse = ", "), "; add a default in get_default_attrs()")
  }
  return(dat)
}

# Birth placement (GlobalMix India #65) -------------------------------------------
#
# On the relationship design, with a population that carries sex, each newborn
# has a mother: a woman aged 15-49 who has not given birth in the last
# birth.interval.min days, drawn with weight
#
#   birth.age.weight[her five-year age group] x exp(birth.hh.coef . x_h),
#
# where x_h is the composition of her household (birth_hh_terms). The age
# weights are the age-specific fertility rates, and the household coefficients
# are calibrated, with an adjustment to the age weights, so that newborns'
# households match the infants' households of the starting population while
# mothers' ages keep the fertility rates' distribution (the project's
# R/A-networks/population/birth_placement.R). The newborn joins its mother's
# household, and the mother's last.birth is the step. Before engine 1.10.0
# newborns joined the household of a random 0-9y member (hh_join_young), which
# drew households in proportion to their young children, so there were no
# first births and child households grew; the legacy four-layer model, and a
# relationship population without sex, keep that rule.

# Household composition terms of birth.hh.coef, in order: members aged 0-9,
# 10-19, women 20-39, men 20-39, 40-59, 60 and over, and no member under 10
birth_hh_terms <- c("n0_9", "n10_19", "f20_39", "m20_39", "n40_59", "n60", "nochild")

# Mothers' ages: five-year groups from 15 to 49, the groups of birth.age.weight
birth_age_breaks <- seq(15, 50, by = 5)

births_placed <- function(dat) {
  activation_on(dat) && has_pop_class(dat)
}

# Checked at initialization: the placement parameters, and a household layer
# whose arrival rule keeps the household place_births() sets
check_birth_placement <- function(dat) {
  if (!births_placed(dat)) return(invisible(TRUE))
  w <- get_param(dat, "birth.age.weight", override.null.error = TRUE)
  if (!(is.numeric(w) && length(w) == length(birth_age_breaks) - 1 && all(is.finite(w)) &&
        all(w >= 0) && any(w > 0))) {
    stop("birth.age.weight must be ", length(birth_age_breaks) - 1, " non-negative weights for the ",
         "mother's age groups 15-19 to 45-49, not all 0 (GlobalMix India #65)", call. = FALSE)
  }
  b <- get_param(dat, "birth.hh.coef", override.null.error = TRUE)
  if (!(is.numeric(b) && identical(names(b), birth_hh_terms) && all(is.finite(b)))) {
    stop("birth.hh.coef must be finite coefficients named ", toString(birth_hh_terms),
         " (GlobalMix India #65)", call. = FALSE)
  }
  f <- get_param(dat, "birth.interval.min", override.null.error = TRUE)
  if (!(is.numeric(f) && length(f) == 1 && is.finite(f) && f >= 0)) {
    stop("birth.interval.min must be one non-negative number of days (GlobalMix India #65)",
         call. = FALSE)
  }
  # Household ids index the composition table by row, so they must be positive
  # whole numbers
  hh.ids <- get_attr(dat, "hh.ids", override.null.error = TRUE)
  if (is.null(hh.ids) || !is.numeric(hh.ids) ||
      any(!is.na(hh.ids) & (hh.ids < 1 | hh.ids != round(hh.ids)))) {
    stop("the birth placement needs hh.ids to be positive whole numbers (GlobalMix India #65)",
         call. = FALSE)
  }
  hh <- hh_layer(dat)
  if (!is.na(hh)) {
    fun <- get_nwparam(dat, network = hh)$arrivals.FUN
    if (is.null(fun) || !identical(deparse(body(fun)), deparse(body(hh_join_birth)))) {
      stop("the household layer must use arrivals.FUN = hh_join_birth, which keeps the ",
           "mother's household that place_births() sets for each newborn (GlobalMix India #65)",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

# Composition of every household among the living (birth_hh_terms), one row per
# household id
birth_hh_features <- function(age, sex, hh, live) {
  H <- max(c(0L, hh[live]), na.rm = TRUE)
  k <- ifelse(age < 10, 1L, ifelse(age < 20, 2L, ifelse(age < 40, ifelse(sex == "F", 3L, 4L),
                                                     ifelse(age < 60, 5L, 6L))))
  sel <- live & !is.na(hh)
  x <- vapply(1:6, function(j) tabulate(hh[sel & k == j], nbins = H), numeric(H))
  x <- matrix(x, nrow = H)
  x <- cbind(x, x[, 1] == 0)
  colnames(x) <- birth_hh_terms
  x
}

# Weight of a woman of each age as a mother, from her five-year age group; 0
# outside 15-49
birth_age_wt <- function(age, w) {
  g <- findInterval(age, birth_age_breaks)
  out <- numeric(length(age))
  ok <- g >= 1 & g < length(birth_age_breaks)
  out[ok] <- w[g[ok]]
  out
}

#' @title Birth Placement
#'
#' @description Draws a mother for each newborn, places the newborn in her
#'   household and records the birth (GlobalMix India #65).
#'
#' @param dat Main `netsim_dat` object.
#' @param at Current time step.
#' @param new_ids Positional ids of this step's newborns.
#'
#' @details
#' Mothers are drawn without replacement among the living women aged 15-49
#' whose last birth (`last.birth`) is at least `birth.interval.min` days ago,
#' with weight `birth.age.weight` of their five-year age group times
#' `exp(sum(birth.hh.coef * x))` on their household's composition `x` (members
#' aged 0-9, 10-19, women and men 20-39, 40-59 and 60 and over, and whether no
#' member is under 10). Each newborn takes its mother's `hh.ids`, and the
#' household layer's arrival rule, [hh_join_birth()], keeps it.
#'
#' The step's births by the mother's five-year age group are the epi series
#' `nNew.mage1` to `nNew.mage7` (15-19 to 45-49).
#'
#' @return `dat` with the newborns' `hh.ids` and the mothers' `last.birth` set.
#' @keywords internal
place_births <- function(dat, at, new_ids) {
  active <- get_attr(dat, "active")
  age <- get_attr(dat, "age")
  sex <- get_attr(dat, "sex")
  hh <- get_attr(dat, "hh.ids")
  last.birth <- get_attr(dat, "last.birth")
  gap <- get_param(dat, "birth.interval.min")

  live <- active == 1
  live[new_ids] <- FALSE
  ok <- live & sex == "F" & age >= 15 & age < 50 & !is.na(hh) &
    (is.na(last.birth) | at - last.birth >= gap)
  ids <- which(ok)
  x <- birth_hh_features(age, sex, hh, live)[hh[ids], , drop = FALSE]
  w <- birth_age_wt(age[ids], get_param(dat, "birth.age.weight")) *
    exp(drop(x %*% get_param(dat, "birth.hh.coef")))
  ids <- ids[w > 0]; w <- w[w > 0]
  if (length(ids) < length(new_ids)) {
    stop(length(new_ids), " births at step ", at, " but only ", length(ids), " women can be ",
         "mothers (aged 15-49, in a household, no birth in the last birth.interval.min days)",
         call. = FALSE)
  }
  mothers <- ids[sample.int(length(ids), length(new_ids), prob = w)]
  hh[new_ids] <- hh[mothers]
  last.birth[mothers] <- at
  dat <- set_attr(dat, "hh.ids", hh)
  dat <- set_attr(dat, "last.birth", last.birth)
  by_age <- tabulate(findInterval(age[mothers], birth_age_breaks), nbins = length(birth_age_breaks) - 1)
  for (g in seq_along(by_age)) {
    dat <- set_epi(dat, paste0("nNew.mage", g), by_age[g], at = at)
  }
  dat
}

# Mothers of the starting population's youngest children (GlobalMix India #65).
# Each child younger than birth.interval.min days is attributed to a woman of
# their household who was 15-49 at the birth, drawn by birth.age.weight at her
# age then, youngest child first and at most one child per woman, so that she
# is not a mother again within the interval. A second such child in a household
# (a cousin, or a twin) goes to another woman who was 15-49 at the birth if
# there is one, and otherwise has no mother on record. last.birth is the step
# of the birth: 1 less the child's age in days.
init_last_birth <- function(dat) {
  active <- get_attr(dat, "active")
  age <- get_attr(dat, "age")
  sex <- get_attr(dat, "sex")
  hh <- get_attr(dat, "hh.ids")
  gap <- get_param(dat, "birth.interval.min")
  w <- get_param(dat, "birth.age.weight")
  last.birth <- rep(NA_real_, length(active))
  kids <- which(active == 1 & age * days_per_year < gap & !is.na(hh))
  kids <- kids[order(age[kids])]
  women <- which(active == 1 & sex == "F" & !is.na(hh))
  by_hh <- split(women, hh[women])
  for (k in kids) {
    cand <- by_hh[[as.character(hh[k])]]
    if (is.null(cand)) next
    p <- birth_age_wt(age[cand] - age[k], w) * is.na(last.birth[cand])
    if (!any(p > 0)) next
    m <- cand[sample.int(length(cand), 1, prob = p)]
    last.birth[m] <- 1 - round(age[k] * days_per_year)
  }
  set_attr(dat, "last.birth", last.birth)
}

#' @param new_ids Positional ids of the nodes arriving at this time step.
#' @param network Index of the household clique layer.
#'
#' @details
#' \code{hh_join_birth} is the \code{arrivals.FUN} of the household layer, built
#' with \code{netclique(nw, "hh.ids", arrivals = "join", arrivals.FUN =
#' hh_join_birth)}. A newborn whose household the arrivals module set (its
#' mother's, on the relationship design; GlobalMix India #65) keeps it. Any
#' other newborn joins the household of a uniformly drawn existing 0-9y member,
#' so households are drawn in proportion to their number of 0-9y members (the
#' rule before engine 1.10.0, which the legacy four-layer model keeps). EpiModel
#' then connects the newborn to every active member of that household,
#' including other newborns placed there in the same step.
#'
#' @rdname moduleset-gmc19
#' @export
hh_join_birth <- function(dat, at, new_ids, network) {
  active <- get_attr(dat, "active")
  age.grp <- get_attr(dat, "age.grp")
  hh.ids <- get_attr(dat, "hh.ids")

  placed <- hh.ids[new_ids]
  todo <- is.na(placed)
  if (!any(todo)) {
    return(placed)
  }

  # Existing households that contain a young (0-9y) member; new arrivals (age 0)
  # are placed into one of these. The pool leaves out this step's newborns
  # (their hh.ids is not set yet) and this step's dead, who are still in the
  # attribute list until nwupdate.FUN removes them after this call. Guard the
  # degenerate case (issue #39): if no 0-9y node exists the pool is empty, so
  # return NA and leave the new nodes in no household (isolates on the layer)
  # rather than misassigning them. Index via sample.int so a length-1 pool is
  # drawn from correctly (base sample() would treat a single numeric value as
  # 1:value).
  existing <- active == 1
  existing[new_ids] <- FALSE
  hh_pool <- hh.ids[existing & age.grp == "0-9y" & !is.na(hh.ids)]
  if (length(hh_pool) == 0) {
    return(placed)
  }
  placed[todo] <- hh_pool[sample.int(length(hh_pool), sum(todo), replace = TRUE)]
  placed
}
