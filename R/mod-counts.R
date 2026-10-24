# Instrumentation for the redo estimands (GlobalMix India #77): persons ever
# infected including the dead and the day-1 cohort, years of life lost from
# exact age at death, eta-squared of degree on age band, and a baseline degree
# quartile with random tie-breaking that leaves the run's random stream alone.
# None of these draws from the simulation's random numbers, so a run's seeded
# output is otherwise identical with or without them.

# Departed persons --------------------------------------------------------------------
#
# ever.inf.num and ever.inf.age* count the living, so the people who died after
# an infection leave them while every episode they had stays in se.flow. The
# departures module records, cumulatively over the run, the departed who had
# ever been infected or symptomatic, by age band at death and for infants, and
# the members of the day-1 cohort (entrTime 1) among them; the prevalence
# module adds them to the living, so the person counts never fall when someone
# dies, and reports the cohort on its own.

departed_counts <- function(dat) {
  d <- dat$run$departed
  if (is.null(d)) {
    nAge <- length(get_param(dat, "age.breaks")) - 1
    d <- list(n = 0, inf = 0, inf.age = numeric(nAge), inf.infant = 0, sympt = 0,
              sympt.age = numeric(nAge), cohort = 0, inf.cohort = 0)
  }
  d
}

record_departed <- function(dat, ids) {
  d <- departed_counts(dat)
  if (length(ids) > 0) {
    nAge <- length(d$inf.age)
    inf <- get_attr(dat, "ever.inf")[ids] %in% 1
    sympt <- get_attr(dat, "ever.sympt")[ids] %in% 1
    ag <- cut(get_attr(dat, "age")[ids], get_param(dat, "age.breaks"), labels = FALSE, right = FALSE)
    is_infant <- get_attr(dat, "is_infant", override.null.error = TRUE)
    cohort <- get_attr(dat, "entrTime")[ids] %in% 1
    d$n <- d$n + length(ids)
    d$inf <- d$inf + sum(inf)
    d$inf.age <- d$inf.age + tabulate(ag[inf], nAge)
    if (!is.null(is_infant)) d$inf.infant <- d$inf.infant + sum(inf & is_infant[ids] %in% TRUE)
    d$sympt <- d$sympt + sum(sympt)
    d$sympt.age <- d$sympt.age + tabulate(ag[sympt], nAge)
    d$cohort <- d$cohort + sum(cohort)
    d$inf.cohort <- d$inf.cohort + sum(cohort & inf)
  }
  dat$run$departed <- d
  dat
}

# Years of life lost -----------------------------------------------------------------
#
# yll.le is the remaining life expectancy at exact ages 0, 1, 2, ... years of
# the setting's life table, supplied at deploy like mort.rates. A death at
# exact age a loses the expectancy at a, interpolated linearly between the
# integer ages and held at the last value beyond them (the constant hazard the
# table's open age interval assumes). yll.discount is a continuous discount rate
# r on those years, (1 - exp(-r L)) / r; 0 is no discounting.

yll_remaining <- function(age, le) {
  n <- length(le)
  i <- pmin(pmax(floor(age), 0), n - 1)
  f <- pmin(pmax(age - i, 0), 1)
  i1 <- i + 1
  i2 <- pmin(i + 2, n)
  le[i1] + f * (le[i2] - le[i1])
}

yll_discounted <- function(L, r) {
  if (r == 0) L else (1 - exp(-r * L)) / r
}

check_yll_params <- function(dat) {
  le <- get_param(dat, "yll.le", override.null.error = TRUE)
  if (!is.null(le) && !(is.numeric(le) && length(le) >= 2 && all(is.finite(le)) && all(le >= 0))) {
    stop("yll.le must be the remaining life expectancy at ages 0, 1, 2, ... years, finite and ",
         "non-negative (GlobalMix India #77)", call. = FALSE)
  }
  r <- param_or(dat, "yll.discount", 0.03)
  if (!(is.numeric(r) && length(r) == 1 && is.finite(r) && r >= 0 && r < 1)) {
    stop("yll.discount must be one continuous discount rate in [0, 1) (GlobalMix India #77)", call. = FALSE)
  }
  invisible(TRUE)
}

# Degree on age -----------------------------------------------------------------------

# Share of the variance of x that lies between the groups g (eta-squared): one
# less the within-group sum of squares over the total. NA with fewer than two
# values or no variance.
eta_squared <- function(x, g) {
  ok <- is.finite(x) & !is.na(g)
  x <- x[ok]
  g <- g[ok]
  if (length(x) < 2) return(NA_real_)
  sst <- sum((x - mean(x))^2)
  if (sst == 0) return(NA_real_)
  m <- rowsum(x, g) / tabulate(match(g, sort(unique(g))))
  ssw <- sum((x - m[as.character(g), 1])^2)
  1 - ssw / sst
}

# Quartile of a score among the living, 1 to 4 by rank, NA for the rest; ties
# are broken as `ties` says (rank()'s ties.method)
score_quartile <- function(score, live, ties = "random") {
  q <- rep(NA_integer_, length(score))
  n <- sum(live)
  if (n == 0) return(q)
  r <- rank(score[live], ties.method = ties)
  q[live] <- as.integer(pmin(pmax(ceiling(r / n * 4), 1), 4))
  q
}

# Evaluates expr with its own random stream and puts the run's stream back
# exactly where it was, so a random draw inside an instrument changes nothing
# the simulation draws afterwards. The side stream is seeded from the run's
# state, so simulations with different seeds break ties differently.
with_side_rng <- function(expr) {
  genv <- globalenv()
  had <- exists(".Random.seed", envir = genv, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = genv, inherits = FALSE)
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = genv)
    } else if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
      rm(".Random.seed", envir = genv)
    }
  }, add = TRUE)
  seed <- if (had) sum(as.numeric(old[-1]) %% 1e6) %% 2147483646 + 1 else 77
  set.seed(as.integer(seed))
  force(expr)
}
