# Severity per infection on the twelve severity bands (GlobalMix India #71), for
# the test fixtures that run init_gmc19()
sev_test_breaks <- c(0, 0.5, 1, 5, 10, 20, 30, 40, 50, 60, 70, 80, Inf)
sev_test_params <- function(prop.clinical = rep(0.5, 12), ihr = rep(0.02, 12), ifr = rep(0.002, 12)) {
  list(sev.age.breaks = sev_test_breaks, prop.clinical = prop.clinical, ihr = ihr, ifr = ifr)
}
