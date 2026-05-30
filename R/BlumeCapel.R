# Wrapper around IsingSampler for the Blume-Capel model.
#
# The Blume-Capel model is an Ising model with an additional on-site quadratic
# ("single-ion anisotropy" / crystal-field) term in the Hamiltonian:
#
#   H = - sum_i tau_i s_i - sum_{i<j} omega_ij s_i s_j + sum_i delta_i s_i^2
#
# Compared with IsingSampler() the only new ingredient is 'delta'. It is a true
# extension only with more than two (ordered) response options; the default
# response options are the spin-1 values c(-1, 0, 1). With delta_i > 0 the middle
# category is favoured; with delta_i < 0 the extreme categories are favoured.
#
# This is a thin wrapper: it sets the default responses and a required 'delta'
# and forwards everything to IsingSampler(). Exact sampling ("CFTP") is binary
# only and therefore not offered here (delta has no identifiable effect on binary
# responses anyway); use method = "MH" (default) or "direct".
BlumeCapelSampler <- function(n, graph, thresholds, delta, beta = 1, nIter = 100,
                              responses = c(-1L, 0L, 1L),
                              method = c("MH", "direct"), CFTPretry = 10, constrain)
{
  if (missing(delta))
  {
    stop("'delta' must be supplied for the Blume-Capel model (use IsingSampler() for the delta = 0 / Ising case).")
  }
  method <- match.arg(method)

  args <- list(n = n, graph = graph, thresholds = thresholds, beta = beta,
               nIter = nIter, responses = responses, method = method,
               CFTPretry = CFTPretry, delta = delta)
  if (!missing(constrain)) args$constrain <- constrain

  do.call(IsingSampler, args)
}
