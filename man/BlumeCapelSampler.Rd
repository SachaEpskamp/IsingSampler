\name{BlumeCapelSampler}
\alias{BlumeCapelSampler}
\title{
Sample states from the Blume-Capel model
}
\description{
Samples states from the Blume-Capel model, an Ising model with an additional on-site quadratic ("single-ion anisotropy" or "crystal-field") term. This is a thin wrapper around \code{\link{IsingSampler}} that sets ordered (spin-1) response options by default and requires the quadratic parameter \code{delta}.
}
\usage{
BlumeCapelSampler(n, graph, thresholds, delta, beta = 1, nIter = 100,
    responses = c(-1L, 0L, 1L), method = c("MH", "direct"),
    CFTPretry = 10, constrain)
}
\arguments{
  \item{n}{
Number of states to draw.
}
  \item{graph}{
Square matrix indicating the weights of the network. Must be symmetrical with 0 as diagonal.
}
  \item{thresholds}{
Vector indicating the thresholds (external field). A single value is recycled over all nodes.
}
  \item{delta}{
The quadratic (Blume-Capel) parameter, entering the Hamiltonian as \code{delta_i * s_i^2}. Either a single value (recycled over all nodes) or one value per node. Must be supplied. A positive \code{delta} favours the middle category; a negative \code{delta} favours the extreme categories.
}
  \item{beta}{
Scalar indicating the inverse temperature.
}
  \item{nIter}{
Number of iterations in the Metropolis sampler.
}
  \item{responses}{
Ordered response options, defaulting to the spin-1 values \code{c(-1L, 0L, 1L)}. Other (numeric) response options may be supplied, but the \code{delta} term only has an identifiable effect with more than two options.
}
  \item{method}{
The sampling method, either \code{"MH"} (Metropolis-Hastings, default) or \code{"direct"}. Exact sampling (\code{"CFTP"}) is implemented for two response options only and is therefore not available here; \code{delta} has no identifiable effect on binary responses in any case.
}
  \item{CFTPretry}{
Passed to \code{\link{IsingSampler}}; unused for the methods offered here.
}
  \item{constrain}{
A (number of samples) by (number of nodes) matrix with samples that need be constrained; \code{NA} indicates that the sample is unconstrained.
}
}
\details{
The Blume-Capel Hamiltonian implemented here is
\deqn{H(s) = -\sum_i \tau_i s_i - \sum_{i<j} \omega_{ij} s_i s_j + \sum_i \delta_i s_i^2,}
with thresholds \eqn{\tau}, network weights \eqn{\omega} (\code{graph}) and quadratic terms \eqn{\delta} (\code{delta}); states are drawn with probability proportional to \eqn{\exp(-\beta H(s))}. Setting \code{delta = 0} recovers the ordinary (multi-level) Ising model, so \code{BlumeCapelSampler(..., delta = 0)} is equivalent to \code{\link{IsingSampler}}.
}
\value{
A matrix containing samples of states.
}
\author{
Sacha Epskamp (mail@sachaepskamp.com)
}
\seealso{
\code{\link{IsingSampler}}, \code{\link{IsingLikelihood}}
}
\examples{
\dontrun{
# A small ferromagnetic network:
P <- 5
W <- matrix(0.5, P, P); diag(W) <- 0

# Low crystal field: ordered (mostly +1 / -1):
X0 <- BlumeCapelSampler(1000, W, thresholds = 0, delta = 0, beta = 1)

# High crystal field: disordered "0" phase (mostly the middle category):
X3 <- BlumeCapelSampler(1000, W, thresholds = 0, delta = 3, beta = 1)

mean(X0 == 0)  # low
mean(X3 == 0)  # high
}
}
