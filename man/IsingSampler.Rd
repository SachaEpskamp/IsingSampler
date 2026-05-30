\name{IsingSampler}
\alias{IsingSampler}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample states from the Ising model
}
\description{
This function samples states from the Ising model using one of three methods. See details.
}
\usage{
IsingSampler(n, graph, thresholds, beta = 1, nIter = 100, responses = c(0L, 1L), 
    method = c("MH", "CFTP", "direct"), CFTPretry = 10, constrain)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Number of states to draw
}
  \item{graph}{
Square matrix indicating the weights of the network. Must be symmetrical with 0 as diagonal.
}
  \item{thresholds}{
Vector indicating the thresholds, also known as the external field.
}
  \item{beta}{
Scalar indicating the inverse temperature.
}
  \item{nIter}{
Number of iterations in the Metropolis and exact sampling methods.
}
  \item{responses}{
Response options. Typically set to \code{c(-1L, 1L)} or \code{c(0L, 1L)} (default). Two or more response options are supported, and they may be any numeric values, including non-integers (e.g. \code{c(-1L, 0L, 1L)}, \code{1:5} or \code{seq(-1, 1, by = 0.5)}). With more than two options the conditional distribution of each node is the categorical (softmax) generalization of the binary Ising conditional. The pairwise interaction always enters the Hamiltonian as \code{J_ij * s_i * s_j}, i.e. through the (numeric) product of the response values, and each node has a single threshold. When all response options are integer-valued the result is an integer matrix (as before); otherwise it is a numeric (double) matrix.
}
  \item{method}{
The sampling method to use. Must be \code{"MH"}, \code{"CFTP"} or \code{"direct"}. See details.
}
 \item{CFTPretry}{
 The amount of times a sample from CFTP may be retried. If after 100 couplings from the past the chain still results in \code{NA} values the chain is reset with different random numbers. Be aware that data that requies a lot of CFTP resets might not resemble exact samples anymore.
 }
 \item{constrain}{A (number of samples) by (number of nodes) matrix with samples that need be constrained; \code{NA} indicates that the sample is unconstrained. Defaults to a matrix of \code{NA}s.}
}
\details{
This function uses one of three sampling methods. \code{"MH"} can be used to sample using a Metropolis-Hastings algorithm. The chain is initiated with random values from the response options, then for each iteration each node is redrawn from its full-conditional distribution given all other nodes and parameters. For two response options this is the usual binary update; for more than two response options the node is drawn from the categorical (softmax) full-conditional over all response options. Typically, 100 of such iterations should suffice for the chain to converge, though strongly coupled models with more than two response options may require a larger \code{nIter}.

The second method, \code{"CFTP"} enhances the Metropolis-Hastings algorithm with Coupling from the Past (CFTP; Murray, 2007) to draw exact samples from the distribution. This is slower than the default Metropolis-Hastings but guarantees exact samples. However, it does depend on the graph structure and the number of nodes if these exact samples can be obtained in feasable time. \code{"CFTP"} is currently only implemented for two response options; for more than two response options use \code{"MH"} or \code{"direct"}.

The third option, \code{"direct"}, simply computes for every possibly state the probability and draws samples directly from the distribution of states by using these probabilities. This also guarantees exact samples, but quickly becomes intractible (roughly above 10 nodes).
}
\value{
A matrix containing samples of states.
}
\references{
Murray, I. (2007). Advances in Markov chain Monte Carlo methods.}
\author{
Sacha Epskamp (mail@sachaepskamp.com)
}

\seealso{
\code{\link{IsingSampler-package}} for examples
}
\examples{
## See IsingSampler-package help page
}
