\name{K.MixReparametrized}
\alias{K.MixReparametrized}
\title{
Sample from a Gaussian mixture posterior associated with a noninformative prior and obtained by Metropolis-within-Gibbs sampling
}
\description{
This function returns a sample simulated from the posterior distribution of the parameters of a Gaussian mixture under a non-informative prior. This prior is derived from a mean-variance reparameterisation of the mixture distribution, as proposed by Kamary et al. (2015). The algorithm is a Metropolis-within-Gibbs scheme with an adaptive calibration of the proposal distribution scales. Adaptation is driven by the formally optimal acceptance rates of \eqn{0.44} and \eqn{0.234} in one and larger dimensions, respectively (Roberts et al.,1997). This algorithm monitors the convergence of the MCMC sequences via Gelman's and Rubin's (1992) criterion.
}
\usage{
K.MixReparametrized(xobs, k, alpha0, alpha, Nsim)
}
\arguments{
  \item{xobs}{
vector of the observations or dataset
}
  \item{k}{
number of components in the mixture model
}
\item{alpha0}{
hyperparameter of Dirichlet prior distribution of the mixture model weights which is .5 by default 
}
\item{alpha}{
hyperparameter of beta prior distribution of the radial coordinate which is .5 by default
}
\item{Nsim}{
number of MCMC iterations after calibration step of proposal scales
}
}
\details{
The output of this function contains a simulated sample for each parameter of the mixture distribution, the evolution of
the proposal scales and acceptance rates over the number of iterations during the calibration stage, and their final
values after calibration.
}
\value{
The output of this function is a list of the following variables, where the dimension of the vectors is the number of simulations:
  \item{mean global  }{vector of simulated draws from the conditional posterior of the mixture model mean}
  \item{sigma global }{vector of simulated draws from the conditional posterior of the mixture model standard deviation}
  \item{weights }{matrix of simulated draws from the conditional posterior of the mixture model weights with a number of columns equal to the number of components \eqn{k}}
  \item{angles xi }{matrix of simulated draws from the conditional posterior of the angular coordinates of the component standard deviations with a number of columns equal to \eqn{k-1}}
  \item{phi }{vector of simulated draws from the conditional posterior of the radian coordinate}
  \item{angles varpi }{matrix of simulated draws from the conditional posterior of the angular coordinates of the component means with a number of columns equal to \eqn{k-2}}
  \item{accept rat  }{vector of resulting acceptance rates of the proposal distributions without calibration step of the proposal scales}
  \item{optimal para }{vector of resulting proposal scales after optimisation obtained by adaptive MCMC}
  \item{adapt rat }{ list of acceptance rates of batch of 50 iterations obtained when calibrating the proposal scales by adaptive MCMC. The number of columns depends on the number of proposal distributions.}
  \item{adapt scale }{list of proposal scales calibrated by adaptive MCMC for each batch of 50 iterations with respect to the optimal acceptance rate. The number of columns depends on the number of proposal distribution scales.}
  \item{component means }{matrix of MCMC samples of the component means of the mixture model with a number of columns equal to \eqn{k}}
  \item{component sigmas }{matrix of MCMC samples of the component standard deviations of the mixture model with a number of columns equal to \eqn{k}}
}
\references{
Kamary, K., Lee, J.Y., and Robert, C.P. (2015) Non-informative reparameterisation of location-scale mixtures. arXiv.

Robert, C. and Casella, G. (2009). Introducing Monte Carlo Methods with R. Springer-Verlag.

Roberts, G. O., Gelman, A. and Gilks, W. R. (1997). Weak convergence and optimal scaling of random walk Metropolis algorithms. Ann. Applied Probability, 7, 110--120.

Gelman, A. and Rubin, D. (1992). Inference from iterative simulation using multiple sequences (with discussion). Statistical Science, 457--472.
}
\author{
Kaniav Kamary
}
\note{
If the number of MCMC iterations specified in the input of this function exceeds 15,000, after each 1000 supplementry iterations the convergence of simulated chains is checked using the convergence monitoring technique by Gelman and Rubin (1992).
}

\seealso{
\code{\link[Ultimixt]{Ultimixt}}
}
\examples{
data(faithful)
xobs=faithful[50:100,1]
#estimate=K.MixReparametrized(xobs, k=2, alpha0=.5, alpha=.5, Nsim=1e4)
}

\keyword{mixture distribution}
\keyword{non informative parametrisation}
