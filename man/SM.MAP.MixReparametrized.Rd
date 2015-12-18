\name{SM.MAP.MixReparametrized}
\alias{SM.MAP.MixReparametrized}
\title{
summary of the output produced by K.MixReparametrized
}
\description{
Label switching in a simulated Markov chain produced by K.MixReparametrized is removed by the technique of Marin et al.
(2004). Namely, component labels are reorded by the shortest Euclidian distance between a posterior sample and the maximum a
posteriori (MAP) estimate. Let \eqn{\theta_i} be the \eqn{i}-th vector of computed component means, standard deviations
and weights. The MAP estimate is derived from the MCMC sequence and denoted by \eqn{\theta_{MAP}}. For a permutation
\eqn{\tau \in \Im_k} the labelling of \eqn{\theta_i} is reordered by
\deqn{
	\tilde{\theta}_i=\tau_i(\theta_i)
	}
where \eqn{\tau_i=\arg \min_{\tau \in \Im_k} \mid \mid \tau(\theta_i)-\theta_{MAP}\mid \mid}.

Angular parameters \eqn{\xi_1^{(i)}, \ldots, \xi_{k-1}^{(i)}} and \eqn{\varpi_1^{(i)}, \ldots, \varpi_{k-2}^{(i)}}s are
derived from \eqn{\tilde{\theta}_i}. There exists an unique solution in \eqn{\varpi_1^{(i)}, \ldots, \varpi_{k-2}^{(i)}}
while there are multiple solutions in \eqn{\xi^{(i)}} due to the symmetry of \eqn{\mid\cos(\xi) \mid} and
\eqn{\mid\sin(\xi) \mid}. The output of \eqn{\xi_1^{(i)}, \ldots, \xi_{k-1}^{(i)}} only includes angles on \eqn{[-\pi, \pi]}.

The label of components of \eqn{\theta_i} (before the above transform) is defined by

\deqn{
	\tau_i^*=\arg \min_{\tau \in \Im_k}\mid \mid \theta_i-\tau(\theta_{MAP}) \mid \mid.
}
The number of label switching occurrences is defined by the number of changes in \eqn{\tau^*}.
}
\usage{
SM.MAP.MixReparametrized(estimate, xobs, alpha0, alpha)
}
\arguments{
  \item{estimate}{
Output of K.MixReparametrized
}
  \item{xobs}{
Data set
}
  \item{alpha0}{
Hyperparameter of Dirichlet prior distribution of the mixture model weights
}
  \item{alpha}{
Hyperparameter of beta prior distribution of the radial coordinate
}
}
\details{
Details.}
\value{
  \item{MU }{Matrix of MCMC samples of the component means of the mixture model}
  \item{SIGMA }{Matrix of MCMC samples of the component standard deviations of the
mixture model}
  \item{P }{Matrix of MCMC samples of the component weights of the mixture model}
  \item{Ang_SIGMA }{Matrix of computed \eqn{\xi}'s corresponding to SIGMA}
  \item{Ang_MU }{Matrix of computed \eqn{\varpi}'s corresponding to MU. This output only appears when \eqn{k > 2}.}
  \item{Global_Mean }{Mean, median and \eqn{95\%} credible interval for the global mean parameter}
  \item{Global_Std }{Mean, median and \eqn{95\%} credible interval for the global standard deviation parameter}
  \item{Phi }{Mean, median and \eqn{95\%} credible interval for the radius parameter}
  \item{component_mu }{Mean, median and \eqn{95\%} credible interval of MU}
  \item{component_sigma }{Mean, median and \eqn{95\%} credible interval of SIGMA}
  \item{component_p }{Mean, median and \eqn{95\%} credible interval of P}
  \item{l_stay }{Number of MCMC iterations between changes in labelling}
  \item{n_switch }{Number of label switching occurrences}
}
\references{
Marin, J.-M., Mengersen, K. and Robert, C. P. (2004) Bayesian Modelling
and Inference on Mixtures of Distributions, Handbook of Statistics, Elsevier,
Volume 25, Pages 459--507.
}
\author{
Kate Lee
}
\note{
Note.}
\seealso{
\code{\link[Ultimixt]{K.MixReparametrized}}
}
\examples{
data(faithful)
xobs=faithful[50:100,1]
#estimate=K.MixReparametrized(xobs,k=2,alpha0=0.5,alpha=0.5,Nsim=1e4)
#result=SM.MAP.MixReparametrized(estimate,xobs,alpha0=0.5,alpha=0.5)
}
\keyword{summary statistics}
\keyword{maximum a posteriori probability}