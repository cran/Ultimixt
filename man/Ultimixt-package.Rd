\name{Ultimixt-package}
\alias{Ultimixt-package}
\alias{Ultimixt}
\docType{package}
\title{
set of R functions for estimating the parameters of mixture distribution with a Bayesian non-informative
prior
}
\description{
Despite a comprehensive literature on estimating mixtures of Gaussian distributions, there does not exist a
well-accepted reference Bayesian approach to such models. One reason for the difficulty is the general prohibition
against using improper priors (Fruhwirth-Schnatter, 2006) due to the ill-posed nature of such statistical objects.
Kamary, Lee and Robert (2017) took advantage of a mean-variance reparametrisation of a Gaussian mixture model to propose
improper but valid reference priors in this setting. This R package implements the proposal and computes posterior
estimates of the parameters of a Gaussian mixture distribution. The approach applies with an arbitrary number of
components. The Ultimixt R package contains an MCMC algorithm function and further functions for
summarizing and plotting posterior estimates of the model parameters for any number of components. 
}
\details{
\tabular{ll}{
Package: \tab Ultimixt\cr
Type: \tab Package\cr
Version: \tab 2.1\cr
Date: \tab 2017-03-07\cr
License: \tab GPL (>=2.0)\cr
}
Beyond simulating MCMC samples from the posterior distribution of the Gaussian mixture model, this package also produces
summaries of the  MCMC outputs through numerical and graphical methods. 

Note:  The proposed parameterisation of the Gaussian mixture distribution is given by 
\deqn{
f(x| \mu, \sigma , {\bf p}, \varphi, {\bf \varpi, \xi})=\sum_{i=1}^k p_i f\left(x| \mu + \sigma \gamma_i/\sqrt{p_i}, \sigma \eta_i/\sqrt{p_i}\right)
}
under the non-informative prior \eqn{\pi(\mu, \sigma)=1/\sigma}. Here, the vector of the \eqn{\gamma_i=\varphi
\Psi_i\Big({\bf \varpi}, {\bf p}\Big)_i}'s belongs to an hypersphere of radius \eqn{\varphi} intersecting with an
hyperplane.  It is thus expressed in terms of spherical coordinates within that hyperplane that depend on \eqn{k-2}
angular coordinates \eqn{\varpi_i}. Similarly, the vector of \eqn{\eta_i=\sqrt{1-\varphi^2}\Psi_i\Big({\bf
\xi}\Big)_i}'s can be turned
into a spherical coordinate in a k-dimensional Euclidean space, involving a radial coordinate
\eqn{\sqrt{1-\varphi^2}} and \eqn{k-1} angular coordinates \eqn{\xi_i}. A natural prior for \eqn{\varpi} is made of uniforms, \eqn{\varpi_1, \ldots, \varpi_{k-3}\sim U[0, \pi]} and \eqn{\varpi_{k-2} \sim U[0, 2\pi]}, and for \eqn{\varphi}, we consider a beta prior \eqn{Beta(\alpha, \alpha)}. A reference prior on the angles \eqn{\xi} is \eqn{(\xi_1, \ldots, \xi_{k-1})\sim U[0, \pi/2]^{k-1}} and a Dirichlet prior \eqn{Dir(\alpha_0, \ldots, \alpha_0)} is assigned to the weights \eqn{p_1, \ldots, p_k}. 

For a Poisson mixture, we consider 
\deqn{
f(x|\lambda_1, \ldots, \lambda_k)=\frac{1}{x!}\sum_{i=1}^k p_i \lambda_i^x e^{-\lambda_i}
}
with a reparameterisation as \eqn{\lambda=\bf{E}[X]} and \eqn{\lambda_i=\lambda
\gamma_i/p_i}. In this case, we can use the equivalent to the Jeffreys prior for the Poisson
distribution, namely, \eqn{\pi(\lambda)=1/\lambda}, since it leads to a
well-defined posterior with a single positive observation. 
}
\author{
Kaniav Kamary

Maintainer: \email{kamary@ceremade.dauphine.fr}
}
\references{
Fruhwirth-Schnatter, S. (2006). Finite Mixture and Markov Switching Models. Springer-Verlag, New York, New York. 

Kamary, K., Lee, J.Y., and Robert, C.P. (2017) Weakly informative reparameterisation for location-scale mixtures. arXiv.

}
\keyword{ package }
\seealso{
\code{\link[Ultimixt]{Ultimixt}}
}
\examples{
	#K.MixReparametrized(faithful[,2], k=2, alpha0=.5, alpha=.5, Nsim=10000)
}