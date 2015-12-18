\name{SM.MixReparametrized}
\alias{SM.MixReparametrized}
\title{
summary of the output produced by K.MixReparametrized 
}
\description{
This is a generic function that summarizes the MCMC samples produced by K.MixReparametrized. The function invokes
several estimation methods which choice depends on the unimodality or multimodality of the argument. 
}
\usage{
SM.MixReparametrized(xobs, estimate)
}
\arguments{
  \item{xobs}{
vector of observations 
}
  \item{estimate}{
output of K.MixReparametrized 
}
}
\details{
This function outputs posterior  point estimates for all parameters of the mixture model. They mostly
differ from the generaly useless posterior means. The output summarizes unimodal MCMC samples by computing
measures of centrality, including mean and median, while multimodal outputs require a pre-processing, due to the
label switching phenomenon (Jasra et al., 2005). The summary measures are then computed after performing a multi-dimensional k-means clustering (Hartigan and Wong, 1979) following the suggestion of
Fruhwirth-Schnatter (2006). 
}
\value{
  \item{Mean }{vector of mean and median of simulated draws from the conditional posterior of the mixture model mean}
  \item{Sd }{vector of mean and median of simulated draws from the conditional posterior of the mixture model standard deviation}
  \item{Phi }{vector of mean and median of simulated draws from the conditional posterior of the radial coordinate}
  \item{Angles. 1. }{vector of means of the angular coordinates used for the component means in the mixture distribution}
  \item{Angles. 2. }{vector of means of the angular coordinates used for the component standard deviations in the mixture distribution}
  \item{weight.i }{vector of mean and median of simulated draws from the conditional posterior of the component weights of the mixture distribution; \eqn{i=1, \ldots, k}}
  \item{mean.i }{vector of mean and median of simulated draws from the conditional posterior of the component means of the mixture distribution; \eqn{i=1, \ldots, k}}
  \item{sd.i }{vector of mean and median of simulated draws from the conditional posterior of the component standard deviations of the mixture distribution; \eqn{i=1, \ldots, k}}
  \item{Acc rat }{vector of final acceptance rate of the proposal distributions of the algorithm with no calibration
stage for the proposal scales}
  \item{Opt scale }{vector of optimal proposal scales obtained the by calibration stage}
}
\references{
Jasra, A., Holmes, C. and Stephens, D. (2005). Markov Chain Monte Carlo methods and the label switching problem in
Bayesian mixture modeling. Statistical Science, 20, 50--67.

Hartigan, J. A. and Wong, M. A. (1979). A K-means clustering algorithm. Applied Statistics 28, 100--108.

Fruhwirth-Schnatter, S. (2006). Finite mixture and Markov switching models. Springer-Verlag.
}
\author{
Kaniav Kamary
}
\note{
For multimodal outputs such as the mixture model weights, component means, and component variances, for each MCMC draw,
first the labels of the weights \eqn{p_i, i=1, \ldots, k} and corresponding component means and standard deviations are permuted in
such a way that \eqn{p_1\le \ldots \le p_k}. Then the component means and standard deviations are jointly partitioned into \eqn{k}
clusters by applying a standard k-means algorithm with \eqn{k} clusters to a sample of  size \eqn{Tk} (where \eqn{T} is the number of
iterations), following Fruhwirth-Schnatter (2006) method. For each group,  cluster centers are considered as parameter
estimates. 
}
\seealso{
\code{\link[Ultimixt]{K.MixReparametrized}}
}
\examples{
data(faithful)
xobs=faithful[50:100,1]
#estimate=K.MixReparametrized(xobs, k=2, alpha0=.5, alpha=.5, Nsim=1e4)
#summari=SM.MixReparametrized(xobs,estimate)
}
\keyword{summary statistics}
\keyword{mixture parameters}
\keyword{k-means clustering parameters}