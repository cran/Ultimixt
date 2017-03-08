\name{SM.MixPois}
\alias{SM.MixPois}
\title{
summary of the output produced by K.MixPois
}
\description{
This generic function summarizes the MCMC samples produced by K.MixPois when several estimation methods have been invoked depending on the unimodality or multimodality of the argument.  
}
\usage{
SM.MixPois(estimate, xobs)
}
\arguments{
  \item{estimate}{
output of K.MixPois
}
  \item{xobs}{
vector of observations
}
}
\details{
The output of this function contains posterior point estimates for all parameters of the reparameterized Poisson mixture model.
It summarizes unimodal MCMC samples by computing
measures of centrality, including mean and median, while multimodal outputs require a preprocessing, due to the
label switching phenomenon (Jasra et al., 2005). The summary measures are then computed after performing a multi-dimensional k-means clustering (Hartigan and Wong, 1979) following the suggestion of Fruhwirth-Schnatter (2006).

}
\value{
 \item{lambda }{vector of mean and median of simulated draws from the conditional posterior of the mixture model mean}
  \item{gamma.i }{vector of mean and median of simulated draws from the conditional posterior of the component mean hyperparameters; \eqn{i=1, \ldots, k}}
  \item{weight.i }{vector of mean and median of simulated draws from the conditional posterior of the component weights of the mixture distribution; \eqn{i=1, \ldots, k}}
  \item{lambda.i }{vector of mean and median of simulated draws from the conditional posterior of the component means of the mixture distribution; \eqn{i=1, \ldots, k}}
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
For multimodal outputs such as the mixture model weights, component means, and component mean hyperparameters, for each MCMC draw,
first the labels of the weights \eqn{p_i, i=1, \ldots, k} and corresponding component means are permuted in
such a way that \eqn{p_1\le \ldots \le p_k}. Then the posterior component means are partitioned into \eqn{k}
clusters by applying a standard k-means algorithm with \eqn{k} clusters, following Fruhwirth-Schnatter (2006) method. The obtained classification sequence was then used to reorder and identify the other component-specific parameters, namely component mean hyperparameters and weights.
For each group,  cluster centers are considered as parameter
estimates. 
}
\seealso{
\code{\link[Ultimixt]{K.MixPois}}
}
\examples{
N=500
U =runif(N)                                            
xobs = rep(NA,N)
for(i in 1:N){
    if(U[i]<.6){
        xobs[i] = rpois(1,lambda=1)
    }else{
        xobs[i] = rpois(1,lambda=5)
    }
}
#estimate=K.MixPois(xobs, k=2, alpha0=.5, alpha=.5, Nsim=10000)
#SM.MixPois(estimate, xobs)
#plot(estimate[[8]][,1],estimate[[2]][,1],pch=19,col="skyblue",cex=0.5,xlab="lambda",ylab="p")
#points(estimate[[8]][,2], estimate[[2]][,2], pch=19, col="gold", cex=0.5)
#points(c(1,5), c(0.6,0.4), pch=19, cex=1)
}
\keyword{summary statistics}
\keyword{mixture parameters}
\keyword{k-means clustering method}
