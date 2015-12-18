\name{Plot.MixReparametrized}
\alias{Plot.MixReparametrized}
\title{
plot of the MCMC output produced by K.MixReparametrized 
}
\description{
This is a generic function for a graphical rendering of the MCMC samples produced by K.MixReparametrized function. The
function draws boxplots for unimodal variables and for multimodal arguments after clustering them by applying a k-means
algorithm. It also plots line charts for other variables.
}
\usage{
Plot.MixReparametrized(xobs, estimate)
}
\arguments{
  \item{xobs}{
vector of the observations
}
  \item{estimate}{
output of the K. MixReparametrized function
}
}
\details{
Boxplots are produced using the boxplot.default method.
}
\value{
The output of this function consists of 
\item{boxplot}{ three boxplots for the radial coordinates, the mean and the standard deviation of the mixture distribution, \eqn{k} boxplots for each of the mixture model weights,  component means and component standard deviations.} 
\item{histogram}{ an histogram of the observations against an overlaid curve of the density estimate, obtained by averaging over all mixtures corresponding to the MCMC draws,}
\item{line chart}{ line charts that report the evolution of the proposal scales and of the acceptance rates over the number of batch of 50 iterations.}  
}
\references{
Kamary, K., Lee, J.Y., and Robert, C.P. (2015) Non-informative reparameterisation of location-scale mixtures. arXiv.
}
\author{
Kaniav Kamary
}
\note{
	The mixture density estimate is based on the draws simulated of the parameters obtained by K.MixReparametrized function.
}
\seealso{
\code{\link[Ultimixt]{K.MixReparametrized}}
}
\examples{
data(faithful)
xobs=faithful[50:100,1]
#estimate=K.MixReparametrized(xobs, k=2, alpha0=.5, alpha=.5, Nsim=1e4)
#plo=Plot.MixReparametrized(xobs, estimate)
}
\keyword{plot}
\keyword{density curve}