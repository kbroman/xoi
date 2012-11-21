\name{stahlLoglik}
\alias{stahlLoglik}

\title{Calculate log likelihood for Stahl model}

\description{
  Calculate the log likelihood for the Stahl model for varying
  parameters, with data on crossover locations.
}

\usage{
stahlLoglik(xoloc, chrlen, nu, p, max.conv=25, 
         integr.tol=1e-8, max.subd=1000, min.subd=10)

}
\arguments{
\item{xoloc}{A list of crossover locations (in cM), each component being a
  vector of locations for a different meiotic product.}
\item{chrlen}{Chromosome length (in cM), either of length 1 or the same length
  as \code{xoloc}.}
\item{nu}{A vector of interference parameters (\eqn{\nu}{nu}) at which
  to calculate the log likelihood.}
\item{p}{A vector of parameter values for the proportion of crossovers
  from the no interference pathway.}
\item{max.conv}{Maximum limit for summation in the convolutions to get
  inter-crossover distance distribution from the inter-chiasma distance
  distributions.  This should be greater than the maximum number of
  chiasmata on the 4-strand bundle.}
\item{integr.tol}{Tolerance for convergence of numerical integration.}
\item{max.subd}{Maximum number of subdivisions in numerical
  integration.}
\item{min.subd}{Minimum number of subdivisions in numerical
  integration.}
}

\value{
  A vector of log likelihoods.

  The corresponding values of nu and p are saved as attributes.
}

\details{
  See Housworth and Stahl (2003) and Broman and Weber (2000) for details
  of the method. 

  If neither \code{nu} nor \code{p} has length 1, they both must have
  the same length.  If one has length 1 and the other does not, the one
  with length 1 is repeated so that they both have the same length.
}

\references{
  Housworth, E. A. and Stahl, F. W. (2003) Crossover interference in
  humans. \emph{Am. J. Hum. Genet.} \bold{73}, 188--197.
  
  Broman,  K. W. and Weber, J. L. (2000) Characterization of human
  crossover interference. \emph{Am. J. Hum. Genet.} \bold{66},
  1911--1926. 
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu}}

\examples{
data(bssbsb)
xoloc <- find.breaks(bssbsb, chr=1)

loglik <- stahlLoglik(xoloc, nu=4, p=c(0.05, 0.1, 0.15))
}

\seealso{ \code{\link[qtl]{fitstahl}} } 

\keyword{models}
