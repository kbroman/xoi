\name{stahlcoi}
\alias{stahlcoi}

\title{Coincidence function for the Stahl model}

\description{
  Calculates the coincidence function for the Stahl model.
}

\usage{
stahlcoi(nu, p=0, L=103, x, n=400, max.conv=25)
}
\arguments{
\item{nu}{The interference parameter in the gamma model.}
\item{p}{The proportion of chiasmata coming from the no-interference
  mechanism.}  
\item{L}{Maximal distance (in cM) at which to calculate the
  density. Ignored if \code{x} is specified.}
\item{x}{If specified, points at which to calculate the density.}
\item{n}{Number of points at which to calculate the density.  The points
  will be evenly distributed between 0 and \code{L}. Ignored if \code{x}
  is specified.}
\item{max.conv}{Maximum limit for summation in the convolution.  This
  should be greater than the maximum number of chiasmata on the 4-strand
  bundle.}
}

\value{
  A data frame with two columns: \code{x} is the distance (between 0
  and \code{L}, in cM) at which the coicidence was calculated and
  \code{coincidence}.
}

\details{
  The Stahl model is an extension to the gamma model, in which chiasmata
  occur according to two independent mechanisms.  A proportion \eqn{p}
  come from a mechanism exhibiting no interference, and a proportion
  1-\eqn{p} come from a mechanism in which chiasma locations follow a
  gamma model with interference parameter \eqn{\nu}{nu}.

  Let \eqn{f(x;\nu,\lambda)}{f(x;nu,lambda)} denote the density of a gamma random
  variable with parameters shape=\eqn{\nu}{nu} and rate=\eqn{\lambda}{lambda}.

  The coincidence function for the Stahl model is
  \eqn{C(x;\nu,p) = [p + \sum_{k=1}^{\infty} f(x;k\nu,
    2(1-p)\nu)]/2}{C(x;nu,p) = [p + sum_(k=1 to infty) 
    f(x;k*nu,2*(1-p)nu)]/2}.
}

\references{
  Copenhaver, G. P., Housworth, E. A. and Stahl, F. W. (2002) Crossover
  interference in Arabidopsis. \emph{Genetics} \bold{160}, 1631--1639. 

  Housworth, E. A. and Stahl, F. W. (2003) Crossover interference in
  humans. \emph{Am J Hum Genet} \bold{73}, 188--197.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu}}

\examples{
f1 <- stahlcoi(1, p=0, L=200)
plot(f1, type="l", lwd=2, las=1,
     ylim=c(0,1.25), yaxs="i", xaxs="i", xlim=c(0,200))

f2 <- stahlcoi(2.6, p=0, L=200)
lines(f2, col="blue", lwd=2)

f2s <- stahlcoi(2.6, p=0.1, L=200)
lines(f2s, col="blue", lwd=2, lty=2)

f3 <- stahlcoi(4.3, p=0, L=200)
lines(f3, col="red", lwd=2)

f3s <- stahlcoi(4.3, p=0.1, L=200)
lines(f3s, col="red", lwd=2, lty=2)

f4 <- stahlcoi(7.6, p=0, L=200)
lines(f4, col="green", lwd=2)

f4s <- stahlcoi(7.6, p=0.1, L=200)
lines(f4s, col="green", lwd=2, lty=2)
}

\seealso{ \code{\link{gammacoi}}, \code{\link{location.given.one}},
    \code{\link{first.given.two}},
      \code{\link{distance.given.two}},
	\code{\link{ioden}},
	\code{\link{firstden}},
	\code{\link{xoprob}} } 

\keyword{distribution}
