\name{location.given.one}
\alias{location.given.one}

\title{Location of crossover given there is one}

\description{
  Calculates the density of the location of the crossover on a random
  meiotic product, given that there is precisely one crossover, for the
  gamma model.
}

\usage{
location.given.one(nu, L=103, x, n=400, max.conv=25,
         integr.tol=1e-8, max.subd=1000, min.subd=10)

}
\arguments{
 \item{nu}{The interference parameter in the gamma model.}
\item{L}{The length of the chromsome in cM.}
\item{x}{If specified, points at which to calculate the density.}
\item{n}{Number of points at which to calculate the density.  The points
  will be evenly distributed between 0 and \code{L}. Ignored if \code{x}
  is specified.}
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
  A data frame with two columns: \code{x} is the location (between 0
  and \code{L}, in cM) at which the density was calculated and \code{f} is the
  density. 
}

\details{
  Let \eqn{f(x;\nu)}{f(x;nu)} denote the density of a gamma random
  variable with parameters shape=\eqn{\nu}{nu} and rate=\eqn{2\nu}{2
    nu}, and let \eqn{f_k(x;\nu)}{f_k(x;\nu)} denote the density of a gamma random
  variable with parameters shape=\eqn{k \nu}{k nu} and rate=\eqn{2\nu}{2
    nu}.  

  The distribution of the distance from one crossover to the next is
  \eqn{f^*(x;\nu) = \sum_{k=1}^{\infty} f_k(x;\nu)/2^k}{f*(x;nu) = sum_(k=1 to infty)
    f_k(x;\nu)/2^k}.

  The distribution of the distance from the start of the chromosome to
  the first crossover is \eqn{g^*(x;\nu) = 1 - F^*(x;\nu)}{g*(x;nu) = 1
    - F*(x;nu)} where \eqn{F^*}{F*} is the cdf of \eqn{f^*}{f*}.

  We calculate the distribution of the location of the crossover on a
  product with a single crossover as the convolution of \eqn{g^*}{g*}
  with itself, and then rescaled to integrate to 1 on the interval
  (0,L).  
}

\references{
  Broman,  K. W. and Weber, J. L. (2000) Characterization of human
  crossover interference. \emph{Am. J. Hum. Genet.} \bold{66},
  1911--1926. 

  Broman, K. W., Rowe, L. B., Churchill, G. A. and Paigen, K. (2002)
  Crossover interference in the mouse. \emph{Genetics} \bold{160}, 1123--1131.

  McPeek, M. S. and Speed, T. P. (1995) Modeling interference in genetic
  recombination.  \emph{Genetics} \bold{139}, 1031--1044.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu}}

\examples{
f1 <- location.given.one(1, L=200, n=201)
plot(f1, type="l", lwd=2, las=1,
     ylim=c(0,0.006), yaxs="i", xaxs="i", xlim=c(0,200))

f2 <- location.given.one(2.6, L=200, n=201)
lines(f2, col="blue", lwd=2)

f3 <- location.given.one(4.3, L=200, n=201)
lines(f3, col="red", lwd=2)

f4 <- location.given.one(7.6, L=200, n=201)
lines(f4, col="green", lwd=2)
}

\seealso{ 
    \code{\link{first.given.two}},
      \code{\link{distance.given.two}},
      \code{\link{joint.given.two}},
	\code{\link{ioden}}, \code{\link{firstden}},
	\code{\link{xoprob}},
      \code{\link{gammacoi}} } 

\keyword{distribution}
