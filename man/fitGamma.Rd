\name{fitGamma}
\alias{fitGamma}

\title{Fit Gamma model}

\description{
  Fit the gamma model for crossover interference to data on crossover
  locations. 
}

\usage{
fitGamma(d, censor, nu, lo, hi, 
         se=FALSE, supint=FALSE, rescale=FALSE,
         drop=1.5, tol=1e-5, maxit=1000, max.conv=25, 
         integr.tol=1e-8, max.subd=1000, min.subd=10,
         h=0.1, hstep=1.5)

}
\arguments{
 \item{d}{A vector of inter-crossover distances in cM.  This should
   include distances from start of chromosome to first crossover, last
   crossover to end of chromosome, and chromosome length, if there are
   no crossovers.

   Alternatively, this may be a matrix with the first column being the
   distances and second column being the censoring types (\code{censor}).}
\item{censor}{A vector of the same length as \code{d}, indicating the
  censoring type for each distance.  \code{0} = uncensored, \code{1} =
  right-censored, \code{2} = initial crossover on chromosome, \code{3} = whole
  chromosome.}
\item{nu}{A vector of interference parameters (\eqn{\nu}{nu}) at which
  to calculate the log likelihood.  If missing, \code{lo} and \code{hi}
  must be specified.}
\item{lo}{If \code{nu} is unspecified, \code{lo} 
  indicates the lower value of the interval in which to search
  for the MLE.  If \code{supint=TRUE}, this should be below the lower
  limit of the support interval.}
\item{hi}{If \code{nu} is unspecified, \code{hi}
  indicates the upper value of the interval in which to search
  for the MLE.  If \code{supint=TRUE}, this should be above the upper
  limit of the support interval.}
\item{se}{If TRUE and \code{nu} was not specified, an estimated SE
  (based on the second derivative of the log likelihood) is estimated.}
\item{supint}{If TRUE and \code{nu} was not specified, a likelihood
  support interval is calculated, with \code{drop} being the amount to
  drop in log (base 10).}
\item{rescale}{If TRUE and \code{nu} was specified, re-scale the log
  likelihoods so that the maximum is at 0.}
\item{drop}{If \code{supint} was specified, this indicates the amount to
  drop in log (base 10) for the likelihood support interval.}
\item{tol}{Tolerance for converence to calculate the likelihood, SE, and
  likelihood support interval.}
\item{maxit}{Maximum number of iterations in estimating the SE and
  likelihood support interval.}
\item{max.conv}{Maximum limit for summation in the convolutions to get
  inter-crossover distance distribution from the inter-chiasma distance
  distributions.  This should be greater than the maximum number of
  chiasmata on the 4-strand bundle.}
\item{integr.tol}{Tolerance for convergence of numerical integration.}
\item{max.subd}{Maximum number of subdivisions in numerical
  integration.}
\item{min.subd}{Minimum number of subdivisions in numerical
  integration.}
\item{h}{Step used in estimating the second derivative of the log likelihood.}
\item{hstep}{factor by which \code{h} is decreased in each iteration of
  the estimation of the second derivative of the log likelihood.}
}

\value{
  If \code{nu} is specified, we return a data frame with two columns:
  \code{nu} and the corresponding log (base e) likelihood.  If
  \code{rescale=TRUE}, the maximum log likelihood is subtracted off, so
  that its maximum is at 0.

  If \code{lo} and \code{hi} is specified, the output contains a single
  row with the MLE of \eqn{\nu}{nu} and the corresponding log
  likelihood.  If \code{se=TRUE}, we also include the estimated SE.  If
  \code{supint=TRUE}, we include two additional rows with the lower and
  upper limits of the likelihood support interval.
}

\details{
  See Broman and Weber (2000) for details of the method.

  We use R's \code{\link[stats]{integrate}} function for numerical
  integrals, \code{\link[stats]{optimize}} for optimizing the
  likelihood, and \code{\link[stats]{uniroot}} for identifying the
  endpoints of the likelihood support interval.
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
data(bssbsb)
\dontshow{bssbsb <- bssbsb[,1:50]}

xodist <- convertxoloc(find.breaks(bssbsb, chr=1))

# plot a rough log likelihood curve
\dontrun{out <- fitGamma(xodist, nu=seq(1, 19, by=2))}
\dontshow{out <- fitGamma(xodist, nu=seq(1, 19, by=2), tol=0.001)}
plot(out, type="l", lwd=2)

# get MLE
\dontrun{mle <- fitGamma(xodist, lo=8, hi=12)}
\dontshow{mle <- fitGamma(xodist, lo=8, hi=12, tol=0.001)}
mle

abline(v=mle[1], h=mle[2], col="blue", lty=2)

# get MLE and SE
\dontrun{mle <- fitGamma(xodist, lo=9.5, hi=10.5, se=TRUE)}
\dontshow{mle <- fitGamma(xodist, lo=9.5, hi=10.5, se=TRUE, tol=0.001)}
mle

# get MLE and 10^1.5 support interval
\dontrun{int <- fitGamma(xodist, lo=1, hi=20, supint=TRUE)}
\dontshow{int <- fitGamma(xodist, lo=1, hi=20, supint=TRUE, tol=0.001)}
int
abline(v=mle[2:3,1], h=mle[2:3,2], col="red", lty=2)
}

\seealso{ \code{\link[qtl]{fitstahl}} } 

\keyword{models}
