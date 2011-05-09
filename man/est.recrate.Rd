\name{est.recrate}
\alias{est.recrate}

\title{Estimate recombination rate}

\description{
  Obtain a smoothed estimate of the recombination rate along a
  chromosome, using the cM and Mbp position of markers.
}

\usage{
est.recrate(genmap, phymap, pos, window=5)
}

\arguments{
  \item{genmap}{Vector of cM positions of markers.}
  \item{phymap}{Vector of Mbp positions of markers; same length as
    \code{genmap}.}
  \item{pos}{Vector of positions at which the recombination rate should
    be estimated.  If missing, we use the physical marker positions
    plus a grid with 4 positions per Mbp.}
  \item{window}{Length of sliding window (in Mbp).}
}

\value{
  A data.frame containing the positions and estimate recombination
  rates. 
}

\details{
  We assume constant recombination rate within each marker interval. 
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu}}

%\examples{
%}

\seealso{ \code{\link{est.coi}} }


\keyword{models}
