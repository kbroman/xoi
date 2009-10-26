\name{find.breaks}
\alias{find.breaks}

\title{Estimate crossover locations}

\description{
  Estimate the locations of crossovers in a backcross.
}

\usage{
find.breaks(cross, chr)
}
\arguments{
 \item{cross}{An object of class \code{cross}. (This must be a backcross.) See
   \code{\link[qtl]{read.cross}} for details.  }
 \item{chr}{Optional set of chromosomes on which to look for crossovers.
   If missing, all chromosomes are considered.}
}

\value{
  If only one chromosome is considered, this is a list with one
  component for each individual, with each component being
  \code{numeric(0)} if there were no crossovers or a vector giving the
  crossover locations.  The length of the chromosome (in cM) is saved as
  an attribute.  (Note that the format is the same as the output of
  \code{\link{simStahl}}.)   

  If multiple chromosomes were considered, this is a list with one
  element for each chromosome, each of which is a list with one element
  for each individual, as above.  
}

\details{
  This works only a backcross.  We use the internal function (within
  R/qtl) \code{locate.xo}.  Crossovers are estimated to be at the midpoint
  of the interval between the nearest flanking typed markers.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu}}

\examples{
data(bssbsb)

# crossover locations on chromosome 1
xoloc1 <- find.breaks(bssbsb, chr=1)

# crossover locations on all chromosomes
xoloc <- find.breaks(bssbsb)
}

\seealso{ \code{\link{convertxoloc}}, \code{\link{fitGamma}},
      \code{\link{simStahl}} }


\keyword{utilities}
