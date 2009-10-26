\name{convertxoloc}
\alias{convertxoloc}

\title{Convert format of crossover locations data}

\description{
  Convert the format of data on crossover locations to that needed for
  the function \code{\link{fitGamma}.}
}

\usage{
convertxoloc(breaks)
}
\arguments{
 \item{breaks}{A list of crossover locations, as output by
   \code{\link{find.breaks}} or \code{\link{simStahl}}.}
}

\value{
  A data frame with two columns: the inter-crossover and crossover-to
  chromosome end differences (\code{"distance"}) and indicators of
  censoring type (\code{"censor"}), with 0 =
  distance between crossovers, 1=start of chromosome to first crossover,
  2 = crossover to end of chromosome, and 3 = whole chromosome.
}


\author{Karl W Broman, \email{kbroman@biostat.wisc.edu}}

\examples{
data(bssbsb)

# crossover locations on chromosome 1
xoloc1 <- convertxoloc(find.breaks(bssbsb, chr=1))

# crossover locations on all chromosomes
xoloc <- convertxoloc(find.breaks(bssbsb))
}

\seealso{ \code{\link{find.breaks}}, \code{\link{fitGamma}},
      \code{\link{simStahl}} }


\keyword{utilities}
