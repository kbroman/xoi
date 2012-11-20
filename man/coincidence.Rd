\name{coincidence}
\alias{coincidence}

\title{Estimate coincidence function}

\description{
  Estimate coincidence function for a chromosome.
}

\usage{
coincidence(cross, chr, window=5, ncalc=500)
}
\arguments{
 \item{cross}{Cross object; must be a backcross.  See
   \code{\link[qtl]{read.cross}} for format details.}
 \item{chr}{Chromosome to consider (only one is allowed).  If missing,
  the first chromosome is considered.}
 \item{window}{Window size}
 \item{ncalc}{Total number of points for calculations.}
}

\value{
  Data frame with columns \code{distance} and \code{coincidence}.  The
  input argument \code{window} is kept as an attribute.
}

%\details{
%optionally give some details.
%}

%\references{
  % if necessary
%}

\author{Il youp Kwak}

\examples{
map1 <- sim.map(103, n.mar=104, anchor=TRUE, include.x=FALSE, eq=TRUE)
x <- sim.cross(map1, n.ind=2000, m=6, type="bc")

out <- coincidence(x, ncalc=101)
plot(out, type="l", lwd=2, ylim=c(0, max(out[,2])))
}


\seealso{ \code{\link{intensity}}, \code{\link{est.coi}} }



\keyword{utilities} % keywords in R source : doc/KEYWORDS
