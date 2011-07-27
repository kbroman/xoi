\name{intensity}
\alias{intensity}

\title{Estimate intensity function}

\description{
  Estimate intensity function for a chromosome.
}

\usage{
intensity(cross, chr, window=2.5, ncalc=500)
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
  Data frame with columns \code{position} and \code{intensity}.  The
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

out <- intensity(x)
plot(out, type="l", lwd=2, ylim=c(0, max(out[,2])))
}


\seealso{ \code{\link{coincidence}} }



\keyword{utilities} % keywords in R source : doc/KEYWORDS
