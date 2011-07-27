\name{coincidence}
\alias{coincidence}

\title{Estimate coincidence function}

\description{
  Estimate coincidence function for a chromosome.
}

\usage{
coincidence(cross, chr, window=5, ncalc=100)
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

%\examples{
%
%}


\seealso{ \code{\link{intensity}}, \code{\link{est.coi}} }



\keyword{utilities} % keywords in R source : doc/KEYWORDS
