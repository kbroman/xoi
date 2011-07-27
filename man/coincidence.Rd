\name{coincidence}
\alias{coincidence}

\title{Get coincidence function}

\description{
  Get coincidence function from XO information matrix.
}

\usage{
coincidence(xomat, window, marker, n_ind, N)
}
\arguments{
 \item{xomat}{XO information matrix from identify2(dataset) function. }
 \item{window}{window size for intensity. }
 \item{marker}{marker vector. }
 \item{n_ind}{Total number of index. }
 \item{N}{Total number of points to calculate intensity.  }
}

\value{
coincidence : coincidence value vector
center : location vector 
window_size : window size
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


\seealso{ \code{\link{intensity}} }



\keyword{utilities} % keywords in R source : doc/KEYWORDS
