\name{intensity}
\alias{intensity}

\title{Get intensity function}

\description{
  Get intensity function from XO information matrix.
}

\usage{
intensity(xomat, window, marker, n_ind, N)
}
\arguments{
 \item{xomat}{XO information matrix from identify2(dataset) function. }
 \item{window}{window size for intensity. }
 \item{marker}{marker vector. }
 \item{n_ind}{Total number of index. }
 \item{N}{Total number of points to calculate intensity.  }
}

\value{
val : intensity fcn value at each of marker positions
center : marker position vector (It have length N)
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


\seealso{ \code{\link{coincidence}} }



\keyword{utilities} % keywords in R source : doc/KEYWORDS
