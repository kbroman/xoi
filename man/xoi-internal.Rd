\name{xoi-internal}
\alias{kfunc}
\alias{recrate2scanone}
\alias{fitStahl.sub}
\alias{fitStahl.sub2}
\alias{get_n_xo}
\alias{identify_xo}
\alias{intensitySUB}
\alias{coincidenceSUB}

\title{Internal xoi functions}

\description{
  Internal xoi functions. These are generally not to be called by the
  user.
}

\usage{
kfunc(x, d=seq(0,100,by=0.1), lengths, exclude=0, tol=1e-6)
recrate2scanone(recrate, phymap)
fitStahl.sub(param, xoloc, chrlen, max.conv=25, integr.tol=1e-8,
             max.subd=1000, min.subd=10)
fitStahl.sub2(nu, xoloc, chrlen, max.conv=25, integr.tol=1e-8,
              max.subd=1000, min.subd=10)
get_n_xo(sdat)
identify_xo(sdat)
intensitySUB(xomat, window, marker, n_ind, N)
coincidenceSUB(xomat, window, marker, n_ind, N)
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\keyword{internal}
