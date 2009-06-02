######################################################################
# fitGamma.R
#
# copyright (c) 1999-2006, Karl W Broman
#
# last modified Nov, 2006
# first written ~Nov, 1999
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
#
# Part of the R/xoi package
# Contains: fitGamma
#
######################################################################

######################################################################
# fitGamma
#
# d = inter-xo distances (in cM)
# censor = vector of the same length as d, indicating:
#              0 uncensored
#              1 right-censored (last point)
#              2 initial point
#              3 whole chromosome
# nu     = if specified, we calculate the log likelihood at each of
#          these values
# lo,hi  = if nu unspecified, we'll get the MLE, which is assumed to
#          be in this interval
#
# rescale    = if TRUE (and nu specified), rescale log like to have
#              maximum = 0
# max.conv = maximum number of convolutions
# integr.tol = tolerance for integration
# max.subd   = maximum number of subdivisions in integration
# min.subd   = minimum number of subdivisions in integration
#
# tol        = tolerance for converence
# maxit      = maximum number of interations
######################################################################
fitGamma <-
function(d, censor, nu, lo, hi, 
         se=FALSE, supint=FALSE, rescale=FALSE,
         drop=1.5, tol=1e-5, maxit=1000, max.conv=25, 
         integr.tol=1e-8, max.subd=1000, min.subd=10,
         h=0.1, hstep=1.5)
{
  if(missing(censor) && !is.null(d) && ncol(d)==2) {
    censor <- d[,2]
    d <- d[,1]
  }

  if(any(d <= 0 | is.na(d)))
    stop("d should be positive and not NA")
  if(any(is.na(censor) | (censor != 0 & censor != 1 &
                          censor != 2 & censor != 3)))
    stop("censor should be 0, 1, 2 or 3 and not NA")

  if(length(d) != length(censor))
    stop("d and censor should have the same length.")

  d <- d/100

  if(!missing(nu)) {
    if(!missing(lo) || !missing(hi))
      warning("lo and hi ignored")
    if(se || supint)
      warning("se and support interval not calculated when nu is specified.")

    if(any(nu <= 0 | is.na(nu)))
      stop("nu should be positive and not NA")

    result <- .C("GammaS",
                 as.integer(length(d)),
                 as.double(d),
                 as.integer(censor),
                 as.integer(length(nu)),
                 nu=as.double(nu),
                 loglik=as.double(rep(0,length(nu))),
                 as.integer(max.conv),
                 as.integer(rescale),
                 as.double(integr.tol),
                 as.integer(max.subd),
                 as.integer(min.subd),
                 PACKAGE="xoi")

    return(data.frame(nu=result$nu, loglik=result$loglik))
  }
  else {
    if(missing(lo) || missing(hi)) 
      stop("Need to specify nu or both lo and hi")
      

    if(lo < tol) lo <- tol
    if(hi < tol) hi <- tol

    if(lo < 0 || hi < 0 || lo >= hi)
      stop("Must have lo, hi positive and lo < hi")

    result <- .C("GammaMax",
                 as.integer(length(d)),
                 as.double(d),
                 as.integer(censor),
                 as.double(lo),
                 as.double(hi),
                 nu=as.double(0),
                 loglik=as.double(0),
                 as.integer(max.conv),
                 as.double(tol),
                 as.double(integr.tol),
                 as.integer(max.subd),
                 as.integer(min.subd),
                 PACKAGE="xoi")

    nu <- result$nu
    loglik <- result$loglik
    out <- data.frame(nu=nu, loglik=loglik)

    if(se) {
      seresult <- .C("GammaSE",
                     as.integer(length(d)),
                     as.double(d),
                     as.integer(censor),
                     as.double(nu),
                     se=as.double(0),
                     as.double(0), # sec deriv
                     as.integer(max.conv),
                     as.double(h),
                     as.double(hstep),
                     as.double(tol),
                     as.integer(maxit),
                     as.double(integr.tol),
                     as.integer(max.subd),
                     as.integer(min.subd),
                     PACKAGE="xoi")
      out <- cbind(out, se=seresult$se)
    }
    
    if(supint) {
      intresult <- .C("GammaInterval",
                      as.integer(length(d)),
                      as.double(d),
                      as.integer(censor),
                      as.double(lo),
                      as.double(hi),
                      as.double(nu),
                      int=as.double(rep(0,2)),
                      loglik=as.double(rep(0,2)),
                      as.double(drop*log(10)), # drop in ln lik
                      as.integer(max.conv),
                      as.double(tol),
                      as.integer(maxit),
                      as.double(integr.tol),
                      as.integer(max.subd),
                      as.integer(min.subd),
                      PACKAGE="xoi")
      if(se)
        out <- rbind(out, data.frame(nu=intresult$int,
                                     loglik=intresult$loglik,
                                     se=rep(NA,2)))
      else
        out <- rbind(out, data.frame(nu=intresult$int,
                                     loglik=intresult$loglik))
      rownames(out) <- c("est","low","high")
    }
  }
  out
}

# end of fitGamma.R
