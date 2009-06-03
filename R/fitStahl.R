######################################################################
# fitStahl.R
#
# copyright (c) 2009, Karl W Broman
#
# last modified Jun, 2009
# first written Jun, 2009
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
# Contains: stahlLoglik, fitStahl
#
######################################################################

######################################################################
# stahlLoglik
#
# xoloc =  list of crossover locations; each component being a
#          different meiotic product
#          
# chrlen = chromosome lengths, either of length 1 or of the same 
#          length as xoloc
#
# nu     = interference parameter
#
# p      = proportion of crossovers coming from the no interference 
#          pathway
#
# max.conv   = maximum number of convolutions
# integr.tol = tolerance for integration
# max.subd   = maximum number of subdivisions in integration
# min.subd   = minimum number of subdivisions in integration
#
######################################################################
stahlLoglik <-
function(xoloc, chrlen, nu, p,
         max.conv=25, integr.tol=1e-8, max.subd=1000, min.subd=10)
{
  if(!is.list(xoloc)) xoloc <- list(xoloc)
  if(length(chrlen) == 1) chrlen <- rep(chrlen, length(xoloc))
  else if(length(chrlen) != length(xoloc))
    stop("chrlen should have length 1 or the same as length(xoloc).")

  flag <- 0
  for(i in seq(along=xoloc)) {
    if(any(xoloc[[i]] < 0 | xoloc[[i]] > chrlen[i])) {
      flag <- 1
      break
    }
  }
  if(flag) stop("xoloc should be between 0 and chrlen.")

  # convert to Morgans
  xoloc <- lapply(xoloc, function(a) a/100)
  chrlen <- chrlen/100

  if(length(nu)==1) {
    if(length(p) > 1) nu <- rep(nu, length(p))
  }
  else {
    if(length(p)==1) p <- rep(p, length(nu))
    else if(length(p) != length(nu))
      stop("nu and p should be the same length (though either can have length 1).")
  }

  n <- sapply(xoloc, length)
  n.nu <- length(nu)

  loglik <- .C("R_stahl_loglik",
               as.integer(length(xoloc)),
               as.integer(n),
               as.double(unlist(xoloc)),
               as.double(chrlen),
               as.integer(n.nu),
               as.double(nu),
               as.double(p),
               loglik=as.double(rep(0,n.nu)),
               as.integer(max.conv),
               as.double(integr.tol),
               as.integer(max.subd),
               as.integer(min.subd),
               PACKAGE="xoi")

  data.frame(nu=nu, p=p, loglik=loglik$loglik)
}

######################################################################
# the same, but with the arguments reordered and with p and nu stuck
# together, and assuming they have length 1 
######################################################################
fitStahl.sub <-
function(param, xoloc, chrlen, max.conv=25, integr.tol=1e-8,
         max.subd=1000, min.subd=10)
{
  if(param[1] < 0 || param[2] < 0 || param[2] > 1)
    return(Inf)

  -stahlLoglik(xoloc, chrlen, param[1], param[2], max.conv,
               integr.tol, max.subd, min.subd)[3]
}

# here, to optimize for the model with p=0
fitStahl.sub2 <-
function(nu, xoloc, chrlen, max.conv=25, integr.tol=1e-8,
         max.subd=1000, min.subd=10)
{
  if(nu < 0) return(Inf)

  -stahlLoglik(xoloc, chrlen, nu, 0, max.conv,
               integr.tol, max.subd, min.subd)[,3]
}

# function to optimize for the Stahl model
fitStahl <-
function(xoloc, chrlen, nu=c(1,20), p=0.02, max.conv=25, integr.tol=1e-8,
         max.subd=1000, min.subd=10, verbose=TRUE, ...)
{
  if(length(nu) > 2) {
    warning("nu should have length 2; using the first two values.")
    nu <- nu[1:2]
  }
  if(length(nu) != 2)
    stop("nu should have length 2.")
  if(length(p) > 1) {
    warning("p should have length 1; using the first value.")
    p <- p[1]
  }
  if(length(p) != 1)
    stop("p should have length 1.")
    
  out0 <- optimize(fitStahl.sub2, interval=c(nu[1], nu[2]), xoloc=xoloc, chrlen=chrlen,
                   max.conv=max.conv, integr.tol=integr.tol, max.subd=max.subd,
                   min.subd=min.subd)
  nu <- out0$minimum
  if(verbose) cat("For p=0, nuhat =", nu, "\n       log lik =", -out0$objective, "\n")

  if(verbose>1 && !("control" %in% names(list(...))))
    out <- optim(c(nu, p), fitStahl.sub, xoloc=xoloc, chrlen=chrlen,
                 max.conv=max.conv, integr.tol=integr.tol, max.subd=max.subd,
                 min.subd=min.subd, control=list(trace=verbose-1), ...)
  else
    out <- optim(c(nu, p), fitStahl.sub, xoloc=xoloc, chrlen=chrlen,
                 max.conv=max.conv, integr.tol=integr.tol, max.subd=max.subd,
                 min.subd=min.subd, ...)

  if(verbose)
    cat("\n  nuhat =", out$par[1], "\n   phat =", out$par[2], "\nlog lik =", -out$value, "\n")

  if(out0$objective <= out$value) {
    out <- c(out0$minimum, 0, -out0$objective, out0$minimum, -out0$objective, 0)
    if(verbose) cat("Inferred that p=0\n")
  }
  else {
    out <- c(out$par, -out$value, out0$minimum, -out0$objective, out0$objective - out$value)
    if(verbose) cat("Inferred that p>0\n")
  }
  names(out) <- c("nu", "p", "loglik", "nu0", "loglik0", "ln LR testing p=0")
  out
}



# end of fitStahl.R
