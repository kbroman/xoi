######################################################################
# simStahl.R
#
# copyright (c) 2006, Karl W Broman
#
# last modified Nov, 2006
# first written Nov, 2006
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
# Contains: simStahl
#
######################################################################

######################################################################
#
# simStahl
#
# This function simulates crossover locations under the Stahl model.
#
# n.sim = number of meioses to simulate
# nu    = interference parameter
# p     = proportion of chiasmata from the NI mechanism
# L     = chromosome length (in cM)
#
######################################################################


simStahl <-
function(n.sim, nu=1, p=0, L=103, n.bins4start=10000)
{
  L <- L/100
  if(nu <= 0) stop("nu should be positive.")
  if(p < 0 || p > 1) stop("p should be in [0,1].")
  if(n.sim <= 0) stop("n should be a positive integer.")
  if(L < 0) stop("L should be positive.")
  if(n.bins4start < 1000) {
    warning("n.bins4start should be large.  Using 1000.")
    n.bins4start <- 1000
  }

  max.nxo <- qpois(1-1e-10, L)*10

  out <- .C("simStahl",
            as.integer(n.sim),
            as.double(nu),
            as.double(p),
            as.double(L),
            nxo = as.integer(rep(0,n.sim)), # number of crossovers
            loc = as.double(rep(0,n.sim*max.nxo)),
            as.integer(max.nxo),
            as.integer(n.bins4start),
            PACKAGE="xoi")

  out <- lapply(as.data.frame(rbind(out$nxo, matrix(out$loc*100, nrow=max.nxo))),
                function(a) {if(a[1]==0) return(numeric(0)); a[(1:a[1])+1] })
  attr(out, "L") <- L*100
  out
}

# end of simStahl.R 
