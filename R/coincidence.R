######################################################################
# coincidence.R
#
# copyright (c) 2006-7, Karl W Broman
#
# last modified Apr, 2007
# first written Dec, 2006
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
# Contains: est.coi
#
######################################################################

est.coi <-
function(cross, chr, pos, window=0, 
         fill.method=c("imp", "argmax"), error.prob=1e-10,
         map.function=c("haldane", "kosambi", "c-f", "morgan"))
{
  if(length(class(cross)) < 2 || class(cross)[1] != "bc")
    stop("This function is only prepared for backcrosses.")

  if(!missing(chr)) {
    if(length(chr) != 1)
      stop("You should specify just one chromosome.")
    cross <- subset(cross, chr=chr)
  }
  else cross <- subset(cross, chr=1)
    
  dat <- cross$geno[[1]]$data
  if(any(is.na(dat))) {

    fill.method <- match.arg(fill.method)
    map.function <- match.arg(map.function)
    cross <- fill.geno(cross, method=fill.method, error.prob=error.prob,
                       map.function=map.function)
    dat <- cross$geno[[1]]$data
    if(any(is.na(dat))) {
      warning("Some data still missing.")
      dat[is.na(dat)] <- 0
    }
  }
    
  map <- cross$geno[[1]]$map
  if(!missing(pos)) {
    if(length(pos) != length(map))
      stop("pos must have length ", length(map))
    map <- pos
  }

  ni <- nrow(dat)
  nm <- ncol(dat)
  if(nm < 3) stop("Need at least three markers.")
  
  npair <- choose(nm-1, 2)

  out <- .C("R_est_coi",
            as.integer(ni),
            as.integer(nm),
            as.integer(npair),
            as.double(map),
            as.integer(dat),
            d=as.double(rep(0, npair)), # distances
            coi1=as.double(rep(0, npair)), # smooth top and bottom and then ratio
            coi2=as.double(rep(0, npair)), # ratio then smooth 
            n=as.integer(0), # no. distances to keep
            as.double(window),
            PACKAGE="xoi")
            
  n <- out$n
  data.frame(d=out$d[1:n], coi1=out$coi1[1:n], coi2=out$coi2[1:n])
}



# end of coincidence.R
