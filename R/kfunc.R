#####################################################################
#
# kfunc.h
#
# copyright (c) 2006, Karl W Broman
#
# last modified Nov, 2006
# first written Apr, 2006
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
# Contains: kfunc
#
######################################################################

######################################################################
# code to estimate the 1-d version of Ripley's K function
#
# x = list with sorted locations of the data
# d = values at which to calculate the function
# lengths = lengths of segments studied
# exclude = distance to exclude
######################################################################

kfunc <-
function(x, d=seq(0,100,by=0.1), lengths, exclude=0, tol=1e-6)
{
  npt <- sapply(x, length)
  if(!any(npt>0)) stop("Need to have some points.")

  x <- x[npt > 0]
  if(missing(lengths)) 
    lengths <- sapply(x, max)
  else lengths <- lengths[npt > 0]

  if(length(lengths) != length(x))
    stop("length(lengths) != length(x)")
  if(any(d < exclude)) {
    warning("some d < exclude; these excluded")
    d <- d[d > exclude]
  }

  output <- .C("R_kfunc",
               as.integer(length(x)),
               as.integer(sapply(x, length)),
               as.double(unlist(x)),
               as.double(lengths),
               as.integer(length(d)),
               as.double(d),
               as.double(exclude),
               k=as.double(rep(0,length(d))),
               area=as.double(rep(0,length(d))),
               rate = as.double(0),
               as.double(tol),
               PACKAGE="xoi")

  rate <- output$rate
  k <- output$k
  area = output$area
  result <- data.frame(d=d, k=k, se=1/sqrt(output$area * rate))
  attr(result,"rate") <- rate
  class(result) <- c("kfunc","data.frame")
  result
}
