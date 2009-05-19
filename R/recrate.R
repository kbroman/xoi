######################################################################
# recrate.R
#
# copyright (c) 2008-9, Karl W Broman
#
# last modified May, 2009
# first written Aug, 2008
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
# Contains: est.recrate
#
######################################################################

######################################################################
# est.recrate
#
# Obtain a smoothed estimate of recombination rate across
# a chromosome, using a sliding window of fixed width
#
# genmap = cM positions of markers
# phymap = Mbp positions of markers
# pos    = Mbp positions at which to calculate the cM/Mbp rate
# window = total length of sliding window (in Mbp)
######################################################################
est.recrate <-
function(genmap, phymap, pos, window=5)
{
  if(length(genmap) != length(phymap))
    stop("genmap and phymap should be the same length.")

  if(any(is.na(genmap) | is.na(phymap))) {
    drop <- is.na(genmap) | is.na(phymap)
    genmap <- genmap[!drop]
    phymap <- phymap[!drop]
    if(length(genmap) == 0)
      stop("Need multiple markers with known genetic and physical positions.")
  }

  if(length(unique(phymap)) != length(phymap))
    stop("Can't have markers right on top of each other (in physical distance).")

  if(window < 0)
    stop("window should be > 0")

  if(missing(pos)) {
    pos <- sort(c(phymap, seq(min(phymap), max(phymap), by=0.25)))
    pos <- unique(pos[pos >= min(phymap) & pos <= max(phymap)])
  }

  if(any(diff(genmap)<0))
    stop("genmap should be non-decreasing")
  if(any(diff(phymap)<0))
    stop("phymap should be non-decreasing")
  if(any(diff(pos)<0))
    stop("pos should be non-decreasing")

  if(any(pos < min(phymap) | pos > max(phymap))) {
    warning("pos should be within range of phymap")
    pos <- pos[pos >= min(phymap) & pos <= max(phymap)]
  }

  if(diff(range(phymap)) < window)
    stop("range of phymap should exceed the window length.")

  rate <- .C("R_est_recrate",
             as.integer(length(genmap)),
             as.double(genmap),
             as.double(phymap),
             as.integer(length(pos)),
             as.double(pos),
             rate=as.double(rep(0,length(pos))),
             as.double(window),
             as.double(rep(0, length(genmap)-1)),
             PACKAGE="xoi")$rate
            
  data.frame(pos=pos, rate=rate)
}

# end of recrate.R
