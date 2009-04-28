######################################################################
# recrate.R
#
# copyright (c) 2008, Karl W Broman
#
# last modified Aug, 2008
# first written Aug, 2008
# Licensed under the GNU General Public License version 2 (June, 1991)
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
    pos <- sort(c(phymap, phymap-window/2, phymap+window/2))
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
