#**********************************************************************
#* 
#* intensity.R
#* 
#* Input
#*   xomat : xo matrix
#*   window : window size
#*   marker : marker vector
#*   n_ind : # of observations
#*   N : # of marker positions that we want to calculate
#*
#* Return : 
#*   val : intensity fcn value at each of marker positions
#*   center : marker position vector (It have length N)
#*   window_size : window size 
#*
#**********************************************************************/

intensity <-
function(cross, chr, window=2.5, ncalc=500)
{
  if(!missing(chr)) {
    cross <- subset(cross, chr)
    if(nchr(cross) > 1)
      warning("Considering only chr ", names(cross$geno)[1])
  }

  if(class(cross)[1] != "bc")
    stop("coincidence() currently working only for a backcross.")

  g <- cross$geno[[1]]$data
  g[is.na(g)] <- 0
  map <- cross$geno[[1]]$map

  xomat <- identify_xo(g)$xomat
  n.ind <- nind(cross)

  intensitySUB(xomat, window, map, n.ind, ncalc)
}


intensitySUB <- 
function(xomat, window, marker, n_ind, N)
{

  n_xo <- ncol(xomat)
  n_pos <- length(marker)

  center <- seq(marker[1], marker[n_pos], length=N)

  n_center <- length(center)
  xovec <- c(xomat)

  output <- .C("R_get_intensity",
               as.integer(xovec),
               as.double(window),
               as.double(center),
               as.integer(n_pos),
               as.integer(n_xo),
               as.integer(n_center),
               as.double(marker),
               intensity=as.double(rep(0,n_center)),
               as.integer(n_ind),
               PACKAGE="xoi")

  result <- data.frame(position=center,
                       intensity=output$intensity)
  attr(result, "window") <- window
  result
}

# end of intensity.R
