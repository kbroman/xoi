#**********************************************************************
#* 
#* coincidence.R
#*
#* Input
#*   xomat : xo matrix
#*   window : window size
#*   marker : marker vector
#*   n_ind : # of observations
#*   N : # of marker positions that we want to calculate
#*
#* Return
#*   coincidence : coincidence value vector
#*   center : location vector 
#*   window_size : window size
#*
#**********************************************************************/

coincidence <-
function(cross, chr, window=5, ncalc=500)
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

  coincidenceSUB(xomat, window, map, n.ind, ncalc)
}


coincidenceSUB <-
function(xomat, window, marker, n_ind, N)
{

   inten <- intensitySUB(xomat, window, marker, n_ind, N)
   int_dat <- inten$intensity*window/100
   center = inten$position

   n_pos <- length(marker)
   n_center <- length(center)
   n_xo <- ncol(xomat)

   start_d <- min(which(center >= window/2))-1

   output <- .C("R_get_coincidence",
                   as.integer(xomat),
                   as.double(int_dat),
                   as.double(window),
                   as.double(center),
                   as.integer(n_xo),
                   as.integer(n_pos),
                   as.integer(n_center),
                   as.integer(start_d),
                   as.double(marker),
                   coincidence=as.double(rep(0,n_center)),
                   PACKAGE="xoi")

  result <- data.frame(distance=center, coincidence=output$coincidence/n_ind)
  attr(result, "window") <- window
  result[!is.na(result[,2]),]
}

# end of coincidence.R
