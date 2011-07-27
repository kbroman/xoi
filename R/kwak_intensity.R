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

intensity <- function(xomat, window, marker, n_ind, N)
{

  n_xo <- ncol(xomat)
  n_pos <- length(marker)

  center = NULL
  inc = ( marker[n_pos] - marker[1] ) / (N + 1)
  pt = marker[1] + inc
  
  for(i in 1:N)
    {
        center <- c(center, pt)
        pt = pt + inc
    }

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
               PACKAGE="xoi")
  return(list(val = output$intensity/n_ind/window*100, center = center, window_size = window))
}

# end of intensity.R
