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

coincidence <- function(xomat, window, marker, n_ind, N)
{

   inten <- intensity(xomat, window, marker, n_ind, N)
   int_dat <- inten$val*window/100
   center = inten$center

   n_pos <- length(marker)
   n_center <- length(center)
   n_xo <- ncol(xomat)

   xovec <- c(xomat)

   start_d <- 1
   aa=0
   i = 1
   while(aa == 0)
   {
       if(center[i] > 3)
       {
           start_d <- i-1
           aa = 1
       }
       i = i+1
   }

   output <- .C("R_get_coincidence",
                   as.integer(xovec),
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

   return( list(coincidence = output$coincidence/n_ind, center =
center, window_size = window ))
}

# end of coincidence.R
