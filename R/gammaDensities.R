######################################################################
# gammaDensities.R
#
# copyright (c) 1999-2007, Karl W Broman
#
# last modified Apr, 2007
# first written ~Jan, 2001
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
# Contains: location.given.one, first.given.two, distance.given.two,
#           joint.given.two, xoprob, ioden, firstden, gammacoi,
#           stahlcoi
#
######################################################################

# dist'n of location of XO given exactly one XO 
location.given.one <-
function(nu, L=103, x, n=400, max.conv=25,
         integr.tol=1e-8, max.subd=1000, min.subd=10)
{
  if(nu <= 0) stop("nu should be positive.")

  if(missing(x)) {
    x <- seq(0,L,length=n+1)
    x <- x[-1]-x[2]/2
  }
  if(any(x < 0 || x > L)) {
    x <- x[x >= 0 & x <= L]
    warning("Dropping values outside of [0,L].")
  }
  x <- x/100
  
  result <- .C("location_given_one",
               as.double(nu),
               as.double(x),
               y=as.double(rep(0,length(x))),
               as.integer(length(x)),
               as.double(L/100),
               as.integer(max.conv),
               as.double(integr.tol),
               as.integer(max.subd),
               as.integer(min.subd),
               PACKAGE="xoi")

  data.frame(x=x*100,f=result$y/100)
}


# dist'n of location of first XO given exactly two XOs 
first.given.two <-
function(nu, L=103, x, n=400, max.conv=25,
         integr.tol=1e-8, max.subd=1000, min.subd=10)
{
  if(nu <= 0) stop("nu should be positive.")

  if(missing(x)) {
    x <- seq(0,L,length=n+1)
    x <- x[-1]-x[2]/2
  }
  if(any(x < 0 || x > L)) {
    x <- x[x >= 0 & x <= L]
    warning("Dropping values outside of [0,L].")
  }
  x <- x/100
  
  result <- .C("first_given_two",
               as.double(nu),
               as.double(L/100),
               as.double(x),
               y=as.double(rep(0,length(x))),
               as.integer(length(x)),
               as.integer(max.conv),
               as.double(integr.tol),
               as.integer(max.subd),
               as.integer(min.subd),
               PACKAGE="xoi")

  data.frame(x=x*100, f=result$y/100)
}

  
# dist'n of distance between XOs given two XOs 
distance.given.two <-
function(nu, L=103, x, n=400, max.conv=25,
         integr.tol=1e-8, max.subd=1000, min.subd=10)
{
  if(nu <= 0) stop("nu should be positive.")

  if(missing(x)) {
    x <- seq(0,L,length=n+1)
    x <- x[-1]-x[2]/2
  }
  if(any(x < 0 || x > L)) {
    x <- x[x >= 0 & x <= L]
    warning("Dropping values outside of [0,L].")
  }
  x <- x/100
  
  result <- .C("distance_given_two",
               as.double(nu),
               as.double(L/100),
               as.double(x),
               y=as.double(rep(0,length(x))),
               as.integer(length(x)),
               as.integer(max.conv),
               as.double(integr.tol),
               as.integer(max.subd),
               as.integer(min.subd),
               PACKAGE="xoi")

  data.frame(x=x*100, f=result$y/100)
}

# joint dist'n of locations of XOs given two XOs 
joint.given.two <-
function(nu, L=103, x, y, n=20, max.conv=25,
         integr.tol=1e-8, max.subd=1000, min.subd=10)
{
  if(nu <= 0) stop("nu should be positive.")

  if(missing(x) && missing(y)) { # compute on a grid
    x <- seq(0,L, length=n+1)
    x <- x[-1]-x[2]/2
    y <- x
    x <- rep(x, n)
    y <- rep(y, rep(n,n))
  }
  else if(missing(x)) { 
    x <- seq(0,L,length=n+1)
    x <- x[-1]-x[2]/2
    m <- length(x)
    x <- rep(x, length(y))
    y <- rep(y, rep(m,length(y)))
  }
  else if(missing(y)) {
    y <- seq(0,L,length=n+1)
    y <- y[-1]-y[2]/2
    m <- length(x)
    x <- rep(x, length(y))
    y <- rep(y, rep(m,length(y)))
  }
  else if(length(x) != length(y)) { # both x and y given, but not same length
    m <- length(x)
    x <- rep(x, length(y))
    y <- rep(y, rep(m,length(y)))
  }
    
  if(any(x < 0 | x > L | y < 0 | y > L)) {
    w <- (x >= 0 & x <= L & y>=0 & y <= L)
    x <- x[w]
    y <- y[w]
    warning("Dropping values outside of [0,L].")
  }

  # only include triangle
  if(any(x>y)) {
    w <- (x <= y)
    x <- x[w]
    y <- y[w]
  }

  x <- x/100
  y <- y/100

  result <- .C("joint_given_two",
               as.double(nu),
               as.double(L/100),
               as.double(x),
               as.double(y),
               z=as.double(rep(0,length(x))),
               as.integer(length(x)),
               as.integer(max.conv),
               as.double(integr.tol),
               as.integer(max.subd),
               as.integer(min.subd),
               PACKAGE="xoi")

  data.frame(x=x*100, y=y*100, f=result$z/10000)
}

# calculate probabilities of 0, 1, 2 and >2 XOs 
xoprob <-
function(nu, L=103, max.conv=25,
         integr.tol=1e-8, max.subd=1000, min.subd=10)
{
  if(nu <= 0) stop("nu should be positive.")

  result <- .C("xoprob",
               as.double(nu),
               as.double(L/100),
               prob=as.double(rep(0,4)),
               as.integer(max.conv),
               as.double(integr.tol),
               as.integer(max.subd),
               as.integer(min.subd),
               PACKAGE="xoi")

  result$prob
}

# inter-crossover density 
ioden <- 
function(nu, L=103, x, n=400, max.conv=25)
{
  if(nu <= 0) stop("nu should be positive.")

  if(missing(x)) {
    x <- seq(0,L,length=n+1)
    x <- x[-1]-x[2]/2
  }
  if(any(x < 0)) {
    x <- x[x >= 0]
    warning("Dropping values < 0")
  }
  x <- x/100
  
  result <- .C("ioden",
               as.double(nu),
               as.double(x),
               y=as.double(rep(0,length(x))),
               as.integer(length(x)),
               as.integer(max.conv),
               PACKAGE="xoi")

  y <- result$y

  data.frame(x=x*100, f=y/100)
}

# density for location of first XO
firstden <- 
function(nu, L=103, x, n=400, max.conv=25)
{
  if(nu <= 0) stop("nu should be positive.")

  if(missing(x)) {
    x <- seq(0,L,length=n+1)
    x <- x[-1]-x[2]/2
  }
  if(any(x < 0)) {
    x <- x[x >= 0]
    warning("Dropping values < 0")
  }
  x <- x/100
  
  result <- .C("firstden",
               as.double(nu),
               as.double(x),
               y=as.double(rep(0,length(x))),
               as.integer(length(x)),
               as.integer(max.conv),
               PACKAGE="xoi")

  data.frame(x=x*100, f=result$y/100)
}

# coincidence function for gamma model
gammacoi <-
function(nu, L=103, x, n=400, max.conv=25)
{
  if(nu <= 0) stop("nu should be positive.")

  if(missing(x)) {
    x <- seq(0,L,length=n+1)
    x <- x[-1]-x[2]/2
  }
  if(any(x < 0)) {
    x <- x[x >= 0]
    warning("Dropping values < 0")
  }
  x <- x/100

  result <- .C("GammaCoincidence",
               as.double(nu),
               as.double(x),
               y=as.double(rep(0,length(x))),
               as.integer(length(x)),
               as.integer(max.conv),
               PACKAGE="xoi")

  data.frame(x=x*100, coincidence=result$y)
}  

# coincidence function for gamma model
stahlcoi <-
function(nu, p=0, L=103, x, n=400, max.conv=25)
{
  if(nu <= 0) stop("nu should be positive.")
  if(p < 0 || p > 1) stop("p should be between 0 and 1.")

  if(missing(x)) {
    x <- seq(0,L,length=n+1)
    x <- x[-1]-x[2]/2
  }
  if(any(x < 0)) {
    x <- x[x >= 0]
    warning("Dropping values < 0")
  }
  x <- x/100

  result <- .C("StahlCoincidence",
               as.double(nu),
               as.double(p),
               as.double(x),
               y=as.double(rep(0,length(x))),
               as.integer(length(x)),
               as.integer(max.conv),
               PACKAGE="xoi")

  data.frame(x=x*100, coincidence=result$y)
}  



# end of gammaDensities.R
