######################################################################
# util.R
#
# copyright (c) 1999-2012, Karl W Broman
#
# last modified Nov, 2012
# first written ~Jun, 1999
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
# Contains: find.breaks, find.breaks.F2, inferxoloc.F2, countxo,
#           convertxoloc, addlog
#
######################################################################

# find breakpoints
find.breaks <-
function(cross, chr)
{
  if(length(class(cross)) < 2 || class(cross)[2] != "cross")
    stop("Input should have class \"cross\".")

  type <- class(cross)[1]

  if(!missing(chr)) cross <- subset(cross, chr=chr)
  
  if(type == "f2") return(find.breaks.F2(cross))

  if(type != "bc" && type != "risib" && type != "riself") 
    stop("This works only for a backcross or RIL.")

  v <- vector("list", nchr(cross))
  thechr <- names(v) <- names(cross$geno)
  L <- chrlen(cross)
  for(i in seq(along=thechr)) {
    v[[i]] <- locateXO(subset(cross, chr=thechr[i]))
    attr(v[[i]], "L") <- L[i]
  }

  if(length(v)==1) return(v[[1]])
  v
}             

# find breakpoints in F2
find.breaks.F2 <-
function(cross)
{
  v <- vector("list", nchr(cross))
  names(v) <- thechr <- names(cross$geno)
  L <- chrlen(cross)
  for(i in seq(along=thechr)) {
    v[[i]] <- lapply(locateXO(subset(cross, chr=thechr[i]), full.info=TRUE),
                     inferxoloc.F2)
    attr(v[[i]], "L") <- L[i]
  }

  if(length(v)==1) return(v[[1]])
  v
}

# infer strand-specific XO locations in F2
inferxoloc.F2 <-
function(fullxoinfo)
{
  # no XOs
  if(length(fullxoinfo)==0) return(list(list(numeric(0), numeric(0))))

  # 1 XO
  if(nrow(fullxoinfo) == 1) return(list(list(fullxoinfo[1,1], numeric(0))))

  # drop extraneous rows
  fullxoinfo <- fullxoinfo[c(TRUE, fullxoinfo[-nrow(fullxoinfo),4] != fullxoinfo[-1,4]),, drop=FALSE]

  # make sure we have midpoints
  fullxoinfo[,1] <- (fullxoinfo[,2]+fullxoinfo[,3])/2

  xo <- fullxoinfo[,1]
  gleft <- fullxoinfo[,6]
  gright <- fullxoinfo[,7]

  if(gleft[1]==gright[1]+2 || gleft[1]+2==gright[1]) {
    result <- list(list(xo[1], xo[1]))
    last <- list(0)
  }
  else {
    result <- list(list(xo[1], numeric(0)))
    last <- list(1)
  }

  if(nrow(fullxoinfo)==1) return(result)

  for(i in 2:nrow(fullxoinfo)) {
    if(gleft[i]==gright[i]+2 || gleft[i]+2==gright[i]) { # A-B or B-A
      for(j in seq(along=result)) {
        result[[j]][[1]] <- c(result[[j]][[1]], xo[i])
        result[[j]][[2]] <- c(result[[j]][[2]], xo[i])
      }
    }

    else if(gleft[i]==2 && gright[i]==gleft[i-1]) { # A-H-A or B-H-B
      for(j in seq(along=result)) {
        result[[j]][[last[[j]]]] <- c(result[[j]][[last[[j]]]], xo[i])
      }
    }

    else if(gleft[i]==2 && gright[i]!=gleft[i-1]) { # A-H-B or B-H-A
      for(j in seq(along=result)) {
        last[[j]] <- 3 - last[[j]] # put on opposite strand
        result[[j]][[last[[j]]]] <- c(result[[j]][[last[[j]]]], xo[i])
      }
    }

    else if(gleft[i-1]==2 && gright[i]==2) { # H-B-H or H-A-H
      result2add <- result
      last2add <- last
      for(j in seq(along=result)) {
        result[[j]][[1]] <- c(result[[j]][[1]], xo[i])
        last[[j]] <- 1

        result2add[[j]][[2]] <- c(result2add[[j]][[2]], xo[i])
        last2add[[j]] <- 2
      }

      result <- c(result, result2add)
      last <- c(last, last2add)
    }

    else if(gright[i] == 2) { # A-B-H or B-A-H
      if(any(gleft[1:i] == 2)) { # was a previous H; consider both possibilities
        result2add <- result
        last2add <- last
        for(j in seq(along=result)) {
          result[[j]][[1]] <- c(result[[j]][[1]], xo[i])
          last[[j]] <- 1
          
          result2add[[j]][[2]] <- c(result2add[[j]][[2]], xo[i])
          last2add[[j]] <- 2
        }

        result <- c(result, result2add)
        last <- c(last, last2add)
      }
      else { # arbitrary
        for(j in seq(along=result)) {
          result[[j]][[1]] <- c(result[[j]][[1]], xo[i])
          last[[j]] <- 1
        }
      }
    }

  } # end loop over rows in fullxoinfo

  result
}




# count number of crossovers
countxo <-
function(cross, chr)
{
  if(length(class(cross)) < 2 || class(cross)[2] != "cross")
    stop("Input should have class \"cross\".")

  type <- class(cross)[1]
  if(type != "bc" && type != "risib" && type != "riself") 
    stop("This works only for a backcross or RIL.")

  if(!missing(chr)) cross <- subset(cross, chr=chr)

  br <- find.breaks(cross)

  if(!is.list(br[[1]])) return(sapply(br, length))
  apply(sapply(br, sapply, length),1,sum)
}

convertxoloc <-
function(breaks)
{
  f <- function(x, L) {
    if(length(x)==0) return(rbind(L,3))
    else {
      d <- diff(c(0,x,L))
      cen <- c(2, rep(0,length(x)-1), 1)
      return(rbind(d,cen))
    } }

  if(is.list(breaks[[1]])) {
    v <- vector("list", length(breaks))
    names(v) <- names(breaks)
    for(i in 1:length(breaks)) {
      v[[i]] <- lapply(breaks[[i]], f, attr(breaks[[i]], "L"))
      v[[i]] <- matrix(unlist(v[[i]]), ncol=2, byrow=TRUE)
    }
    for(i in 2:length(v))
      v[[1]] <- rbind(v[[1]],v[[i]])
    v <- v[[1]]
  }
  else {
    v <- lapply(breaks, f, attr(breaks, "L"))
    v <- matrix(unlist(v), ncol=2, byrow=TRUE)
  }
  v <- as.data.frame(v)
  names(v) <- c("distance", "censor")
  v
}  


# addlog: calculates log(sum(exp(input))
addlog <- function(..., threshold=200)
{
  a <- unlist(list(...))
  if(length(a)<=1) return(a)
  x <- a[1]
  a <- a[-1]
  for(i in seq(along=a)) {
    if(a[i] > x + threshold) x <- a[i]
    else if(x < a[i] + threshold)
      x <- a[i] + log1p(exp(x-a[i]))
  }
  x
}


# end of util.R
