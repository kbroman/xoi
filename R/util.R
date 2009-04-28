######################################################################
# util.R
#
# copyright (c) 1999-2008, Karl W Broman
#
# last modified Aug, 2008
# first written ~Jun, 1999
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/xoi package
# Contains: find.breaks, countxo
#
######################################################################

# find breakpoints
find.breaks <-
function(cross, chr)
{
  if(length(class(cross)) < 2 || class(cross)[2] != "cross")
    stop("Input should have class \"cross\".")

  type <- class(cross)[1]
  if(type != "bc" && type != "risib" && type != "riself") 
    stop("This works only for a backcross or RIL.")

  if(!missing(chr)) cross <- subset(cross, chr=chr)
  
  v <- vector("list", nchr(cross))
  names(v) <- names(cross$geno)
  L <- chrlen(cross)
  thechr <- names(cross$geno)
  for(i in seq(along=thechr)) {
    v[[i]] <- locate.xo(subset(cross, chr=thechr[i]))
    attr(v[[i]], "L") <- L[i]
  }

  if(length(v)==1) return(v[[1]])
  v
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
  v
}  


# end of util.R
