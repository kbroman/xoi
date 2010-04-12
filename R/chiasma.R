######################################################################
# chiasma.R
#
# copyright (c) 1999-2010, Karl W Broman
#
# last modified Apr, 2010
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
# Contains: chiasma
#
######################################################################

######################################################################
# given a vector of crossover counts, fits the Poisson model,
# truncated Poisson model (with obligate chiasma),
# freely varying and obligate chiasma, all assuming no chromatid
# interference
######################################################################

chiasma <-
function(xo, max.chiasma=max(xo)*2+5, n.iter=10000, tol=1e-6, verbose=FALSE)
{
  n.xo <- length(xo)

  res <- .C("chiasma",
            as.integer(xo),
            as.integer(n.xo),
            as.integer(max.chiasma),
            p.ch = as.double(rep(0,(max.chiasma+1)*4)),
            p.xo = as.double(rep(0,(max.chiasma+1)*4)),
            lambda = as.double(c(0,0)),
            as.double(rep(0,(max.chiasma+1)*n.xo+2+2*(max.chiasma+1))),
            n.iter = as.integer(c(n.iter,0,0,0,0)),
            as.double(tol),
            PACKAGE="xoi")

  if(verbose)
    cat("Done!  number of iterations = ",
        paste(as.character(res$n.iter[-1]),collapse=" "), "\n")

  xo.table <- rbind(table(factor(xo,levels=0:max.chiasma))/length(xo),
                    matrix(res$p.xo,nrow=4,byrow=TRUE))
  ch.table <- matrix(res$p.ch,nrow=4,byrow=TRUE)
  colnames(ch.table) <- 0:(ncol(ch.table)-1)

  rownames(ch.table) <- c("truncPois","Pois","oblchi","free")
  rownames(xo.table) <- c("observed", "truncPois","Pois","oblchi","free")

  out <- list(xo.table=xo.table*n.xo,ch.table=ch.table,
              lambda = res$lambda)
  attr(out, "n.iter") <-  res$n.iter[-1]
  out
}

# end of chiasma.R
