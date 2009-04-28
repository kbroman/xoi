######################################################################
# chiasma.R
#
# copyright (c) 1999-2006, Karl W Broman
#
# last modified Nov, 2006
# first written ~Jun, 1999
# Licensed under the GNU General Public License version 2 (June, 1991)
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
function(xo, max.chiasma=max(xo)*2+5, n.iter=10000, tol=1e-6)
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
