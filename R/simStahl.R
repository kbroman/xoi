## simStahl.R

#' Simulate crossover locations under the Stahl model
#' 
#' Simulate crossover locations under the Stahl model.
#' 
#' The Stahl model is an extension to the gamma model, in which chiasmata occur
#' according to two independent mechanisms.  A proportion \eqn{p} come from a
#' mechanism exhibiting no interference, and a proportion 1-\eqn{p} come from a
#' mechanism in which chiasma locations follow a gamma model with interference
#' parameter \eqn{\nu}{nu}.
#' 
#' @param n.sim Number of meiotic products to simulate.
#' @param nu The interference parameter in the gamma model.
#' @param p The proportion of chiasmata coming from the no-interference
#' mechanism.
#' @param L Chromosome length (in cM).
#' @param n.bins4start We approximate the distribution of the location of the
#' first crossover from the mechanism exhibiting interference using a even grid
#' with this many bins.
#' @return A vector of length \code{n.sim}, each element being empty (for
#' products with no crossovers) or a vector of crossover locations, in cM.  An
#' attribute, \code{L}, contains the chromosome length in cM.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{fitGamma}}, \code{\link[qtl]{sim.cross}}
#' @references Copenhaver, G. P., Housworth, E. A. and Stahl, F. W. (2002)
#' Crossover interference in Arabidopsis. \emph{Genetics} \bold{160},
#' 1631--1639.
#' 
#' Housworth, E. A. and Stahl, F. W. (2003) Crossover interference in humans.
#' \emph{Am J Hum Genet} \bold{73}, 188--197.
#' @keywords datagen
#' @examples
#' 
#' # simulations with no interference
#' xoNI <- simStahl(100, nu=1, p=0, L=200)
#' 
#' # simulations under gamma model with nu=7.6
#' xogamma <- simStahl(100, nu=7.6, p=0, L=200)
#' 
#' # simulations under Stahl model with nu=7.6, p=0.1
#' xostahl <- simStahl(100, nu=7.6, p=0.1, L=200)
#' 
#' @importFrom stats qpois
#' @useDynLib xoi
#' @export
simStahl <-
function(n.sim, nu=1, p=0, L=103, n.bins4start=10000)
{
  L <- L/100
  if(nu <= 0) stop("nu should be positive.")
  if(p < 0 || p > 1) stop("p should be in [0,1].")
  if(n.sim <= 0) stop("n should be a positive integer.")
  if(L < 0) stop("L should be positive.")
  if(n.bins4start < 1000) {
    warning("n.bins4start should be large.  Using 1000.")
    n.bins4start <- 1000
  }

  max.nxo <- qpois(1-1e-10, L)*10

  out <- .C("simStahl",
            as.integer(n.sim),
            as.double(nu),
            as.double(p),
            as.double(L),
            nxo = as.integer(rep(0,n.sim)), # number of crossovers
            loc = as.double(rep(0,n.sim*max.nxo)),
            as.integer(max.nxo),
            as.integer(n.bins4start),
            PACKAGE="xoi")

  out <- lapply(as.data.frame(rbind(out$nxo, matrix(out$loc*100, nrow=max.nxo))),
                function(a) {if(a[1]==0) return(numeric(0)); a[(1:a[1])+1] })
  attr(out, "L") <- L*100
  out
}
