######################################################################
# print the installed version of R/xoi
######################################################################
xoiversion <-
function()
{
  u <- strsplit(library(help=xoi)[[3]][[1]][2]," ")[[1]]
  u[length(u)]
}
