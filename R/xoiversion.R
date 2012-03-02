######################################################################
# print the installed version of R/xoi
######################################################################
xoiversion <-
function()
{
  version <- unlist(packageVersion("xoi"))

  # make it like #.#-#
  paste(c(version,".","-")[c(1,4,2,5,3)], collapse="")
}
