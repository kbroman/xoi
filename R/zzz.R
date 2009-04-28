######################################################################
#
# zzz.R
#
# copyright (c) 2006, Karl W Broman
# written Nov, 2006
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/xoi package
#
# .First.lib is run when the package is loaded with library(xoi)
#
######################################################################

.First.lib <- function(lib, pkg) library.dynam("xoi", pkg, lib)

# end of zzz.R
