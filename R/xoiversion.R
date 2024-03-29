## xoiversion.R

#' Installed version of R/xoi
#'
#' Print the version number of the currently installed version of R/xoi.
#'
#'
#' @return A character string with the version number of the currently
#' installed version of R/xoi.
#' @author Karl W Broman, \email{broman@@wisc.edu}
#' @keywords print
#' @examples
#' xoiversion()
#'
#' @export xoiversion
xoiversion <-
    function()
{
    version <- unlist(utils::packageVersion("xoi"))

    if(length(version) == 3) {
        # make it like #.#-#
        return( paste(c(version, ".", "-")[c(1,4,2,5,3)], collapse="") )
    }

    paste(version, collapse=".")
}
