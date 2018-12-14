# MEANDIFF = BIVARIATE #################################

#' @export
roundMP.meandiff <- function(x, y = NULL, deltax = NULL, assumptions = NULL, verbose = FALSE) {
    # validation and conversion
    if (is.null(deltax)||(length(deltax)>1)) stop("deltax must receive an integer value")
    if (MP.rowLengths(x) !=1) stop("x not a vector or a one-column matrix/data.frame")
    if (MP.rowLengths(x)+MP.rowLengths(y) !=2) stop("input x and y must be vectors or both one-column matrix/data.frame")
    x          = MP.flatten(x)
    y          = MP.flatten(y)

    # statistic computations
    dmn        = mean(x)-mean(y)
    pr         = 2 * deltax

    # precision computations
    pr         = pr + pr/10000 # avoid rounding errors
    rd         = round(dmn, -log10(pr)+0.5)
    assumptext = "assumption-free"

    # output results
    if (verbose) MP.showVerbose("meandiff", dmn, deltax, pr, rd, assumptext)
    return(setNames( c(dmn,deltax,pr,rd),
        c("meandiff","delta_x","precision","rounded mean")
    ))
}

