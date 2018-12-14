# MEAN = UNIVARIATE ####################################

#' @export
roundMP.mean <- function(x, deltax = NULL, assumptions = NULL, verbose = FALSE) {
    # validation and conversion
    if (is.null(deltax)||(length(deltax)>1)) stop("deltax must receive an integer value")
    if (MP.rowLengths(x) !=1) stop("input not a vector or a one-column matrix/data.frame")
    x          = MP.flatten(x)
    # statistic computations
    mn         = mean(x)
    # precision computations
    pr         = deltax 
    pr         = pr + pr/10000 # avoid rounding errors
    rd         = round(mn, -log10(pr)+0.5)
    assumptext = "assumption-free"
    # output results
    if (verbose) MP.showVerbose("mean", mn, deltax, pr, rd, assumptext)
    return(setNames( c(mn,deltax,pr,rd),
        c("mean","delta_x","precision","rounded mean") ) 
    )
}
