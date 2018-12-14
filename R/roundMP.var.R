# VAR = UNIVARIATE #####################################

#' @export
roundMP.var <- function(x, deltax = NULL, assumptions = TRUE, verbose = FALSE) {
    # validation and conversion
    if (is.null(deltax)||(length(deltax)>1)) stop("deltax must receive an integer value")
    if (MP.rowLengths(x) !=1) stop("input not a vector or a one-column matrix/data.frame")
    x = MP.flatten(x)
    # statistic computations
    var            = var(x)
    n              = length(x)
    # precision computations
    if (assumptions) {
        pr         = (2*n/(n-1)) * deltax * sqrt(2/pi) * sqrt(var)
        assumptext = "based on the normality assumption"
    } else {
        pr         = (2*n/(n-1)) * deltax * MP.absoluteCentralMoment(x)
        assumptext = "assumption-free"
    }
    pr             = pr + pr/10000 # avoid rounding errors
    rd             = round(var, -log10(pr)+0.5)
    # output results
    if (verbose) MP.showVerbose("var", var, deltax, pr, rd, assumptext)
    return(setNames( c(var,deltax,pr,rd),
        c("var","delta_x","precision","rounded var")
    ))
}

