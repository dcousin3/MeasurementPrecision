# SD = UNIVARIATE ######################################

#' @export
roundMP.sd <- function(x, deltax = NULL, assumptions = TRUE, verbose = FALSE) {
    # validation and conversion
    if (is.null(deltax)||(length(deltax)>1)) stop("deltax must receive an integer value")
    if (MP.rowLengths(x) !=1) stop("input not a vector or a one-column matrix/data.frame")
    x              = MP.flatten(x)
    # statistic computations
    sd             = sd(x)
    n              = length(x)
    # precision computations
    if (assumptions) {
        pr         = (n/(n-1)) * deltax * sqrt(2/pi)
        assumptext = "based on the normality assumption"
    } else {
        pr         = (n/(n-1)) * deltax * MP.absoluteCentralMoment(x)/sd
        assumptext = "assumption-free"
    }
    pr             = pr + pr/10000 # avoid rounding errors
    rd             = round(sd, -log10(pr)+0.5)
    # output results
    if (verbose) MP.showVerbose("sd", sd, deltax, pr, rd, assumptext)
    return(setNames(
        c(sd,deltax,pr,rd),
        c("sd","delta_x","precision","rounded sd")
    ))
}

