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
    # a. Extrinsinc precision
    prEP       = sd/sqrt(2*(n-1))
    rdEP       = round(sd, -log10(prEP)+0.5)
    # b. Worst-case intrinsinc precision
    if (assumptions) {
        prWC       = (n/(n-1)) * deltax * sqrt(2/pi)
        assumptext = "based on the normality assumption"
    } else {
        prWC       = (n/(n-1)) * deltax * MP.absoluteCentralMoment(x)/sd
        assumptext = "assumption-free"
    }
    rdWC           = round(sd, -log10(prWC * 1.0001 )+0.5)
    # c. Best-case instrinsinc precision
    prBC       = deltax/sqrt(n-1)
    rdBC       = round(sd, -log10(prBC)+0.5)
    # d. Middle-ground intrinsinc precision
    prMG       = (prWC + prBC)/2
    rdMG       = round(sd, -log10(prMG)+0.5)

    # output results
    if (verbose) MP.showVerbose("sd", sd, deltax, prEP, rdEP, prWC, rdWC, prBC, rdBC, prMG, rdMG, assumptext)
    return(setNames( c(sd, rdEP, rdWC, rdBC, rdMG),
        c("sd","EXrounded", "WCrounded", "BCrounded","MGrounded") ) 
    )
}

