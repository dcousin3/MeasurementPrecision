# SEmean = UNIVARIATE ###################################

#' @export
roundMP.semean <- function(x, deltax = NULL, assumptions = TRUE, verbose = FALSE) {
    # validation and conversion
    if (is.null(deltax)||(length(deltax)>1)) stop("deltax must receive an integer value")
    if (MP.rowLengths(x) !=1) stop("input not a vector or a one-column matrix/data.frame")
    x              <- MP.flatten(x)
    # statistic computations
    sd             <- sd(x)
    n              <- length(x)
    sem            <- sd / sqrt(n)

    # precision computations
    # a. Extrinsinc precision
    prEP           <- sd/sqrt(2 * n * (n-1) )
    rdEP           <- round(sem, -log10(prEP)+0.5)
    # b. Worst-case intrinsinc precision
    if (assumptions) {
        prWC       <- (n/(n-1)) * deltax * sqrt(2/pi) / sqrt(n)
        assumptext <- "based on the normality assumption"
    } else {
        prWC       <- (1/(n-1)) * deltax * MP.absoluteCentralMoment(x) / sem 
        assumptext <- "assumption-free"
    }
    rdWC           <- round(sem, -log10(prWC * 1.0001 )+0.5)
    # c. Best-case instrinsinc precision
    prBC           <- deltax/sqrt(n * (n-1) )
    rdBC           <- round(sem, -log10(prBC)+0.5)
    # d. Middle-ground intrinsinc precision
    prMG           <- (prWC/2 + prBC/2)/2
    rdMG           <- round(sem, -log10(prMG)+0.5)

    # output results
    if (verbose) MP.showVerbose("semean", sem, deltax, prEP, rdEP, prWC, rdWC, prBC, rdBC, prMG, rdMG, assumptext)
    return(setNames( c(sem, rdEP, rdWC, rdBC, rdMG),
        c("semean","EXrounded", "WCrounded", "BCrounded","MGrounded") ) 
    )
}

