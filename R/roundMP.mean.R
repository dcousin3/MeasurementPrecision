# MEAN = UNIVARIATE ####################################

#' @export
roundMP.mean <- function(x, deltax = NULL, assumptions = NULL, verbose = FALSE) {
    # validation and conversion
    if (is.null(deltax)||(length(deltax)>1)) stop("deltax must receive an integer value")
    if (MP.rowLengths(x) !=1) stop("input not a vector or a one-column matrix/data.frame")
    x          <- MP.flatten(x)
    # statistic computations
    mn         <- mean(x)

    # precision computations
    # a. Extrinsinc precision
    prEP       <- sd(x)/sqrt(length(x))
    rdEP       <- round(mn, -log10(prEP)+0.5)
    # b. Worst-case intrinsinc precision
    prWC       <- deltax * 1.0001 # avoids rounding error
    rdWC       <- round(mn, -log10(prWC)+0.5)
    assumptext <- "assumption-free"
    # c. Best-case instrinsinc precision
    prBC       <- deltax / sqrt(length(x)) * 1.0001 # avoids rounding error
    rdBC       <- round(mn, -log10(prBC)+0.5)
    # d. Middle-ground intrinsinc precision
    prMG       <- (prWC/2 + prBC/2)/2
    rdMG       <- round(mn, -log10(prMG)+0.5)
    
    # output results
    if (verbose) MP.showVerbose("mean", mn, deltax, prEP, rdEP, prWC, rdWC, prBC, rdBC, prMG, rdMG, assumptext)
    return(setNames( c(mn, rdEP, rdWC, rdBC, rdMG),
        c("mean","EXrounded", "WCrounded", "BCrounded","MGrounded") ) 
    )
}
