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
    # a. Extrinsinc precision
    prEP       = sqrt(2 /(n-1)) * var
    rdEP       = round(var, -log10(prEP)+0.5)
    # b. Worst-case intrinsinc precision
    if (assumptions) {
        prWC         = (2*n/(n-1)) * deltax * sqrt(2/pi) * sqrt(var)
        assumptext = "based on the normality assumption"
    } else {
        prWC         = (2*n/(n-1)) * deltax * MP.absoluteCentralMoment(x)
        assumptext = "assumption-free"
    }
    rdWC           = round(var, -log10(prWC * 1.0001 )+0.5)
    # c. Best-case instrinsinc precision
    prBC       = 2 / sqrt(n-1) * deltax
    rdBC       = round(var, -log10(prBC)+0.5)
    # d. Middle-ground intrinsinc precision
    prMG       = (prWC + prBC)/2
    rdMG       = round(var, -log10(prMG)+0.5)

    # output results
    if (verbose) MP.showVerbose("var", var, deltax, prEP, rdEP, prWC, rdWC, prBC, rdBC, prMG, rdMG, assumptext)
    return(setNames( c(var, rdEP, rdWC, rdBC, rdMG),
        c("var","EXrounded", "WCrounded", "BCrounded","MGrounded") ) 
    )
}

