# MEANDIFF = BIVARIATE #################################

#' @export
roundMP.meandiff <- function(x, y = NULL, deltax = NULL, assumptions = NULL, verbose = FALSE) {
    # validation and conversion
    if (is.null(deltax)||(length(deltax)>1)) stop("deltax must receive an integer value")
    if (MP.rowLengths(x) !=1) stop("x not a vector or a one-column matrix/data.frame")
    if (MP.rowLengths(x)+MP.rowLengths(y) !=2) stop("input x and y must be vectors or both one-column matrix/data.frame")
    x          <- MP.flatten(x)
    y          <- MP.flatten(y)

    # statistic computations
    dmn        <- mean(x)-mean(y)
    ns         <- c(length(x), length(y))
    nh         <- 1/mean(1/ns)
    sds        <- c(sd(x), sd(y))
    sdp        <- sqrt(sum((ns - 1) * sds^2)/(sum(ns) - length(ns))) 

    # precision computations
    # a. Extrinsinc precision
    prEP       <- sqrt(2) * sdp / sqrt(nh) * 1.0001 # avoids rounding error
    rdEP       <- round(dmn, -log10(prEP)+0.5)
    # b. Worst-case intrinsinc precision
    prWC       <- 2 * deltax * 1.0001 # avoids rounding error
    rdWC       <- round(dmn, -log10(prWC)+0.5)
    assumptext <- "assumption-free"
    # c. Best-case instrinsinc precision
    prBC       <- sqrt(2) * deltax / sqrt(nh)
    rdBC       <- round(dmn, -log10(prBC)+0.5)
    # d. Middle-ground intrinsinc precision
    prMG       <- (prWC/2 + prBC/2)/2
    rdMG       <- round(dmn, -log10(prMG)+0.5)

    # output results
    if (verbose) MP.showVerbose("meandiff", dmn, deltax, prEP, rdEP, prWC, rdWC, prBC, rdBC, prMG, rdMG, assumptext)
    return(setNames( c(dmn, rdEP, rdWC, rdBC, rdMG),
        c("meandiff","EXrounded", "WCrounded", "BCrounded","MGrounded") ) 
    )
}

