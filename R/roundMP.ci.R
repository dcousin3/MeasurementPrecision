#' @export

# CI = UNIVARIATE ######################################
roundMP.ci <- function(x, deltax = NULL, gamma = 0.95, assumptions = TRUE, verbose = FALSE) {
    # validation and conversion
    if (is.null(deltax)||(length(deltax)>1)) stop("deltax must receive an integer value")
    if (MP.rowLengths(x) !=1) stop("input not a vector or a one-column matrix/data.frame")
    if ((gamma<0)||(gamma>1)) stop("the coverage factor gamma must be between 0 and 1")
    x              = MP.flatten(x)
    # statistic computations
    mn             = mean(x)
    sd             = sd(x)
    n              = length(x)
    se             = sd / sqrt(n)
    tc             = qt(1/2 + gamma/2, n-1)
    ci             = c(mn-tc*se, mn+tc*se)
    
    # precision computations
    # a. Extrinsinc precision
    prEP       = sd/sqrt(n) * sqrt(1+ tc^2 / (n-1) )
    rdEP       = round(ci, -log10(prEP)+0.5)
    # b. Worst-case intrinsinc precision
    if (assumptions) {
        prWC       = deltax * (1 + tc * (n/(n-1)) * sqrt(2/pi) / sqrt(n) )
        assumptext = "based on the normality assumption"
    } else {
        prWC       = deltax * (1 + tc * (1/(n-1)) * MP.absoluteCentralMoment(x) / se )
        assumptext = "assumption-free"
    }
    rdWC           = round(ci, -log10(prWC * 1.0001 )+0.5)
    # c. Best-case instrinsinc precision
    prBC       = deltax/sqrt(n) * sqrt(1+ tc^2 / (n-1) )
    rdBC       = round(ci, -log10(prBC)+0.5)
    # d. Middle-ground intrinsinc precision
    prMG       = (prWC + prBC)/2
    rdMG       = round(ci, -log10(prMG)+0.5)

    # output results
    if (verbose) MP.showVerbose("ci", c(ci[1]," ",ci[2]), deltax, prEP, c(rdEP[1]," ",rdEP[2]), prWC, c(rdWC[1]," ",rdWC[2]), prBC, c(rdBC[1]," ",rdBC[2]), prMG, c(rdMG[1]," ",rdMG[2]), assumptext)
    return(setNames( c(ci, rdEP, rdWC, rdBC, rdMG),
        c("ci-low","ci-high","EXrounded-lo", "EXrounded-hi", "WCrounded-lo","WCrounded-hi","BCrounded-lo","BCrounded-hi","MGrounded-lo","MGrounded-hi") ) 
    )
}

