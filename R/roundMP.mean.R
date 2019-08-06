# MEAN = UNIVARIATE ####################################

#' @export
roundMP.mean <- function(deltax = NULL, assumptions = NULL, verbose = FALSE, fromStatistics = NULL, fromData = NULL) {
    # validation and conversion
    if (is.null(deltax)||(length(deltax)>1)) stop("deltax must be a single real value")
    if ((is.null(fromData))&&(is.null(fromStatistics))) stop("you must use fromStatistics or fromData")

    if (is.null(fromStatistics)) {
        dta <- MP.getData(fromData, "1")
        mn  <- mean(dta)
        sd  <- sd(dta)
        n   <- length(dta)
    } else {
        sts <- MP.vfyStat(fromStatistics, c("mean","sd","n"))
        mn  <- sts[["mean"]]
        sd  <- sts[["sd"]]
        n   <- sts[["n"]]
    }
    
    # precision computations
    # a. Extrinsinc precision
    prEP       <- sd/sqrt(n)
    rdEP       <- round(mn, -log10(prEP)+0.5)
    # b. Systematic (worst-case) intrinsinc precision
    prWC       <- deltax * 1.0001 # avoids rounding error
    rdWC       <- round(mn, -log10(prWC)+0.5)
    assumptext <- "assumption-free"
    # c. Non-systematic (best-case) instrinsinc precision
    prBC       <- deltax / sqrt(n) * 1.0001 # avoids rounding error
    rdBC       <- round(mn, -log10(prBC)+0.5)
    
    # output results
    if (verbose) MP.showVerbose("mean", mn, deltax, prEP, rdEP, prWC, rdWC, prBC, rdBC, assumptext)
    res <- setNames( c(mn, rdEP, rdWC, rdBC),
        c("machine.precision","extrinsic","systematic","non.systematic") ) 
    return(as.data.frame(t(res)))
}
