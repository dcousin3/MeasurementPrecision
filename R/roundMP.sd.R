# SD = UNIVARIATE ######################################

#' @export
roundMP.sd <- function(deltax = NULL, assumptions = TRUE, verbose = FALSE, fromStatistics = NULL, fromData = NULL) {
    # validation and conversion
    if (is.null(deltax)||(length(deltax)>1)) stop("deltax must be a single real value")
    if ((is.null(fromData))&&(is.null(fromStatistics))) stop("you must use fromStatistics or fromData")

    if (is.null(fromStatistics)) {
        dta <- MP.getData(fromData, "1")
        sd  <- sd(dta)
        n   <- length(dta)
    } else {
        sts <- MP.vfyStat(fromStatistics, c("sd","n"))
        sd  <- sts[["sd"]]
        n   <- sts[["n"]]
    }

    # precision computations
    # a. Extrinsinc precision
    prEP           <- sd/sqrt(2*(n-1))
    rdEP           <- round(sd, -log10(prEP)+0.5)
    # b. Systematic (worst-case) intrinsinc precision
    if (assumptions) {
        prWC       <- (n/(n-1)) * deltax * sqrt(2/pi)
        assumptext <- "based on the normality assumption"
    } else {
        prWC       <- (n/(n-1)) * deltax * MP.absoluteCentralMoment(x)/sd
        assumptext <- "assumption-free"
    }
    rdWC           <- round(sd, -log10(prWC * 1.0001 )+0.5)
    # c. Non-systematic (best-case) instrinsinc precision
    prBC           <- deltax/sqrt(n-1)
    rdBC           <- round(sd, -log10(prBC)+0.5)

    # output results
    if (verbose) MP.showVerbose("sd", sd, deltax, prEP, rdEP, prWC, rdWC, prBC, rdBC, assumptext)
    res <- setNames( c(sd, rdEP, rdWC, rdBC),
        c("machine.precision","extrinsic","systematic","non.systematic") ) 
    return(as.data.frame(t(res)))
}

