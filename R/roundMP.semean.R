# SEmean = UNIVARIATE ###################################

#' @export
roundMP.semean <- function(deltax = NULL, assumptions = TRUE, verbose = FALSE, fromStatistics = NULL, fromData = NULL) {
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
    sem            <- sd / sqrt(n)

    # precision computations
    # a. Extrinsinc precision
    prEP           <- sd/sqrt(2 * n * (n-1) )
    rdEP           <- round(sem, -log10(prEP)+0.5)
    # b. Systematic (worst-case) intrinsinc precision
    if (assumptions) {
        prWC       <- (n/(n-1)) * deltax * sqrt(2/pi) / sqrt(n)
        assumptext <- "based on the normality assumption"
    } else {
        prWC       <- (1/(n-1)) * deltax * MP.absoluteCentralMoment(x) / sem 
        assumptext <- "assumption-free"
    }
    rdWC           <- round(sem, -log10(prWC * 1.0001 )+0.5)
    # c. Non-systematic (best-case) instrinsinc precision
    prBC           <- deltax/sqrt(n * (n-1) )
    rdBC           <- round(sem, -log10(prBC)+0.5)

    # output results
    if (verbose) MP.showVerbose("semean", sem, deltax, prEP, rdEP, prWC, rdWC, prBC, rdBC, assumptext)
    res <- setNames( c(sem, rdEP, rdWC, rdBC),
        c("machine.precision","extrinsic","systematic","non.systematic") ) 
    return(as.data.frame(t(res))[getOption("roundMP.selectedScenario")])
}

