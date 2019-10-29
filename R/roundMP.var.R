# VAR = UNIVARIATE #####################################

#' @export
roundMP.var <- function(deltax = NULL, assumptions = TRUE, verbose = FALSE, fromStatistics = NULL, fromData = NULL) {
    # validation and conversion
    if (is.null(deltax)||(length(deltax)>1)) stop("deltax must be a single real value")
    if ((is.null(fromData))&&(is.null(fromStatistics))) stop("you must use fromStatistics or fromData")

    if (is.null(fromStatistics)) {
        dta <- MP.getData(fromData, "1")
        var  <- var(dta)
        n   <- length(dta)
    } else {
        sts <- MP.vfyStat(fromStatistics, c("var","n"))
        var  <- sts[["var"]]
        n   <- sts[["n"]]
    }
    
    # precision computations
    # a. Extrinsinc precision
    prEP           <- sqrt(2 /(n-1)) * var
    rdEP           <- round(var, -log10(prEP)+0.5)
    # b. Systematic (worst-case) intrinsinc precision
    if (assumptions) {
        prWC       <- (2*n/(n-1)) * deltax * sqrt(2/pi) * sqrt(var)
        assumptext <- "based on the normality assumption"
    } else {
        prWC       <- (2*n/(n-1)) * deltax * MP.absoluteCentralMoment(x)
        assumptext <- "assumption-free"
    }
    rdWC           <- round(var, -log10(prWC * 1.0001 )+0.5)
    # c. Non-systematic (best-case) instrinsinc precision
    prBC           <- 2 / sqrt(n-1) * deltax
    rdBC           <- round(var, -log10(prBC)+0.5)

    # output results
    if (verbose) MP.showVerbose("var", var, deltax, prEP, rdEP, prWC, rdWC, prBC, rdBC, assumptext)
    res <- setNames( c(var, rdEP, rdWC, rdBC),
        c("machine.precision","extrinsic","systematic","non.systematic") ) 
    return(as.data.frame(t(res))[getOption("roundMP.selectedScenario")])
}

