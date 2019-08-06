#' @export

# CI = UNIVARIATE ######################################
roundMP.ci <- function(deltax = NULL, gamma = 0.95, assumptions = TRUE, verbose = FALSE, fromStatistics = NULL, fromData = NULL) {
    # validation and conversion
    if (is.null(deltax)||(length(deltax)>1)) stop("deltax must be a single real value")
    if ((is.null(fromData))&&(is.null(fromStatistics))) stop("you must use fromStatistics or fromData")

    if (is.null(fromStatistics)) {
        dta <- MP.getData(fromData, "1")
        mn             <- mean(dta)
        sd             <- sd(dta)
        n              <- length(dta)
    } else {
        sts <- MP.vfyStat(fromStatistics, c("mean","sd","n"))
        mn  <- sts[["mean"]]
        sd  <- sts[["sd"]]
        n   <- sts[["n"]]
    }
    # additional statistic computations
    se             <- sd / sqrt(n)
    tc             <- qt(1/2 + gamma/2, n-1)
    ci             <- c(mn-tc*se, mn+tc*se)
    
    # precision computations
    # a. Extrinsinc precision
    prEP           <- sd/sqrt(n) * sqrt(1+ tc^2 / (n-1) )
    rdEP           <- round(ci, -log10(prEP)+0.5)
    # b. Systematic (worst-case) intrinsinc precision
    if (assumptions) {
        prWC       <- deltax * (1 + tc * (n/(n-1)) * sqrt(2/pi) / sqrt(n) )
        assumptext <- "based on the normality assumption"
    } else {
        prWC       <- deltax * (1 + tc * (1/(n-1)) * MP.absoluteCentralMoment(x) / se )
        assumptext <- "assumption-free"
    }
    rdWC           <- round(ci, -log10(prWC * 1.0001 )+0.5)
    # c. Non-systematic (best-case) instrinsinc precision
    prBC           <- deltax/sqrt(n) * sqrt(1+ tc^2 / (n-1) )
    rdBC           <- round(ci, -log10(prBC)+0.5)

    # output results
    if (verbose) MP.showVerbose("ci", c(ci[1]," ",ci[2]), deltax, prEP, c(rdEP[1]," ",rdEP[2]), prWC, c(rdWC[1]," ",rdWC[2]), prBC, c(rdBC[1]," ",rdBC[2]), assumptext)
    res <- setNames( c(ci, rdEP, rdWC, rdBC),
        c("machine.precision.low","machine.precision.high","extrinsic.low","extrinsic.high",
          "systematic.low","systematic.high","non.systematic.low","non.systematic.high") ) 
    return(as.data.frame(t(res)))
}

