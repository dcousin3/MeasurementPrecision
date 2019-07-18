# SDpool = MULTIVARIATE ####################################

#' @export
roundMP.sdpool <- function(deltax = NULL, assumptions = TRUE, verbose = FALSE, fromStatistics = NULL, fromData = NULL) {
    # validation and conversion
    if (is.null(deltax)||(length(deltax)>1)) stop("deltax must be a single real value")
    if ((is.null(fromData))&&(is.null(fromStatistics))) stop("you must use fromStatistics or fromData")

    if (is.null(fromStatistics)) {
        dta <- MP.getData(fromData, "any")
        sds            <- unlist(lapply(dta, sd))
        ns             <- unlist(lapply(dta, length))
    } else {
        sts  <- MP.vfyStat(fromStatistics, c("sd1","n1","sd2","n2"))
        ngrp <- 1
        sds  <- c(sts[["sd1"]],sts[["sd2"]])
        ns   <- c(sts[["n1"]],sts[["n2"]])
    }
    # additional statistic computations
    sdp            <- sqrt(sum((ns - 1) * sds^2)/(sum(ns) - length(ns))) 

    # precision computations
    # a. Extrinsinc precision
    prEP           <- sdp/sqrt(2*(sum(ns)-length(ns)))
    rdEP           <- round(sdp, -log10(prEP)+0.5)
    # b. Worst-case intrinsinc precision
    if (assumptions) {
        prWC       <- sum(ns)/(sum(ns)-length(ns)) * deltax * sqrt(2/pi)  
        assumptext <- "based on the normality assumption and the homogeneity of variance across groups"
    } else {
        fsum       <- 1:length(ns)
        for (i in 1:length(ns)) 
            fsum[i]<- ns[i] * MP.absoluteCentralMoment(dta[[i]])
        tsum       <- sum(fsum) 
        prWC       <- tsum * 1/(sum(ns)-length(ns)) * deltax /sdp
        assumptext <- "assumption-free"
    }
    rdWC           <- round(sdp, -log10(prWC * 1.0001 )+0.5)
    # c. Best-case instrinsinc precision
    prBC           <- deltax / sqrt(sum(ns) - length(ns))
    rdBC           <- round(sdp, -log10(prBC * 1.0001 )+0.5)    

    # output results
    if (verbose) MP.showVerbose("sdpool", sdp, deltax, prEP, rdEP, prWC, rdWC, prBC, rdBC, assumptext)
    return(setNames( c(sdp, rdEP, rdWC, rdBC),
        c("sdpool","EXrounded", "WCrounded", "BCrounded") ) 
    )
}
