# MEANDIFF = BIVARIATE #################################

#' @export
roundMP.meandiff <- function(deltax = NULL, assumptions = NULL, verbose = FALSE, fromStatistics = NULL, fromData = NULL) {
    # validation and conversion
    if (is.null(fromStatistics)) {
        dta <- MP.getData(fromData, "2")
        x1  <- dta[[1]]
        y1  <- dta[[2]]
        dmn            <- mean(x1)-mean(y1)
        args           <- list(x1,y1)
        sds            <- unlist(lapply(args, sd))
        ns             <- unlist(lapply(args, length))
    } else {
        sts  <- MP.vfyStat(fromStatistics, c("mean1","sd1","n1","mean2","sd2","n2"))
        dmn  <- sts[["mean1"]]-sts[["mean2"]]
        sds  <- c(sts[["sd1"]],sts[["sd2"]])
        ns   <- c(sts[["n1"]],sts[["n2"]])
    }
    
    # statistic computations
    nh         <- 1/mean(1/ns)
    sdp        <- sqrt(sum((ns - 1) * sds^2)/(sum(ns) - length(ns))) 

    # precision computations
    # a. Extrinsinc precision
    prEP       <- sqrt(2) * sdp / sqrt(nh) * 1.0001 # avoids rounding error
    rdEP       <- round(dmn, -log10(prEP)+0.5)
    # b. Systematic (worst-case) intrinsinc precision
    prWC       <- 2 * deltax * 1.0001 # avoids rounding error
    rdWC       <- round(dmn, -log10(prWC)+0.5)
    assumptext <- "assumption-free"
    # c. Non-systematic (best-case) instrinsinc precision
    prBC       <- sqrt(2) * deltax / sqrt(nh)
    rdBC       <- round(dmn, -log10(prBC)+0.5)

    # output results
    if (verbose) MP.showVerbose("meandiff", dmn, deltax, prEP, rdEP, prWC, rdWC, prBC, rdBC, assumptext)
    res <- setNames( c(dmn, rdEP, rdWC, rdBC),
        c("machine.precision","extrinsic","systematic","non.systematic") ) 
    return(as.data.frame(t(res)))
}

