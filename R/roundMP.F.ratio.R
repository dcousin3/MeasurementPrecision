## FRATIO = MULTIVARIATE ################################
# the best-case intrinsinc scenario remains to be found
# as well as the standard error of an F ratio

#' @export
roundMP.F.ratio <- function(deltax = NULL, assumptions = TRUE, verbose = FALSE, fromStatistics = NULL, fromData = NULL) {
    # validation and conversion
    if (is.null(deltax)||(length(deltax)>1)) stop("deltax must be a single real value")
    if ((is.null(fromData))&&(is.null(fromStatistics))) stop("you must use fromStatistics or fromData")

    if (is.null(fromStatistics)) {
        dta  <- MP.getData(fromData, "any")
        mns  <- unlist(lapply(dta, mean))
        sds  <- unlist(lapply(dta, sd))
        ns   <- unlist(lapply(dta, length))
        ntilde <- 1/sum(1/ns)*length(ns)
    } else {
        stop ("cannot round F ratio statistic from statistics.")
    }
    # additional statistic computations
    sdp        <- sqrt(sum((ns - 1) * sds^2)/(sum(ns) - length(ns))) 
    GM         <- mean(mns)
    sdm        <- sd(mns)
    ssa        <- sum(ntilde * (mns-GM)^2 )
    ss         <- 1:length(ns)
    for (i in 1:length(ns)) {
      ss[i]    <- sum((dta[[i]]-mns[i])^2)      
    }
    sse        <- sum(ss)
    fratio     <- (ssa / (length(ns)-1)) / (sse/(sum(ns)-length(ns)))
    
    # precision computations
    # a. Extrinsinc precision
    prEP       <- NA
    rdEP       <- NA
    # b. Systematic (worst-case) intrinsinc precision
    fsum       <- 1:length(ns)
    for (i in 1:length(ns)) {
        fsum[i]<- sum(abs(ssa/sse*(dta[[i]]-mns[i]) -(mns[i]-GM) ))
    }
    tsum       <- sum(fsum) 
    pr         <- tsum * 2/(length(ns)-1) * deltax /sdp^2
    assumptext <- "assumption-free"
    prWC       <- pr * 1.00001 # avoid rounding errors
    rdWC       <- round(fratio, -log10(prWC)+0.5)
    # c. Non-systematic (best-case) instrinsinc precision
    prBC       <- NA
    rdBC       <- NA

        
    # output results
    if (verbose) MP.showVerbose("F.ratio", fratio, deltax, prEP, rdEP, prWC, rdWC, prBC, rdBC, assumptext)
    res <- setNames( c(fratio, rdEP, rdWC, rdBC),
        c("machine.precision","extrinsic","systematic","non.systematic") ) 
    return(as.data.frame(t(res))[getOption("roundMP.selectedScenario")])
}

