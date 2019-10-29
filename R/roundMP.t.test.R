# T.TEST = UNIVARIATE or BIVARIATE #####################

#' @export
roundMP.t.test <- function(
                        deltax = NULL, mu0 = NULL, assumptions = TRUE,
                        verbose = FALSE, fromStatistics = NULL, 
                        fromData = NULL, fromObject = NULL
) {
    # validation and conversion
    if (is.null(deltax)||(length(deltax)>1)) stop("deltax must be a single real value")
    if ((is.null(fromData))&&(is.null(fromStatistics))&&(is.null(fromObject))) 
      stop("you must use fromStatistics or fromData or fromObject")

    if (!is.null(fromData)) {
        dta <- MP.getData(fromData, "1or2")
        if(MP.rowLengths(dta) == 1) {
            if(is.null(mu0)) stop ("for one-sample t-test, you must provide a value to mu0")
            ngrp <- 1
            x1   <- dta[[1]]
            y1   <- mu0
            y2   <- NULL
        } else {
            ngrp     <- 2
            x1       <- dta[[1]]
            y1 <- y2 <- dta[[2]]
        }
        dmn            <- mean(x1)-mean(y1)
        args           <- list(x1,y2)
        args[sapply(args, is.null)] <- NULL # removing nulls
        sds            <- unlist(lapply(args, sd))
        ns             <- unlist(lapply(args, length))
    } else if (!is.null(fromStatistics)) {
        if (is.null(mu0)) {
            sts  <- MP.vfyStat(fromStatistics, c("mean1","sd1","n1","mean2","sd2","n2"))
            ngrp <- 2
            dmn  <- sts[["mean1"]]-sts[["mean2"]]
            sds  <- c(sts[["sd1"]],sts[["sd2"]])
            ns   <- c(sts[["n1"]],sts[["n2"]])
        } else {
            sts  <- MP.vfyStat(fromStatistics, c("mean","sd","n"))
            ngrp <- 1
            dmn  <- sts[["mean"]]-mu0
            sds  <- sts[["sd"]]
            ns   <- sts[["n"]]        
        }
    } else {
        if (class(fromObject)=="htest") {
            if (fromObject$method == " Two Sample t-test") {
                ngrp <- 2
                dmn  <- fromObject$estimate[1] - fromObject$estimate[2]
                sds  <- c(fromObject$sds[1], fromObject$sds[2])
                ns   <- c(fromObject$ns[1], fromObject$ns[2])
            } else if (fromObject$method == "One Sample t-test") {
                ngrp <- 1
                dmn  <- fromObject$estimate[1] - mu0
                sds  <- fromObject$sds[1]
                ns   <- fromObject$ns[1]
            
            } else 
                stop("Not a regular t.test object (no Welch test accepted)")
        } else 
            stop("Not an adequate object passed to fromObject")
    }
    
    # additional statistic computations
    sdp            <- sqrt(sum((ns - 1) * sds^2)/(sum(ns) - length(ns))) 
    nh             <- 1/mean(1/ns)
    ttest          <- dmn /(sqrt(ngrp) * sdp / sqrt(nh) ) 
    nu             <- sum(ns) - length(ns)
    J              <- if (nu>300) 1 else gamma(nu/2) / (sqrt(nu/2) * gamma((nu-1)/2) )
    
    # precision computations
    # a. Extrinsinc precision
    prEP           <- sqrt(nu/(nu-2)* (1+ ttest^2) - ttest^2/J^2)
    rdEP           <- round(ttest, -log10(prEP * 1.0001 )+0.5)
    # b. Systematic (worst-case) intrinsinc precision
    if (assumptions) {
        prWC       <- sqrt(ngrp) * sqrt(nh) * deltax / sdp
        assumptext <- "based on the normality assumption, an absence of effect and the homogeneity of variance if two groups"
    } else {
        d          <- abs(dmn)/sdp # cohen's d
        #dta        <- list(x,y)
        fsum       <- 1:ngrp
        for (i in 1:ngrp) 
            fsum[i]<- sum(abs(1/ns[i] -1/(sum(ns)-ngrp) * (unlist(dta[i])-mean(unlist(dta[i])))/sd(unlist(dta[i]))*d))
        tsum       <- sum(fsum) 
        prWC       <- tsum * deltax / (sqrt(ngrp) * sdp / sqrt(nh) )
        assumptext <- "assumption-free"
    }
    rdWC           <- round(ttest, -log10(prWC * 1.0001 )+0.5)
    # c. Non-systematic (best-case instrinsinc precision
    prBC           <- deltax / sdp
    rdBC           <- round(ttest, -log10(prBC * 1.0001 )+0.5)

    # output results
    if (verbose) MP.showVerbose("t.test", ttest, deltax, prEP, rdEP, prWC, rdWC, prBC, rdBC, assumptext)
    res <- setNames( c(ttest, rdEP, rdWC, rdBC),
        c("machine.precision","extrinsic","systematic","non.systematic") ) 
    return(as.data.frame(t(res))[getOption("roundMP.selectedScenario")])
}

