# COHEND = UNIVARIATE or BIVARIATE #####################

#' @export
roundMP.cohen.d <- function(deltax = NULL, mu0 = NULL, assumptions = TRUE, verbose = FALSE, fromStatistics = NULL, fromData = NULL) {
    # validation and conversion
    if (is.null(deltax)||(length(deltax)>1)) stop("deltax must be a single real value")
    if ((is.null(fromData))&&(is.null(fromStatistics))) stop("you must use fromStatistics or fromData")

    if (is.null(fromStatistics)) {
        dta <- MP.getData(fromData, "1or2")
        if(MP.rowLengths(dta) == 1) {
            if(is.null(mu0)) stop ("for one-sample d_1, you must provide a value to mu0")
            ngrp <- 1
            x1   <- dta
            y1   <- mu0
            y2   <- NULL
        } else {
            ngrp     <- 2
            x1       <- dta[,1]
            y1 <- y2 <- dta[,2]
        }
        dmn            <- mean(x1)-mean(y1)
        args           <- list(x1,y2)
        args[sapply(args, is.null)] <- NULL # removing nulls
        sds            <- unlist(lapply(args, sd))
        ns             <- unlist(lapply(args, length))
    } else {
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
    }
    # additional statistic computations
    sdp            <- sqrt(sum((ns - 1) * sds^2)/(sum(ns) - length(ns))) 
    nh             <- 1/mean(1/ns)
    d              <- abs(dmn)/sdp # cohen's d
    nu            <- sum(ns) - length(ns)
    J              <- if (nu>300) 1 else gamma(nu/2) / (sqrt(nu/2) * gamma((nu-1)/2) )
    g              <-  d * J  # this is unbiased Cohen's d (Goulet-Pelletier & Cousineau, 2018)

    # precision computations
    # a. Extrinsinc precision
    prEP           <- sqrt(nu/(nu-2)*2/nh * (1+ g^2 *nh/2)- g^2/J^2)
    rdEP           <- round(g, -log10(prEP * 1.0001)+0.5)
    # b. Worst-case intrinsinc precision
    if (assumptions) {
        prWC       <- ngrp * deltax / sdp
        assumptext <- "based on the normality assumption, an absence of effect and the homogeneity of variance if two groups"
    } else {
        dta        <- list(x,y)
        fsum       <- 1:ngrp
        for (i in 1:ngrp) 
            fsum[i]<- sum(abs(1/ns[i] -1/(sum(ns)-ngrp) * (unlist(dta[i])-mean(unlist(dta[i])))/sd(unlist(dta[i]))*d))
        tsum       <- sum(fsum) 
        prWC       <- tsum * deltax / sdp
        assumptext <- "assumption-free"
    }
    rdWC           <- round(g, -log10(prWC * 1.0001 )+0.5)
    # c. Best-case instrinsinc precision
    prBC           <- sqrt(ngrp)/sqrt(nh) *deltax / sdp
    rdBC           <- round(g, -log10(prBC * 1.0001)+0.5)

    # output results
    if (verbose) MP.showVerbose("cohen.d(unbs'd)", g, deltax, prEP, rdEP, prWC, rdWC, prBC, rdBC, assumptext)
    return(setNames( c(g, rdEP, rdWC, rdBC),
        c("g","EXrounded", "WCrounded", "BCrounded") ) 
    )
}

