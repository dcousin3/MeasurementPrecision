# T.TEST = UNIVARIATE or BIVARIATE #####################

#' @export
roundMP.t.test <- function(x, y = NULL, mu0 = NULL, deltax = NULL, assumptions = TRUE, verbose = FALSE) {
    # validation and conversion
    if (is.null(deltax)||(length(deltax)>1)) stop("deltax must receive an integer value")
    if (MP.rowLengths(x) !=1) stop("x not a vector or a one-column matrix/data.frame")
    if (!(MP.rowLengths(x)+MP.rowLengths(y) %in% c(1,2)) ) stop("input not containing 1 (single-group) or 2 (2 groups) columns of data")
    if ((MP.rowLengths(x)+MP.rowLengths(y) ==1)&&(is.null(mu0))) stop("for a one-sample t.test, you must provide the parameter mu0")
    x              <- MP.flatten(x)
    y              <- MP.flatten(y)
    if (is.null(y))               {ngrp <- 1; x1 <- x; y1 <- mu0}
    if (is.null(x))               {ngrp <- 1; x1 <- y; y1 <- mu0}
    if (!is.null(x)&&!is.null(y)) {ngrp <- 2; x1 <- x; y1 <- y}

    # statistic computations
    dmn            <- mean(x1)-mean(y1)
    args           <- list(x,y)
    args[sapply(args, is.null)] <- NULL # removing nulls
    sds            <- unlist(lapply(args, sd))
    ns             <- unlist(lapply(args, length))
    sdp            <- sqrt(sum((ns - 1) * sds^2)/(sum(ns) - length(ns))) 
    nh             <- 1/mean(1/ns)
    ttest          <- dmn /(sqrt(ngrp) * sdp / sqrt(nh) ) 
    eta            <- sum(ns) - length(ns)
    J              <- if (eta>300) 1 else gamma(eta/2) / (sqrt(eta/2) * gamma((eta-1)/2) )
    
    # precision computations
    # a. Extrinsinc precision
    prEP           <- sqrt(eta/(eta-2)* (1+ ttest^2) - ttest^2/J^2)
    rdEP           <- round(ttest, -log10(prEP * 1.0001 )+0.5)
    # b. Worst-case intrinsinc precision
    if (assumptions) {
        prWC       <- sqrt(ngrp) * sqrt(nh) * deltax / sdp
        assumptext <- "based on the normality assumption, an absence of effect and the homogeneity of variance if two groups"
    } else {
        d          <- abs(dmn)/sdp # cohen's d
        dta        <- list(x,y)
        fsum       <- 1:ngrp
        for (i in 1:ngrp) 
            fsum[i]<- sum(abs(1/ns[i] -1/(sum(ns)-ngrp) * (unlist(dta[i])-mean(unlist(dta[i])))/sd(unlist(dta[i]))*d))
        tsum       <- sum(fsum) 
        prWC       <- tsum * deltax / (sqrt(ngrp) * sdp / sqrt(nh) )
        assumptext <- "assumption-free"
    }
    rdWC           <- round(ttest, -log10(prWC * 1.0001 )+0.5)
    # c. Best-case instrinsinc precision
    prBC           <- deltax / sdp
    rdBC           <- round(ttest, -log10(prBC * 1.0001 )+0.5)
    # d. Middle-ground intrinsinc precision
    prMG           <- (prWC/2 + prBC/2)/2
    rdMG           <- round(ttest, -log10(prMG)+0.5)

    # output results
    if (verbose) MP.showVerbose("t.test", ttest, deltax, prEP, rdEP, prWC, rdWC, prBC, rdBC, prMG, rdMG, assumptext)
    return(setNames( c(ttest, rdEP, rdWC, rdBC, rdMG),
        c("t.test","EXrounded", "WCrounded", "BCrounded","MGrounded") ) 
    )
}

