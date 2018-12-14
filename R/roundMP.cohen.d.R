# COHEND = UNIVARIATE or BIVARIATE #####################

#' @export
roundMP.cohen.d <- function(x, y = NULL, mu0 = NULL, deltax = NULL, assumptions = TRUE, verbose = FALSE) {
    # validation and conversion
    if (is.null(deltax)||(length(deltax)>1)) stop("deltax must receive an integer value")
    if (MP.rowLengths(x) !=1) stop("x not a vector or a one-column matrix/data.frame")
    if (!(MP.rowLengths(x)+MP.rowLengths(y) %in% c(1,2)) ) stop("input not containing 1 (single-group) or 2 (2 groups) columns of data")
    if ((MP.rowLengths(x)+MP.rowLengths(y) ==1)&&(is.null(mu0))) stop("for a one-sample t.test, you must provide the parameter mu0")
    x          = MP.flatten(x)
    y          = MP.flatten(y)
    if (is.null(y))               {ngrp = 1; x1 = x; y1 = mu0 }
    if (is.null(x))               {ngrp = 1; x1 = y; y1 = mu0 }
    if (!is.null(x)&&!is.null(y)) {ngrp = 2; x1 = x; y1 = y }

    # statistic computations
    dmn        = mean(x1)-mean(y1)
    args       <- list(x,y)
    args[sapply(args, is.null)] <- NULL # removing nulls
    sds        = unlist(lapply(args, sd))
    ns         = unlist(lapply(args, length))
    sdp        = sqrt(sum((ns - 1) * sds^2)/(sum(ns) - length(ns))) 
    nh         = 1/mean(1/ns)
    d          = abs(dmn)/sdp # cohen's d

    # precision computations
    if (assumptions) {
        pr         = ngrp * deltax / sdp
        assumptext = "based on the normality assumption, an absence of effect and the homogeneity of variance if two groups"
    } else {
        dta        = list(x,y)
        fsum       = 1:ngrp
        for (i in 1:ngrp) 
            fsum[i]= sum(abs(1/ns[i] -1/(sum(ns)-ngrp) * (unlist(dta[i])-mean(unlist(dta[i])))/sd(unlist(dta[i]))*d))
        tsum       = sum(fsum) 
        pr         = tsum * deltax / sdp
        assumptext = "assumption-free"
    }
    pr             = pr + pr/10000 # avoid rounding errors
    rd             = round(d, -log10(pr)+0.5)
    # output results
    if (verbose) MP.showVerbose("Cohen d", d, deltax, pr, rd, assumptext)
    return(setNames( c(d,deltax,pr,rd),
        c("Cohen d","delta_x","precision","rounded Cohen d")
    ))
}

