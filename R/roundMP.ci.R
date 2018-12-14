#' @export

# CI = UNIVARIATE ######################################
roundMP.ci <- function(x, deltax = NULL, gamma = 0.95, assumptions = TRUE, verbose = FALSE) {
    # validation and conversion
    if (is.null(deltax)||(length(deltax)>1)) stop("deltax must receive an integer value")
    if (MP.rowLengths(x) !=1) stop("input not a vector or a one-column matrix/data.frame")
    if ((gamma<0)||(gamma>1)) stop("the coverage factor gamma must be between 0 and 1")
    x              = MP.flatten(x)
    # statistic computations
    mn             = mean(x)
    sd             = sd(x)
    n              = length(x)
    se             = sd / sqrt(n)
    tc             = qt(1/2 + gamma/2, n-1)
    ci             = c(mn-tc*se, mn+tc*se)
    # precision computations
    if (assumptions) {
        pr         = deltax * (1 + tc * (n/(n-1)) * sqrt(2/pi) / sqrt(n) )
        assumptext = "based on the normality assumption"
    } else {
        pr         = deltax * (1 + tc * (1/(n-1)) * MP.absoluteCentralMoment(x) / se )
        assumptext = "assumption-free"
    }
    pr             = pr + pr/10000 # avoid rounding errors
    rd             = round(ci, -log10(pr)+0.5)
    # output results
    if (verbose) MP.showVerbose("ci", ci, deltax, pr, rd, assumptext)
    return(setNames( c(ci,deltax,pr,rd),
        c("ci.low","ci.up","delta_x","precision","rounded ci.low","rounded ci.up")
    ))
}

