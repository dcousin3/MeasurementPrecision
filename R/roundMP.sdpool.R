# SDpool = MULTIVARIATE ####################################

#' @export
roundMP.sdpool <- function(..., deltax = NULL, assumptions = TRUE, verbose = FALSE) {
    # validation and conversion
    if (is.null(deltax)||(length(deltax)>1)) stop("deltax must receive an integer value")
    args          <- list( ... )
    args[sapply(args, is.null)] <- NULL # removing nulls
    if (!length(args)) stop("sdpool needs at least one vector argument")
    args          = lapply(args, as.matrix)
    # unpack all the columns into a single list of columns
    args          = c(lapply( args,
        function(part) as.vector(split(part, col(part)))
    ))
    args          = unlist(args,recursive = F)

    # statistic computations
    sds            = unlist(lapply(args, sd))
    ns             = unlist(lapply(args, length))
    sdp            = sqrt(sum((ns - 1) * sds^2)/(sum(ns) - length(ns))) 

    # precision computations
    if (assumptions) {
        pr         = sum(ns)/(sum(ns)-length(ns)) * deltax * sqrt(2/pi)  
        assumptext = "based on the normality assumption and the homogeneity of variance across groups"
    } else {
        fsum       = 1:length(ns)
        for (i in 1:length(ns)) 
            fsum[i]= ns[i] * MP.absoluteCentralMoment(unlist(args[i]))
        tsum       = sum(fsum) 
            pr         = tsum * 1/(sum(ns)-length(ns)) * deltax /sdp
        assumptext = "assumption-free"
    }
    pr             = pr + pr/10000 # avoid rounding errors
    rd             = round(sdp, -log10(pr)+0.5)
    # output results
    if (verbose) MP.showVerbose("sdpool", sdp, deltax, pr, rd, assumptext)
    return(setNames( c(sdp,deltax,pr,rd),
        c("sdpool","delta_x","precision","rounded sdpool")
    ))
}
