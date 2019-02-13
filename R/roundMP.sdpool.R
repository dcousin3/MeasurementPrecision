# SDpool = MULTIVARIATE ####################################

#' @export
roundMP.sdpool <- function(..., deltax = NULL, assumptions = TRUE, verbose = FALSE) {
    # validation and conversion
    if (is.null(deltax)||(length(deltax)>1)) stop("deltax must receive an integer value")
    args          <- list( ... )
    args[sapply(args, is.null)] <- NULL # removing nulls
    if (!length(args)) stop("sdpool needs at least one vector argument")
    args          <- lapply(args, as.matrix)
    # unpack all the columns into a single list of columns
    args          <- c(lapply( args,
        function(part) as.vector(split(part, col(part)))
    ))
    args          <- unlist(args,recursive = F)

    # statistic computations
    sds            <- unlist(lapply(args, sd))
    ns             <- unlist(lapply(args, length))
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
            fsum[i]<- ns[i] * MP.absoluteCentralMoment(unlist(args[i]))
        tsum       <- sum(fsum) 
        prWC       <- tsum * 1/(sum(ns)-length(ns)) * deltax /sdp
        assumptext <- "assumption-free"
    }
    rdWC           <- round(sdp, -log10(prWC * 1.0001 )+0.5)
    # c. Best-case instrinsinc precision
    prBC           <- deltax / sqrt(sum(ns) - length(ns))
    rdBC           <- round(sdp, -log10(prBC * 1.0001 )+0.5)    
    # d. Middle-ground intrinsinc precision
    prMG           <- (prWC/2 + prBC/2)/2
    rdMG           <- round(sdp, -log10(prMG)+0.5)

    # output results
    if (verbose) MP.showVerbose("sdpool", sdp, deltax, prEP, rdEP, prWC, rdWC, prBC, rdBC, prMG, rdMG, assumptext)
    return(setNames( c(sdp, rdEP, rdWC, rdBC, rdMG),
        c("sdpool","EXrounded", "WCrounded", "BCrounded","MGrounded") ) 
    )
}
