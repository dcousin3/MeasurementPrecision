## FRATIO = MULTIVARIATE ################################

#' @export
roundMP.F.ratio <- function(..., deltax = NULL, assumptions = TRUE, verbose = FALSE) {
    # validation and conversion
    if (is.null(deltax)||(length(deltax)>1)) stop("deltax must receive an integer value")
    args      <- list( ... )
    args[sapply(args, is.null)] <- NULL  # removing nulls
    if (!length(args)) stop("F.ratio needs at least one vector argument")
    args       = lapply(args, as.matrix)
    # unpack all the columns into a single list of columns
    args       = c(lapply( args,
        function(part) as.vector(split(part, col(part)))
    ))
    args       = unlist(args,recursive = F)

    # statistic computations
    sds        = unlist(lapply(args, sd))
    ns         = unlist(lapply(args, length))
    sdp        = sqrt(sum((ns - 1) * sds^2)/(sum(ns) - length(ns))) 
    mns        = unlist(lapply(args, mean))    
    GM         = mean(unlist(args))
    sdm        = sd(mns)
    ssa        = sum(ns * (mns-GM)^2 )
    ss         = 1:length(ns)
    for (i in 1:length(ns)) {
      ss[i]    = sum((args[[i]]-mns[i])^2)      
    }
    sse        = sum(ss)
    fratio     = (ssa / (length(ns)-1)) / (sse/(sum(ns)-length(ns)))

    # precision computations
    fsum       = 1:length(ns)
    for (i in 1:length(ns)) {
        fsum[i]= sum(abs(ssa/sse*(args[[i]]-mns[i]) -(mns[i]-GM) ))
    }
    tsum       = sum(fsum) 
    pr         = tsum * 2/(length(ns)-1) * deltax /sdp^2
    assumptext = "assumption-free"
    pr         = pr + pr/10000 # avoid rounding errors
    rd         = round(fratio, -log10(pr)+0.5)
    # output results
    if (verbose) MP.showVerbose("fratio", fratio, deltax, pr, rd, assumptext)
    return(setNames( c(fratio,deltax,pr,rd),
        c("fratio","delta_x","precision","rounded fratio")
    ))
}

