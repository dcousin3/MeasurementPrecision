#' @title Measurement Precision Toolkit
#'
#' @description
#' The measurement precision toolkit returns the value of descritive statistics
#' rounded according to the measurement precision. If measurements are performed 
#' with a certain precision, called delta_x, then the statistics derived from  
#' those measurements cannot have more than a certain precisions, computed 
#' according to the formulas underlying those statistics. 
#' The descriptive statistics for which an expression of the precision is known
#' are: 
#' For univariate statistics: 
#'   mean, sd (standard deviation),
#'   semean (standard error of the mean), ci (confidence interval),
#'   cohen.d (one-sample Cohen's d, d_1), var (variance), 
#'   t.test (one-sample t-test);
#' For bivariate statistics: 
#'   cohen.d (two-sample Cohen's d, d_p), meandiff (mean difference),
#'   t.test (two-sample t-test);
#' For multivariables: 
#'   sdpool (pooled standard deviation), F.ratio (only worst-case scenario).
#'
#' Three scenarios are considered: 
#' - Extrinsinc precision: precision is estimated
#'   according to a population point of view (uses standard error of 
#'   the statistic);
#' - Intrinsinc precision (worst-case): precision is estimated assuming
#'   systematic measurement errors and the maximal impact it can have on 
#'   the statistic;
#' - Intrinsinc precision (best-case): precision is estimated assuming non-
#'   systematic measurement errors and the root-mean-squared impact in can have;
#'
#' @param fromStatistics  a list of already computed statistics; use if you do not provide fromData;
#' @param fromData        a vector, a matrix or a dataframe containing raw data; use if you do not provide fromStatistics;
#' @param deltax          the precision of the instrument;
#' @param assumptions     boolean (TRUE to assume relevant symplifying assumptions);
#' @param verbose         boolean (TRUE to display a human-readable output);
#' @param gamma           for confidence intervals, the coverage level  (default if omitted 95\%);
#' @param mu0             for the one-sample cohen.d and one-sample t.test, the mean of reference;
#'
#' @usage
#' roundMP.mean    (fromStatistics||fromData, deltax, ...)
#' roundMP.sd      (fromStatistics||fromData, deltax, ...)
#' roundMP.var     (fromStatistics||fromData, deltax, ...)
#' roundMP.semean  (fromStatistics||fromData, deltax, ...)
#' roundMP.ci      (fromStatistics||fromData, deltax, gamma, ...)
#' roundMP.cohen.d (fromStatistics||fromData, deltax, mu0, ...)
#' roundMP.t.test  (fromStatistics||fromData, deltax, mu0, ...)
#' roundMP.meandiff(fromStatistics||fromData, deltax, ...)
#' roundMP.cohen.d (fromStatistics||fromData, deltax, ...)
#' roundMP.t.test  (fromStatistics||fromData, deltax, ...)
#' roundMP.sdpool  (fromStatistics||fromData, deltax, ...)
#' roundMP.F.ratio (fromData, deltax, ...)
#'
#' @return            the summary statistic and its value rounded based 
#'                    on the measurement precision
#'
#' @details
#' These functions returns a summary statistic which is rounded
#' according to the measurement's precision.
#'
#' @author Denis Cousineau, \email{denis.cousineau@@uottawa.ca}
#' @references \url{https://.../...}
#' @keywords Measurement precision; rounding
#'
#' @examples
#' # define a vector (it could be a 1-colum matrix or a one-column data.frame)
#' x1 <- c(3,4,5)
#'
#' # get the rounded mean assuming that the instrument is precise to +or- 1
#' roundMP.mean(fromData = x1, deltax = 1)
#' roundMP.mean(fromStatistics = list(mean = 4, sd = 1, n = 3), deltax = 1, verbose = TRUE)
#'
#' # get the rounded standard error, the rounded confidence intervals
#' roundMP.semean(fromData = x1, deltax = 1)
#' roundMP.ci(fromData = x1, deltax = 1)
#'
#' # get the rounded mean difference between two vectors;
#' x2 <- c(5,7,9)
#' roundMP.meandiff(fromData = cbind(x1,x2), deltax = 1)
#'
#' # get the rounded F ratio, the rounded Cohen's d, and the rounded t-test
#' # for the last, do not assume symplifiying assumptions
#' roundMP.F.ratio(fromData = cbind(x1,x2), deltax = 1)
#' roundMP.cohen.d(fromData = cbind(x1,x2), deltax = 1)
#' roundMP.t.test( fromData = cbind(x1,x2), deltax = 1, assumptions = FALSE)
#'
#' # The F ratio and the pooled standard deviation take any number of columns
#' x3 <- c(2,5,9,11,25)
#' roundMP.sdpool( fromData = list(x1,x2,x3), deltax = 1)
#' roundMP.F.ratio(fromData = list(x1,x2,x3), deltax = 1)
#' # the F ratio only works with fromData.


# the code is long but it can be run with vectors, 
# column matrices or data.frames


#########################################################
# some subsidiary functions 
#########################################################

# extract column(s) from the "fromData" parameter
# reporting everything into a matrix
MP.getData <- function(dta, number) {
    switch(EXPR = number,
        "1" = {
            # looking for a single vector
            if (MP.rowLengths(dta) != 1) stop("there is no single vector or one-column matrix/dataframe provided... Exiting.")
            if (is.vector(dta))          res <- as.matrix(dta)
            else if (is.matrix(dta))     res <- as.matrix(dta[,1])
            else if (is.data.frame(dta)) res <- as.matrix(dta[,1])
            else if (is.list(dta))       res <- as.matrix(dta[[1]])
            else stop("Weird: unknown data structure... Exiting.") 
            },
         "2" = {
            # looking for two vectors
            if (MP.rowLengths(dta) != 2) stop("there is no two-column matrix/dataframe provided... Exiting.")
            if (is.matrix(dta))          res <- list(dta[,1],dta[,2])
            else if (is.data.frame(dta)) res <- list(dta[,1],dta[,2])
            else if (is.list(dta))       res <- dta
            else stop("Weird: unknown data structure... Exiting.") 
        },
        "1or2" = {
            # looking for 1 or 2 vectors
            if ((MP.rowLengths(dta) > 2)||(MP.rowLengths(dta) < 1 )) stop("there is not one or two column(s) in the data provided... Exiting.")
            if (MP.rowLengths(dta) == 1) {
                if (is.vector(dta))          res <- list(dta)
                else if (is.matrix(dta))     res <- list(dta[,1])
                else if (is.data.frame(dta)) res <- list(dta[,1])
                else if (is.list(dta))       res <- dta
                else stop("Weird: unknown data structure... Exiting.") 
            } else {
                if (is.matrix(dta))          res <- list(dta[,1],dta[,2])
                else if (is.data.frame(dta)) res <- list(dta[,1],dta[,2])
                else if (is.list(dta))       res <- dta
                else stop("Weird: unknown data structure... Exiting.") 
            }
        },
        "any" = {
            # looking for unspecified number of vectors
            if (MP.rowLengths(dta) < 1)  stop("there is no data provided... Exiting.")
            if (is.matrix(dta)) {
                res <- list()
                for(i in 1:MP.rowLengths(dta)) res[[i]] <- dta[,i]
            } else if (is.data.frame(dta)) {
                res <- list()
                for(i in 1:MP.rowLengths(dta)) res[[i]] <- dta[,i]
            } else if (is.list(dta))         res <- dta
            else stop("Weird: unknown data structure... Exiting.") 
        },
        stop("Weird: unknown number of columns required...")
    )
    res
}

# verify that the named statistics are in the list 
# or else issue a message
MP.vfyStat <- function(statlist, statname) {
    if (!(all(statname %in% names(statlist))))
      stop(paste("The list of statistics is incomplete, we need: ", 
                 paste(statname, collapse = " "), sep = "")
           )
    statlist
}

# show some description of the results
MP.showVerbose <- function(fct, mn, deltax, prEP, rdEP, prWC, rdWC, prBC, rdBC, assumptxt) {
    s1<-function() {rep(" ",18-nchar(fct))}
    s2<-function() {rep(" ",5)}
    cat(rep("-",70),"\n",sep="")
    cat(fct," of input is:             ",   s1(),mn,"\n",sep="")
    cat("delta_x of instrument is:              ",s2(),deltax,"\n",sep="")
    cat("EXTRINSINC PRECISION:  (this result is based on the standard error of the",fct,")\n")
    cat("  - precision for ",fct," is:    ",s1(), prEP,"\n",sep="")
    cat("  - rounded ",fct," of input is: ",s1(),rdEP,"\n",sep="")
    cat("WORST-CASE INTRINSINC PRECISION: (this result is", assumptxt,")","\n")
    cat("  - precision for ",fct," is:    ",s1(), prWC,"\n",sep="")
    cat("  - rounded ",fct," of input is: ",s1(),rdWC,"\n",sep="")
    cat("BEST-CASE INTRINSINC PRECISION: (this result is", assumptxt,")","\n")
    cat("  - precision for ",fct," is:    ",s1(), prBC,"\n",sep="")
    cat("  - rounded ",fct," of input is: ",s1(),rdBC,"\n",sep="")
    cat(rep("-",70),"\n",sep="")
}

# the absolute central moment of order 1
MP.absoluteCentralMoment <- function(x) {
        sum(abs(x-mean(x)))/length(x)
}

# a function to get the number of columns whatever the data type
MP.rowLengths <- function(x) { 
    if (is.null(x)) {
        0
    } else if (is.list(x)) { #must be done early as lists are also vectors...
        sum(unlist(lapply(x, MP.rowLengths)))
    } else if (is.vector(x)) {
        1
    } else if (is.matrix(x)) {
        dim(x)[2] 
    } else if (is.data.frame(x)) {
        dim(x)[2]
    } else stop("Weird: unknown data structure (not vector and not one-column matrix/data.frame) and not list of these...") 
}


