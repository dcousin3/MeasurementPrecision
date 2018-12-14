#' @title Measurement Precision Toolkit
#'
#' @description
#' The measurement precision toolkit returns the value of descritive statistics
#' rounded according to the measurement precision. If measurements are performed 
#' with a certain precision, called delta_x, then the statistics derived from those 
#' measurements cannot have more than a certian precisions, computed according to 
#' the formulas underlying those statistics. 
#' The descriptive statistics for which an expression of the precision is known are:
#' For univariate statistics: mean, sd (standard deviation), semean (standard error of the mean)
#' ci (confidence interval), cohen.d (one-sample Cohen's d, d_1), var (variance), and 
#' t.test (one-sample t-test);
#' For bivariate statistics: cohen.d (two-sample Cohen's d, d_p), meandiff (mean difference),
#' and t.test (two-sample t-test);
#' For multivariables: sdpool (pooled standard deviation) and F.ratio (one-way F test).
#'
#' @param x           a vector of numbers;
#' @param y           (optional) a second vector for bivariate statistics;
#' @param ...         (optional) any number of vectors for multivariate statistics;
#' @param fct         the summary statistic function between quotes: mean, sd, semean, cohen.d, ci, var, t.test, meandiff, F.ratio, sdpool;
#' @param deltax      the precision of the instrument;
#' @param assumptions boolean (TRUE to assume relevant symplifying assumptions);
#' @param verbose     boolean (TRUE to display a human-readable output);
#' @param gamma       the coverage level for confidence intervals (default if omitted 95\%);
#' @param mu0         for the one-sample cohen.d and one-sample t.test, provide the mean of reference;
#'
#' @usage
#' roundMP         (x, fct, deltax, ...)
#' roundMP.mean    (x, deltax, ...)
#' roundMP.sd      (x, deltax, ...)
#' roundMP.var     (x, deltax, ...)
#' roundMP.semean  (x, deltax, ...)
#' roundMP.ci      (x, deltax, gamma, ...)
#' roundMP.cohen.d (x, deltax, mu0, ...)
#' roundMP.t.test  (x, deltax, mu0, ...)
#' roundMP.meandiff(x, y, deltax, ...)
#' roundMP.cohen.d (x, y, deltax, ...)
#' roundMP.t.test  (x, y, deltax, ...)
#' roundMP.sdpool  (x, y, ..., deltax, ...)
#' roundMP.F.ratio (x, y, ..., deltax, ...)
#'
#' @return            the summary statistic and its value rounded based on the measurement precision
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
#' roundMP(x1, 'mean', deltax = 1)
#' roundMP.mean(x1, deltax = 1)
#' roundMP.mean(x1, deltax = 1, verbose = TRUE)
#'
#' # get the rounded standard error, the rounded confidence intervals
#' roundMP.semean(x1, deltax = 1)
#' roundMP.ci(x1, deltax = 1)
#'
#' # get the rounded mean difference between two vectors;
#' x2 <- c(5,7,9)
#' roundMP.meandiff(x1, y = x2, deltax = 1)
#'
#' # get the rounded F ratio, the rounded Cohen's d, and the rounded t-test
#' # for the last, do not assume symplifiying assumptions
#' roundMP.F.ratio(x1, y = x2, deltax = 1)
#' roundMP.cohen.d(x1, y = x2, deltax = 1)
#' roundMP.t.test( x1, y = x2, deltax = 1, assumptions = FALSE)
#'
#' # The F ratio and the pooled standard deviation take any number of columns
#' x3 <- c(2,5,9,11,25)
#' roundMP.sdpool( x1, x2, x3, deltax = 1)
#' roundMP.F.ratio(x1, x2, x3, deltax = 1)
#'


#########################################################
# documentation meant for roxygen2
#########################################################

# yaurai aussi: usage (des exemples) arguments (les param?es)





# the code is long because it can be run with vectors, 
# column matrices or data.frame with a single column...



#########################################################
# the wrapper function
#########################################################

#' @export
roundMP <- function(
                x,                  # a vector/matrix/data.frame
                fct,                # name of the summary statistic
                deltax,             # the precision of the instrument
                y           = NULL, # may be used for bivariate statistics
                assumptions = TRUE, # uses simplifying assumptions?
                verbose     = FALSE,# provides detailed output?
                ...                 # some requires extra parameters
) { 
    # The list of summary statistics for which the
    # impact of measurement precision has been found
    MP.listfcts = list(
        "univariable"  = c("mean","semean","ci","cohen.d","sd","var","t.test"),
        "bivariable"   = c("meandiff","cohen.d","t.test"),
        "multivariable"= c("sdpool", "F.ratio")
    )

    # check that the function is in the list of known functions
    if (!(fct %in% unlist(MP.listfcts))) 
        stop(paste("no known solution for function",fct, "or unknown function"))

    # if so, call the right function based on the shape of the input
    # and undo data.frame structure
    m = sapply(MP.listfcts, function(t) fct %in% t )
    if ((m["univariable"])&&(MP.rowLengths(x)==1)) {
        do.call(paste('roundMP',fct, sep="."),
            list(x=x, deltax=deltax, ..., assumptions=assumptions, verbose = verbose) )
    } else if ((m["bivariable"])&&((MP.rowLengths(x) == 2))) {
        do.call(paste('roundMP',fct, sep="."),
            list(x=x[,1], y=x[,2], deltax=deltax, ..., assumptions=assumptions, verbose = verbose) )
    } else if ((m["bivariable"])&&((MP.rowLengths(x)==1)&&(MP.rowLengths(y)==1))) {
        do.call(paste('roundMP',fct, sep="."),
            list(x=x, y=y, deltax=deltax, ..., assumptions=assumptions, verbose = verbose) )    
    } else if ((m["multivariable"])&&(MP.rowLengths(x)>1)) {
        do.call(paste('roundMP',fct, sep="."),
            list(x=x, y=y, ..., deltax=deltax, assumptions=assumptions, verbose = verbose) )        
    } else {
        stop("it seems that the number of columns of the data does not match the function... Exiting")
    }
}


#########################################################
# some subsidiary functions 
#########################################################

# show some description of the results
MP.showVerbose <- function(fct, mn, deltax, pr, rd, assumptxt) {
    s1<-function() {rep(" ",12-nchar(fct))}
    s2<-function() {rep(" ",5)}
    cat(rep("-",60),"\n",sep="")
    cat(fct," of input is:         ",   s1(),mn,"\n",sep="")
    cat("delta_x of instrument is:    ",s2(),deltax,"\n",sep="")
    cat("precision for ",fct," is:    ",s1(), pr,"\n",sep="")
    cat("rounded ",fct," of input is: ",s1(),rd,"\n",sep="")
    cat("this results is ", assumptxt,".\n",sep="")
    cat(rep("-",60),"\n",sep="")
}

# the absolute central moment of order 1
MP.absoluteCentralMoment <- function(x) {
        sum(abs(x-mean(x)))/length(x)
}

# a function to get the number of columns whatever the data type
MP.rowLengths <- function(x) { 
    if (is.null(x)) {
        0
    } else if (is.vector(x)) {
        1
    } else if (is.matrix(x)) {
        dim(x)[2] 
    } else if (is.data.frame(x)) {
        dim(x)[2]
    } else {stop("Weird: unknown data structure (not vector and not one-column matrix/data.frame)") }
}

# a function that revert to vector whatever the data type containing the column
MP.flatten <- function(x) {
    if (is.null(x)) {
        NULL
    } else if (is.vector(x)) {
        x
    } else if (is.matrix(x)) {
        as.vector(x) 
    } else if (is.data.frame(x)) {
        as.vector(as.matrix(x))
    } else {stop("Weird: unknown data structure (not vector and not one-column matrix/data.frame)") }
}

