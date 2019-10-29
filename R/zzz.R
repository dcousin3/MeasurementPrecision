.onLoad <- function(libname, pkgname) {
    # Because the object htest returned by t.test does not contain
    # all the relevant information, we hack this function.

    # make a copy of the original t test function
    t.test.old = stats::t.test

    # redefine the t test to add additional information
    t.test.new <- function(x, y = NULL,
           alternative = c("two.sided", "less", "greater"),
           mu = 0, paired = FALSE, var.equal = FALSE,
           conf.level = 0.95, ...
    ) {
        # run the t.test as requested
        temp <- t.test.old(x, y, alternative, mu, paired, var.equal, conf.level,...)
        
        # add  attributes ss and sds to the object for true t tests
        if ((var.equal == TRUE)&(paired == FALSE)) {
            if (is.null(y)) {
                temp$ss  <- length(x)
                temp$sds <- sd(x)
            } else {
                temp$ss  <- c(length(x), length(y))
                temp$sds <- c(sd(x), sd(y))        
            }
        }
        temp
    }

    unlockBinding("t.test", as.environment("package:stats"))
    assign("t.test", t.test.new, as.environment("package:stats"))
    lockBinding("t.test", as.environment("package:stats"))

}