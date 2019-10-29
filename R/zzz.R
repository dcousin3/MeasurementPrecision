.onLoad <- function(libname, pkgname) {
    
    # set the default scenario displayed to all:
    options(roundMP.selectedScenario = c("machine.precision", "extrinsic", "systematic", "non.systematic"))
    
    # Because the object htest returned by t.test does not contain
    # all the relevant information, we hack this function:

    # 1-make a copy of the original t test function
    t.test.old = stats::t.test

    # 2-redefine the t test to add additional information
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
                temp$ns  <- length(x)
                temp$sds <- sd(x)
            } else {
                temp$ns  <- c(length(x), length(y))
                temp$sds <- c(sd(x), sd(y))        
            }
        }
        temp
    }

    # 3-replace the function with the redefined function
    unlockBinding("t.test", as.environment("package:stats"))
    assign("t.test", t.test.new, as.environment("package:stats"))
    lockBinding("t.test", as.environment("package:stats"))

}