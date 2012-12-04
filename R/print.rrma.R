#' @S3method print rrma
print.rrma <- function(x, digits = max(3, getOption("digits") - 3)){
	cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
    cat(paste("QE =", round(x$Qe, digits)), "\n\n")
    print(x$est)
    cat("\n", paste("tau^2 =", round(x$tau_sq_est, digits)), "\n", sep="")
    cat(paste("sample size =", x$n), "\n")
    cat(paste("Number of studies =", x$k), "\n\n")
}

