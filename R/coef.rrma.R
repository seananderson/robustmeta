#' @S3method coef rrma
coef.rrma <- function(x) {
	ret <- x$est[,2]	
	names(ret) <- x$est[,1]
	ret
}

