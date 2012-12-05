#' @title Coefficient plot of an rrma object
#'
#' @description
#' Creates a coefficient plot of an \code{rrma} object.
#' 
#' @param x An \code{rrma} meta-analysis object.
#' @return A coefficient plot.
#' @export
#' @examples
#' data(broad)
#' m <- rrma(formula = lnorReg ~ d18OresidualMean.cent, data =
#' broad, study_id = study.ID, var_eff = vlnorReg, rho = 0.5) 
#' coefplot(m)
coefplot <- function(x) {
# coefficient plot from rrma
  est <- x$est$estimate
  ci.lb <- x$est$ci.lb
  ci.ub <- x$est$ci.ub
  coef.names <- x$est$beta
  xlim <- range(ci.lb, ci.ub)
  plot(est, as.numeric(coef.names), xlim = xlim, xlab = "Estimate", ylab = "", yaxt = "n")
  segments(ci.lb, as.numeric(coef.names), ci.ub, as.numeric(coef.names))
  abline(v = 0, lty = 2, col = "grey50")
  axis(2, at = as.numeric(coef.names), labels = abbreviate(coef.names, minlength = 5), las = 1)
}

