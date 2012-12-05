#' @title Funnel plot
#'
#' @description
#' Creates a funnel plot of an \code{rrma} object. The y-axis shows
#' the study-level standard error (square root of the study-level
#' variances) and the x-axis shows the residuals.
#' 
#' @details
#' Residuals are coloured according to the study ID. For now the colours are
#' hardcoded to Dark2 from RColorBrewer and will only appear if there
#' are 9 or fewer studies.
#' 
#' @param x An \code{rrma} meta-analysis object.
#'
#' @return A funnel plot.
#'
#' @examples
#' data(broad)
#' m <- rrma(formula = lnorReg ~ d18OresidualMean.cent, data =
#' broad, study_id = study.ID, var_eff = vlnorReg, rho = 0.5) 
#' funnel(m)
#'
#' @export
funnel <- function(x) {
  require(RColorBrewer)
  se <- sqrt(x$vi)
  res <- x$res
  ylim <- rev(range(se))
  n.study <- length(unique(x$mf$study_id)) # will need to check that n.study isn't > 9
  if(n.study <= 9)
    col <- paste(brewer.pal(n.study, "Dark2"), "90", sep = "")
  else  # too many colours for the palette, for now, give up
    col <- rep("grey35", n.study)
  col_studies <- col[as.numeric(x$mf$study_id)]
  plot(res, se, ylim = ylim, xlab = "Residual", ylab = "Standard Error", col = col_studies, pch = 19, axes = FALSE)
  axis(1, col = "grey50")
  axis(2, col = "grey50")
  box(col = "grey50")
  abline(v = 0, lty = 2, col = "grey50")
  legend("topleft", legend = unique(x$mf$study_id), col = col, pch = 19, bty = "n", cex = 0.8, text.col = "grey20")
}

#load("/Users/seananderson/Dropbox/NESCent-extinction/meta-analysis/extinction-meta-analysis/hedges_etal_2010/broadData_noNas.rda")
#m <- rrma(formula = lnorReg ~ d18OresidualMean.cent, data = broadData_noNAs, study_id = study.ID, var_eff = vlnorReg, rho = 0.5)

