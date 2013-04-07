\name{plot.rrma}
\alias{plot.rrma}
\title{Forest plot of an rrma object}
\usage{
  plot.rrma(x, xlab = "Effect size", xlim = NULL,
    mar = c(4, 11, 1, 1), cex.axis = 1,
    bg_col = "#00000020", bg_line_col = "white", pch = 21,
    zero_lty = 2, diamond_col = "grey40")
}
\description{
  Creates a forest plot of an rrma object, which shows the
  effect sizes and confidence intervals grouped by study
  along with the robust meta-analytic mean and confidence
  intervals.
}
\examples{
data(broad)
m <- rrma(formula = lnorReg ~ d18OresidualMean.cent, data =
broad, study_id = study.ID, var_eff = vlnorReg, rho = 0.5)
plot(m)
}
\references{
  Hedges, L.V., Tipton, E. & Johnson, M.C. (2010). Robust
  variance estimation in meta-regression with dependent
  effect_size estimates. Res. Synth. Method., 1, 39-65.
}

