\name{predict.rrma}
\alias{predict.rrma}
\title{Predict method for robust meta-regression fits}
\usage{
  predict.rrma(object, newdata = NULL, se.fit = FALSE,
    na.action = na.pass, level = 0.95,
    interval = c("none", "confidence"))
}
\arguments{
  \item{object}{A meta-regression object of class
  \code{rrma}}

  \item{newdata}{Optional new data frame}

  \item{se.fit}{A switch indicating if standard errors are
  required.}

  \item{na.action}{function determining what should be done
  with missing values in \code{newdata}. The default is to
  predict \code{NA}.}

  \item{level}{Tolerance/confidence interval}

  \item{interval}{Return a confidence interval?}
}
\value{
  For prediction without standard errors: a vector with
  predicted values; for prediction with standard errors: a
  list with predicted values and standard error values;
  four prediction with confidence intervals: a data frame
  with the column names \code{fit, lwr, upr}.
}
\description{
  \code{rrma} returns predicted values and error from
  robust variance meta-regression models
}
\details{
  predict.rrma is used to generate new predicted variables
  and error from a data set. Note that incorporating study
  and data point level variation is not yet implemented.
}
\examples{
data(broad)
m <- rrma(formula = lnorReg ~ d18OresidualMean.cent, data =
broad, study_id = study.ID, var_eff = vlnorReg, rho = 0.5)

pred <- predict(m, interval = "confidence")

plot(lnorReg ~ d18OresidualMean.cent, data=broad)
matplot(broad$d18OresidualMean.cent, pred$fit, col="red", lwd=2,
add=TRUE, type="l")
idx <- sort(broad$d18OresidualMean.cent, index.return=TRUE)$ix
polygon(c(broad$d18OresidualMean.cent[idx],
rev(broad$d18OresidualMean.cent[idx])), c(pred$lwr[idx],
rev(pred$upr[idx])), col = "#00000020", border = NA)
}
\author{
  Jarrett Byrnes and Sean Anderson
}
\references{
  Hedges, L.V., Tipton, E. & Johnson, M.C. (2010). Robust
  variance estimation in meta-regression with dependent
  effect_size estimates. Res. Synth. Method., 1, 39-65.
}

