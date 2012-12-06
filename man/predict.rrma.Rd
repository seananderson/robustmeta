\name{predict.rrma}
\alias{predict.rrma}
\title{Predict method for robust meta-regression fits}
\usage{
  predict.rrma(object, newdata = NULL, se.fit = FALSE,
    interval = c("none", "confidence", "prediction"),
    na.action = na.pass, level = 0.95)
}
\arguments{
  \item{object}{A meta-regression object of class
  \code{rrma}}

  \item{newdata}{Optional new data frame}

  \item{se.fit}{A switch indicating if standard errors are
  required.}

  \item{interval}{Type of interval calculation.}

  \item{na.action}{function determining what should be done
  with missing values in \code{newdata}. The default is to
  predict \code{NA}.}

  \item{level}{Tolerance/confidence level}
}
\value{
  A vector with predicted values or a data frame with
  predicted and error values
}
\description{
  \code{rrma} returns predicted values and error from
  robust variance meta-regression models
}
\details{
  predict.rrma is used to generate new predicted variables
  and error from a data set.  Note incorporating study and
  data point level variation is not yet implemented.
}
\examples{
data(broad)
m <- rrma(formula = lnorReg ~ d18OresidualMean.cent, data =
broad, study_id = study.ID, var_eff = vlnorReg, rho = 0.5)
m

pred <- predict(m, se.fit=TRUE, interval="confidence")
plot(lnorReg ~ d18OresidualMean.cent, data=broad)
matplot(broad$d18OresidualMean.cent, pred$fit, col="red", lwd=2,
add=TRUE, type="l")

idx <- sort(broad$d18OresidualMean.cent, index.return=TRUE)$ix
matplot(broad$d18OresidualMean.cent[idx], pred[idx,3:4],
        type="l",lty=2, add=TRUE, col="black")
}
\author{
  Jarrett Byrnes and Sean Anderson
}
\references{
  Hedges, L.V., Tipton, E. & Johnson, M.C. (2010). Robust
  variance estimation in meta-regression with dependent
  effect_size estimates. Res. Synth. Method., 1, 39-65.
}

