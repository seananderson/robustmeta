\name{coefplot}
\alias{coefplot}
\title{Coefficient plot of an rrma object}
\usage{
  coefplot(x)
}
\arguments{
  \item{x}{An \code{rrma} meta-analysis object.}
}
\value{
  A coefficient plot.
}
\description{
  Creates a coefficient plot of an \code{rrma} object.
}
\examples{
data(broad)
m <- rrma(formula = lnorReg ~ d18OresidualMean.cent, data =
broad, study_id = study.ID, var_eff = vlnorReg, rho = 0.5)
coefplot(m)
}

