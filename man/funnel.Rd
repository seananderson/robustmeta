\name{funnel}
\alias{funnel}
\title{Funnel plot}
\usage{
  funnel(x)
}
\arguments{
  \item{x}{An \code{rrma} meta-analysis object.}
}
\value{
  A funnel plot.
}
\description{
  Creates a funnel plot of an \code{rrma} object. The
  y-axis shows the study-level standard error (square root
  of the study-level variances) and the x-axis shows the
  residuals.
}
\details{
  Residuals are coloured according to the study ID. For now
  the colours are hardcoded to Dark2 from RColorBrewer and
  will only appear if there are 9 or fewer studies.
}
\examples{
data(broad)
m <- rrma(formula = lnorReg ~ d18OresidualMean.cent, data =
broad, study_id = study.ID, var_eff = vlnorReg, rho = 0.5)
funnel(m)
}

