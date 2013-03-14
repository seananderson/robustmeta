\name{rrma}
\alias{rrma}
\title{Robust random-effect meta-regression with dependent effect sizes}
\usage{
  rrma(formula, data, study_id, var_eff_size, rho)
}
\arguments{
  \item{formula}{The meta-regression formula of the form y
  ~ x1 + x2..., where y is the effect size.}

  \item{data}{The input data frame.}

  \item{study_id}{The study IDs. Can be in any form
  (character, numeric, factor). Will be converted to factor
  and then numeric form internally.}

  \item{var_eff_size}{The variance on each effect size.}

  \item{rho}{The assumed correlation between effect sizes
  of the same study including correlation induced by the
  random effects. Hedges et al. (2010) suggest conducting a
  sensitivity analysis on \code{rho}.}
}
\value{
  A list with the robust covariance matrix, the Qe
  statistic, the tau_sq estimate (variance between
  studies), and a data frame of the robust coefficient
  estimates and standard errors.
}
\description{
  \code{rrma} returns robust variance estimates from random
  effect meta-regression with dependent effect-size
  estimates.
}
\details{
  This is an implementation based on the \code{robust.se}
  function from Hedges et al. (2010).
}
\examples{
data(broad)
m <- rrma(formula = lnorReg ~ d18OresidualMean.cent, data =
broad, study_id = study.ID, var_eff = vlnorReg, rho = 0.5)
m
}
\author{
  Main body of code is from the appendix of Hedges et al.
  (2010) and written by the paper authors. Code adapted to
  an R package with a formula interface and various
  convenience functions by Sean Anderson and Jarrett
  Byrnes.
}
\references{
  Hedges, L.V., Tipton, E. & Johnson, M.C. (2010). Robust
  variance estimation in meta-regression with dependent
  effect_size estimates. Res. Synth. Method., 1, 39-65.
}

