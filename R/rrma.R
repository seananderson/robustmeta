#' @title Robust random-effect meta-regression with dependent effect sizes
#'
#' @description
#' \code{rrma} returns robust variance estimates from random effect
#' meta-regression with dependent effect-size estimates. 
#' 
#' @details
#' This is an implementation based on the \code{robust.se} function from
#' Hedges et al. (2010).
#'
#' @references
#' Hedges, L.V., Tipton, E. & Johnson, M.C. (2010). Robust variance
#' estimation in meta-regression with dependent effect_size estimates.
#' Res. Synth. Method., 1, 39-65.  
#'
#' @author Main body of code is from the appendix of Hedges et al.
#' (2010) and written by the paper authors.
#' Code adapted to an R package with a formula interface and various
#' convenience functions by Sean Anderson and Jarrett Byrnes.
#'
#' @param formula The meta-regression formula of the form \code{y ~ x1 +
#' x2...}, where \code{y} is the effect size. An intercept only meta-analysis
#' can be performed with \code{y ~ 1}. The intercept can be excluded with
#' \code{y ~ -1 + x1 + x2 ...}
#' @param data The input data frame. 
#' @param study_id The study IDs. Can be in any form (character,
#' numeric, factor). Will be converted to factor and then numeric form
#' internally.
#' @param var_eff_size The variance on each effect size.
#' @param rho The assumed correlation between effect sizes of the same
#' study including correlation induced by the random effects. Hedges
#' et al. (2010) suggest conducting a sensitivity analysis on
#' \code{rho}.
#' 
#' @export
#' @return A list with the robust covariance matrix, the Qe statistic,
#' the tau_sq estimate (variance between studies), a data frame of
#' the robust coefficient estimates and standard errors, and a number
#' of diagnostics and input values.
#'
#' @examples
#' data(broad)
#' (m <- rrma(formula = lnorReg ~ d18OresidualMean.cent, data =
#' broad, study_id = study.ID, var_eff = vlnorReg, rho = 0.5))
#'
rrma <- function(formula, data, study_id, var_eff_size, rho) {
# [r]obust [r]andom effect [m]eta-[a]nalysis regression Hedges et al. 2010
# code to parse the formula and create the required data structure:
  require(plyr) # using plyr to save time because I already coded it this way
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "study_id", "var_eff_size"), 
    names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mf <- data.frame(effect_size = mf[,1], model.matrix(formula, mf), 
  					 study_id=mf[["(study_id)"]],
  					 var_eff_size=mf[["(var_eff_size)"]])
  names(mf) <- gsub("\\(", "", names(mf))
  names(mf) <- gsub("\\)", "", names(mf))
  names(mf) <- gsub("X.Intercept.", "intercept", names(mf)) #for later use
  ind_var_names <- names(mf)[-c(1, match(c("study_id", "var_eff_size"), names(mf)))]
  mf <- ddply(mf, "study_id", transform, k = length(study_id), 
    mean_v = mean(var_eff_size), s = sqrt(var_eff_size))
  mf <- ddply(mf, "study_id", transform, weights = 1/(k * mean_v))
  #mf$intercept <- 1
  names(mf)[1] <- "effect_size"
  mf$study <- as.numeric(factor(mf$study_id))
  names_design_matrix <- c("study",  ind_var_names)
  names_data_matrix <- c("study", "k", "mean_v", "weights", 
    "effect_size", "var_eff_size", "s")

# now create 2 data frames: input_data and X.full
  input_data <- mf[,names_data_matrix]
  X_full <- mf[,names_design_matrix]

# and proceed with the function from Hedges et al. 2010:
  p <- ncol(X_full) - 2 # number of non-intercept betas
  N <- max(input_data$study) # number of studies
  sumXWX <- 0
  sumXWy <- 0
  sumXWJWX <- 0
  sumXWVWX <- 0
  sumXW.sig.m.v.WX <- 0
  for (i in (1:N)) { # loop through each study
    tab <- input_data[input_data$study == i, ] # get the meta-analysis data for this study
    W <- diag(tab$weights, tab$k[1]) # weights are 1 / (k * mean_v)
    tab2 <- X_full[X_full$study == i, ] # design matrix for this study
    tab3 <- cbind(tab2[-c(1)]) # remove study column
    X <- data.matrix(tab3) # convert to matrix
    dimnames(X) <- NULL # clean column and row names
    y <- cbind(tab$effect_size) # create effect_size vector
    one <- cbind(rep(1, tab$k[1])) # a matrix of ones k times
    J <- one %*% t(one) # a k x k matrix of ones
    sigma <- (tab$s %*% t(tab$s)) # variance-covariance
    vee <- diag(tab$s^2, tab$k[1]) # variance diagonals
    SigmV <- sigma - vee # the covariance elements
    sumXWX <- sumXWX + t(X) %*% W %*% X # first half of equation 3
    sumXWy <- sumXWy + t(X) %*% W %*% y # second half of equation 3
    sumXWJWX <- sumXWJWX + t(X) %*% W %*% J %*% W %*% X # equation 4; needed for denominator of equation 15
    sumXWVWX <- sumXWVWX + t(X) %*% W %*% vee %*% W %*% X  # needed for equation 15
    sumXW.sig.m.v.WX <- sumXW.sig.m.v.WX + t(X) %*% W %*% SigmV %*% W %*% X # needed for equation 15
  }
  b <- solve(sumXWX) %*% sumXWy # solve for the (not-robust) beta means
  X <- data.matrix(X_full[-c(1)])
  dimnames(X) <- NULL
  input_data$pred <- X %*% b # predictions
  input_data$e <- input_data$effect_size - input_data$pred # e = residuals
  W <- diag(input_data$weights) # matrix with weights on diagonal
  sumW <- sum(input_data$weights) # sum of all weights across all effects
  Qe <- t(input_data$e) %*% W %*% input_data$e # equation 14, the Qe statistic
  denom <- sumW - sum(diag(solve(sumXWX) %*% sumXWJWX)) # denominator of equation 15
  term1 <- (Qe - N + sum(diag(solve(sumXWX) %*% sumXWVWX)))/denom # part of equation 15
  term2 <- (sum(diag(solve(sumXWX) %*% sumXW.sig.m.v.WX)))/denom # part of equation 15
  tau.sq1 <- term1 + rho * term2 # the full equation 15
  tau.sq <- ifelse(tau.sq1 < 0, 0, tau.sq1) # don't allow negative tau
  input_data$r.weights <- 1/(input_data$k * (input_data$mean_v + tau.sq)) # weights adjusted for tau.sq; [r]obust weights; equation 17
  sumXWX.r <- 0
  sumXWy.r <- 0
  for (i in (1:N)) { # now repeat beta estimation with robust weights
    tab <- input_data[input_data$study == i, ]
    W <- diag(tab$r.weights, tab$k[1])
    tab2 <- X_full[X_full$study == i, ]
    tab3 <- cbind(tab2[-c(1)])
    X <- data.matrix(tab3)
    dimnames(X) <- NULL
    y <- cbind(tab$effect_size)
    sumXWX.r <- sumXWX.r + t(X) %*% W %*% X
    sumXWy.r <- sumXWy.r + t(X) %*% W %*% y
  }
  b.r <- solve(sumXWX.r) %*% sumXWy.r # solve for robust beta means
  X <- data.matrix(X_full[-c(1)])
  dimnames(X) <- NULL
  input_data$pred.r <- X %*% b.r # predictions from robust betas
  input_data$e.r <- cbind(input_data$effect_size) - input_data$pred.r # residuals from robust betas
  sumXWeeWX.r <- 0
  for (i in (1:N)) { # now get estimated variance-covariance matrix using the robust residuals
    tab <- input_data[input_data$study == i, ]
    sigma.hat.r  <- tab$e.r %*% t(tab$e.r)
    W <- diag(tab$r.weights, tab$k[1])
    tab2 <- X_full[X_full$study == i, ]
    tab3 <- cbind(tab2[-c(1)])
    X <- data.matrix(tab3)
    dimnames(X) <- NULL
    sumXWeeWX.r <- sumXWeeWX.r + t(X) %*% W %*% sigma.hat.r %*% 
      W %*% X
  }
  VR.r <- solve(sumXWX.r) %*% sumXWeeWX.r %*% solve(sumXWX.r) # the robust (approximated) variance-covariance matrix
  SE <- c(rep(0, p + 1))
  for (i in (1:(p + 1))) { # now convert the variances (the diagonal) to standard errors
    SE[i] <- sqrt(VR.r[i, i]) * sqrt(N/(N - (p + 1)))
  }
  labels <- names(X_full)[-1]
  output <- data.frame(beta = labels, estimate = b.r, SE = SE) 
  
  #fill out the table with z-scores and CI information to match rma
    crit <- qt(0.05/2, df=(N-(p+1)), lower.tail=FALSE )

  output <- cbind(output, with(output, data.frame(
  	tval = estimate/SE,
  	pval = pt(abs(estimate/SE), df = (N - (p + 1)), lower.tail=F)*2,
  	ci.lb = estimate - crit*SE,
  	ci.ub = estimate + crit*SE
  )))
  
  rownames(VR.r) <- colnames(VR.r) <- output[,1]

  ret <- list(VR_r = VR.r, Qe = Qe, `tau_sq_est` = tau.sq, `est` =
              output, W=W, X_full = X_full, input_data=input_data,
              rho=rho, mf=mf, call=cl, n=nrow(input_data), pred =
              input_data$pred.r, res = input_data$e.r, k = N, yi =
              input_data$effect_size, vi = input_data$var_eff_size,
              df = N-(p+1), crit = crit)
  
  class(ret) <- "rrma"
  ret
}

