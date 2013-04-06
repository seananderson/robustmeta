#' @name predict.rrma
#' @title Predict method for robust meta-regression fits
#' @export
#' @S3method predict rrma
#'
#' @description
#' \code{predict.rrma} returns predicted values and error
#' from robust variance meta-regression models 
#' 
#' @details predict.rrma is used to generate new predicted variables
#' and error from a data set. Note that incorporating study and data
#' point level variation is not yet implemented.
#'
#' @references
#' Hedges, L.V., Tipton, E. & Johnson, M.C. (2010). Robust variance
#' estimation in meta-regression with dependent effect_size estimates.
#' Res. Synth. Method., 1, 39-65.  
#'
#' @author Jarrett Byrnes and Sean Anderson 
#'
#' @param object A meta-regression object of class \code{rrma}
#' @param newdata Optional new data frame 
#' @param se.fit A switch indicating if standard errors are required.
#' @param na.action function determining what should be done with 
#' missing values in \code{newdata}. The default is to predict \code{NA}.
#' @param level Tolerance/confidence interval
#' @param interval Return a confidence interval?
#' 
#' @return For prediction without standard errors: a vector with
#' predicted values;  for prediction with standard errors: a list with
#' predicted values and standard error values; four prediction with
#' confidence intervals: a data frame with the column names \code{fit,
#' lwr, upr}.
#'
#' @examples
#' data(broad)
#' m <- rrma(formula = lnorReg ~ d18OresidualMean.cent, data =
#' broad, study_id = study.ID, var_eff = vlnorReg, rho = 0.5) 
#' 
#' pred <- predict(m, interval = "confidence")
#' 
#' plot(lnorReg ~ d18OresidualMean.cent, data=broad)
#' matplot(broad$d18OresidualMean.cent, pred$fit, col="red", lwd=2,
#' add=TRUE, type="l")
#' idx <- sort(broad$d18OresidualMean.cent, index.return=TRUE)$ix
#' polygon(c(broad$d18OresidualMean.cent[idx],
#' rev(broad$d18OresidualMean.cent[idx])), c(pred$lwr[idx],
#' rev(pred$upr[idx])), col = "#00000020", border = NA)

predict.rrma <- function(object, newdata=NULL, se.fit=FALSE, na.action=na.pass, level = 0.95, interval = c("none", "confidence")) {

  ret <- list()
  form <- as.formula(object$call$formula)
  
  if (missing(newdata) || is.null(newdata)) newdata <- object$mf
  
  newdata[[ as.character(form[[2]]) ]] <- NA
  mm <- model.frame(form, newdata, na.action=na.action)
  X <- model.matrix(form, mm)
  	
  beta <- object$est[,2]
  	
  pred <- as.numeric(X %*% beta)
  names(pred) <- row.names(newdata)
  #if(!se.fit) return(pred)
  if(!se.fit & "none" %in% interval) return(pred)
  	
  se.fit <- as.numeric(diag(as.matrix(X) %*% tcrossprod(vcov(object), X)))
  	
  if("none" %in% interval) 
    return(list(fit = pred, se.fit = se.fit, df = object$df))
  	
  tfrac <- qt(1-(1-level)/2, nrow(object$input_data))
  lwr <- pred - tfrac* sqrt(se.fit)
  upr <- pred + tfrac* sqrt(se.fit)
    
  #if("prediction" %in% interval) {
    #lwr <- lwr - tfrac*sqrt(object$tau_sq_est)[1,1] #added matrix notation to avoid odd error
    #upr <- upr + tfrac*sqrt(object$tau_sq_est)[1,1]
   #}
   return(data.frame(fit = pred, se.fit = se.fit, lwr = lwr, upr = upr))
}

