#' @name predict.rrma
#' @title Predict method for robust meta-regression fits
#' @S3method predict rrma
#'
#' @description
#' \code{rrma} returns predicted values and error
#' from robust variance meta-regression models 
#' 
#' @details predict.rrma is used to generate new predicted variables and error from a data set.  Note
#' incorporating study and data point level variation is not yet implemented.
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
#' @param interval Type of interval calculation.
#' @param na.action function determining what should be done with 
#' missing values in \code{newdata}. The default is to predict \code{NA}.
#' @param level Tolerance/confidence level
#' 
#' @export
#' @return A vector with predicted values or a data frame with predicted and error values
#'
#' @examples
#'
#' data(broad)
#' m <- rrma(formula = lnorReg ~ d18OresidualMean.cent, data =
#' broad, study_id = study.ID, var_eff = vlnorReg, rho = 0.5) 
#' m
#'
#' pred <- predict(m, se.fit=TRUE, interval="confidence")
#' plot(lnorReg ~ d18OresidualMean.cent, data=broad)
#' matplot(broad$d18OresidualMean.cent, pred$fit, col="red", lwd=2,
#' add=TRUE, type="l")
#'
#'idx <- sort(broad$d18OresidualMean.cent, index.return=TRUE)$ix
#' matplot(broad$d18OresidualMean.cent[idx], pred[idx,3:4], 
#'         type="l",lty=2, add=TRUE, col="black")

predict.rrma <- function(object, newdata=NULL, 
				se.fit=FALSE, interval=c("none", "confidence", "prediction"),
				na.action=na.pass, level=0.95) {

  ret <- list()
  form <- as.formula(object$call$formula)
  
  if (missing(newdata) || is.null(newdata)) newdata <- object$mf
  
  newdata[[ as.character(form[[2]]) ]] <- NA
  mm <- model.frame(form, newdata, na.action=na.action)
  X <- model.matrix(form, mm)
  	
  	
  beta <- object$est[,2]
  	
  pred <- X %*% beta
  if(!se.fit & "none" %in% interval) return(pred)
  	
  se.fit <- diag(as.matrix(X) %*% tcrossprod(vcov(object), X))
  	
  if("none" %in% interval) return(data.frame(fit = pred, se.fit = se.fit))
  	
  tfrac <- qt(1-(1-level)/2, nrow(object$input_data))
  lwr <- pred - tfrac* sqrt(se.fit)
  upr <- pred + tfrac* sqrt(se.fit)
    
  if("prediction" %in% interval) {
    lwr <- lwr - tfrac*sqrt(object$tau_sq_est)[1,1] #added matrix notation to avoid odd error
    upr <- upr + tfrac*sqrt(object$tau_sq_est)[1,1]
    	
   }
    
   return(data.frame(fit = pred, se.fit = se.fit, lwr = lwr, upr = upr))
}

# this should be re-written to also accept new data

