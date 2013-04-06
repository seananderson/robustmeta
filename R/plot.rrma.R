#' @title  Forest plot of an rrma object
#'
#' @description
#' Creates a forest plot of an rrma object, which shows the
#' effect sizes and confidence intervals grouped by study along with
#' the robust meta-analytic mean and confidence intervals.
#'
#' @references
#' Hedges, L.V., Tipton, E. & Johnson, M.C. (2010). Robust variance
#' estimation in meta-regression with dependent effect_size estimates.
#' Res. Synth. Method., 1, 39-65.  
#' 
#' @export plot rrma
#' @S3method plot rrma
#' 
#' @examples
#' data(broad)
#' m <- rrma(formula = lnorReg ~ d18OresidualMean.cent, data =
#' broad, study_id = study.ID, var_eff = vlnorReg, rho = 0.5)
#' plot(m)

plot.rrma <- function(x, xlab = "Effect size", xlim = NULL, mar = c(4, 11, 1, 1), cex.axis = 1, bg_col = "#00000020", bg_line_col = "white", pch = 21, zero_lty = 2, diamond_col = "grey40") {

  d <- x$input_data

  par(mar = mar)
  #par(mgp = c(2, 0.5, 0))
  #par(tck = -0.03)

  require(plyr)

  d <- ddply(d, "study", transform, mean_eff = mean(effect_size))
  d <- d[order(-d$mean_eff, d$effect_size),] 
  d <- transform(d, label_show = c(1, diff(as.numeric(d$study))), label = study_id)
  d$label <- as.character(d$label)
  for(i in 1:(nrow(d)-1)) {
    if(d$label[i] == d$label[i+1]) d$label[i] <- ""
  }
  d <- d[order(-d$mean_eff, -d$effect_size),] 
  d <- transform(d, plot_ID = 1:nrow(d))
  l <- d$effect_size - 1.96*sqrt(d$var_eff_size)
  u <- d$effect_size + 1.96*sqrt(d$var_eff_size)
  if(is.null(xlim)) xlim <- c(min(l), max(u))
  with(d, plot.default(effect_size, plot_ID, xlim = xlim, ylim = c(-2, nrow(d)+0.5), axes = FALSE, ylab = "", xlab = "Effect size", yaxs = "i", type = "n"))

# background lines and rectangles:
  abline(h = 1:nrow(d), col = bg_line_col)
  rects <- subset(d, label != "")[,c("plot_ID")]
  odd <- seq(1, nrow(d)+1, 2) # 19 is random... just needs to be big enough
  if(length(rects) %% 2 == 1) rects <- c(rects, nrow(d) + 1) # if odd, add rect on end
  for(i in odd) rect(par("usr")[1], rects[i]-0.5, par("usr")[2], rects[i+1]-0.5, col = bg_col, border = FALSE)
  with(d, segments(l, plot_ID, u , plot_ID))
  with(d, points(effect_size, plot_ID, pch = pch))

  box()
  axis(1, cex = cex.axis)
  abline(v = 0, lty = zero_lty)

# center citation labels:
  study_labs <- subset(d, label != "")[,c("label", "plot_ID")]
  study_labs$lab_loc <- NA
  for(i in 1:(nrow(study_labs) - 1)) study_labs$lab_loc[i] <- mean(c(study_labs$plot_ID[i], study_labs$plot_ID[i+1]))
  study_labs$lab_loc[nrow(study_labs)] <- study_labs$plot_ID[nrow(study_labs)] + 0.5
  study_labs$lab_loc <- study_labs$lab_loc - 0.5

  with(study_labs, axis(2, at = lab_loc, labels = label, las = 1, tick = FALSE, cex.axis = cex.axis))


  if("intercept" %in% x$est[,1]) {
   ci.lb.n <- which(names(x$est) == "ci.lb")
   ci.ub.n <- which(names(x$est) == "ci.ub")

  est <- x$est[1, 2]
  ci.lb <- x$est[1, ci.lb.n]
  ci.ub <- x$est[1, ci.ub.n]

  polygon(c(est, ci.lb, est, ci.ub), -c(0.7, 1, 1.3, 1), col = diamond_col)
  axis(2, at = -0.75, labels = "Meta-analytic mean", las = 1, tick = FALSE, cex.axis = cex.axis)
  }else{
    warning("No meta-analytic mean (intercept) to plot.")
  }
}

