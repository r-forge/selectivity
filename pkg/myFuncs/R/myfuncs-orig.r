#-------------------------------------------------------------------------------
#
# function  : plot.effects
#
# purpose   : plots the levels of a factor in a nice plot.
#
# arguments : NOT ALL COMPLETE YET!!
#             link - the link function of the fits and se
#             names - the names of the levels
#             ylim - optional upper/lower limits on y axis
#             mult - 2* or 1.96* se ?
#             yform - the format of the y axis labels - defaults to 0 dp
#
#-------------------------------------------------------------------------------
plot.effects <-
function(dat, link="logit.percent", names=NULL, ylim=c(NA,NA), mult=2, yform="%4.0f",...)
{
  H <- if (link=="logit.percent") function(x) ((1+exp(x))^-1)*100
  n <- nrow(dat)
  names <- if(is.null(names)) dat$names else names
  dat$ciu <- dat$fit - mult*dat$se
  dat$cil <- dat$fit + mult*dat$se
  dat <- dat[,c("cil","fit","ciu")]
  dat <- H(dat)
  ylim <- sapply(1:2, function(i) if (is.na(ylim[i])) range(dat)[i] else ylim[i])
  plot(0,0,type="n",axes=FALSE,ylim=ylim, xlim=c(0,1),xlab="",ylab="")
  x <- 1:n/n - 1/(2*n)
  segments(x, dat[,1], x, dat[,3])
  points(x, dat[,2], pch=16)
  at <- pretty(c(as.vector(dat),ylim))
  axis(2, at=at, labels=sprintf(yform, at), las=1, tck=-0.015)
  mtext(names, 1, outer=FALSE, at=x, font=2)
  hasylab <- function(...) !all(is.na(pmatch(names(list(...)), c("ylab","main"))))
  if (hasylab(...)) title(...)
}

#-------------------------------------------------------------------------------
#
# function  : polygon.fit (dat, mult)
#
# purpose   : plots the levels of a factor in a nice plot.
#
# arguments : NOT ALL COMPLETE YET!!
#             link - the link function of the fits and se
#             names - the names of the levels
#             ylim - optional upper/lower limits on y axis
#             mult - 2* or 1.96* se ?
#             yform - the format of the y axis labels - defaults to 0 dp
#
#-------------------------------------------------------------------------------
polygon.fit <-
function (dat, x=NULL,  mult=2, link="", line=TRUE)
{
  H <- if (link=="logit.percent") {
         function(x) ((1+exp(x))^-1)*100
       } else {
         function(x) x
       }
  if (!is.data.frame(dat)) dat <- as.data.frame(dat)
  n <- nrow(dat)
  dat$ciu <- dat$fit - mult*dat$se.fit
  dat$cil <- dat$fit + mult*dat$se.fit
  dat <- dat[,c("cil","fit","ciu")]
  dat <- H(dat)
  if (is.null(x)) x <- 1:n
  if (line==TRUE) {
     lines(x, dat$cil, lty=2)
     lines(x, dat$ciu, lty=2)
  } else {
    polygon(c(x, rev(x)), c(dat$cil, rev(dat$ciu)), col=grey(0.8), border=NA)
  }
  lines(x, dat$fit)
}


#-------------------------------------------------------------------------------
#
# function  : cov.glm
#
# purpose   : gets the covariance matrix of the estimated effects from a glm
#
# arguments : object - the glm object
#
#-------------------------------------------------------------------------------
cov.glm <-
function (object)
{
  df.r <- object$df.residual
  dispersion <-
    if (object$family$family %in% c("poisson", "binomial"))
    {
      1
    } else
    if (df.r > 0)
    {
      if (any(object$weights == 0))
      {
        warning("observations with zero weight not used for calculating dispersion")
      }
      sum((object$weights * object$residuals^2)[object$weights > 0])/df.r
    } else
    {
      NaN
    }
  p <- object$rank
  if (p > 0)
  {
    p1 <- 1:p
    Qr <- object$qr
    coef.p <- object$coefficients[Qr$pivot[p1]]
    covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
    dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
    covmat <- dispersion * covmat.unscaled
  } else
  {
    covmat <- matrix(,0,0)
  }
  return(covmat)
}


