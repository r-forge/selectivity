taylor_2order <-
function (x0, indx, x_vec, loglFunc, loglFunc_arg)
#(double *a, double *b, double *c, double d, double x0, int indx,
#			  double *x_vec, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, double *step_len)
{
	#/*
	# * compute a,b,c in the taylor expansion around x0 of d*loglFunc(x0,...)
	# *
	# * a + b*(x-x0) + 0.5*c*(x-x0)^2
	# *
	# */

	xx <- double(3)
  xx[1] <- xx[2] <- xx[3] <- x0
  f <- loglFunc(xx, 3, indx, x_vec, loglFunc_arg)

	list (a = f[1], b = f[2], c = f[3])
}

approx_2order <-
function (x0, indx, x_vec, loglFunc, loglFunc_arg, ...)
#(double *a, double *b, double *c, double d, double x0, int indx,
#			  double *x_vec, GMRFLib_logl_tp * loglFunc, void *loglFunc_arg, double *step_len)
{
	#/*
	# * compute the second order approximation to the function around x0 of d*loglFunc(x0,...) as
	# *
	# * a + b*x -0.5*c*x^2
	# *
	# */

	xx <- double(3)
  xx[1] <- xx[2] <- xx[3] <- x0
  f <- loglFunc(xx, 3, indx, x_vec, loglFunc_arg, ...)
  df <- f[2]
	ddf <- f[3]

	a <- NULL # f[1] - df * x0 + 0.5 * ddf * x0^2
	b <- df - x0 * ddf
	c <- -ddf

  list (a = a, b = b, c = c)
}

approx_2order_fast <- 
function (xx, n, loglFunc, dloglFunc, ddloglFunc, loglFunc_arg, ...)
{
	#/*
	# * compute the second order approximation to the function around x0 of d*loglFunc(x0,...) as
	# *
	# * a + b*x -0.5*c*x^2
	# *
	# */

  f <- loglFunc(xx, n, loglFunc_arg, ...)
  df <- dloglFunc(xx, n, loglFunc_arg, ...)
	ddf <- ddloglFunc(xx, n, loglFunc_arg, ...)

	a <- NULL # f - df * xx + 0.5 * ddf * xx^2
	b <- df - xx * ddf
	c <- -ddf

  list (a = a, b = b, c = c)
}