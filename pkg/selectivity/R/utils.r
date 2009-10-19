mvrnormQ <-
function (mu, Q, tol = 1e-06)
{
  # /*
  #  * Sample from a GMRF with mean mu and precision Q
  #  */
  p <- length(mu)
  if (!all(dim(Q) == c(p, p)))
      stop("incompatible arguments")
      
  eQ <- eigen(Q, symmetric = TRUE, EISPACK = TRUE)

  ev <- eQ$values
  if (!all(ev >= -tol * abs(ev[1L])))
      stop("'Q' is not positive definite")
  X <- matrix(rnorm(p), 1)
  evinv <- 1/ev
  evinv[ev<tol] <- 0

  X <- drop(mu) + (eQ$vectors %*% diag(sqrt(evinv), p)) %*% t(X)
  ldet <- sum(log(ev[ev>0]))
  x <- X - mu
  ldens <- 0.5 * ldet - 0.5 * t(x) %*% (Q %*% x)
  list(x=X, ldens = ldens, ldet = ldet)
}


scale_proposal <-
function (A)
{
  # /*
  #  * realization of S where Pr(S=s) \propto 1/K, where K \in (1/A, A)
  #  */
	if (A <= 1.0) {
		return (1.0)
	} else {
		len <- A - 1 / A

		if (runif(1) < len / (len + 2 * log(A)))
    {
			return (1 / A + len * runif(1))
		} else
    {
			return (A^(2.0 * runif(1) - 1.0))
		}
	}
}

fprintf <-
function(fp, ...)
{
  txt <- sprintf(...)
  cat(txt, file=fp)
  invisible(TRUE)
}

printf <-
function(...)
{
  txt <- sprintf(...)
  cat(txt, file="")
  invisible(TRUE)
}

cpu <- function(gcFirst = FALSE)
{
  if (gcFirst) gc(FALSE)
  proc.time()[1L]
}

plot.cor <-
function(cor_mat)
{
  nc <- ncol(cor_mat)
  ylim <- c(1, nc + 1)
  plot(0, 0, ylim = ylim, xlim = ylim, ann = F, axes = F, type = "n")
  colour <- colorRampPalette(c("DarkBlue", "white", "Red"))(201)
  square <-
  function (i, j, col)
  {
    polygon(i + c(0, 1, 1, 0), j + c(0, 0, 1, 1), density = -1, col = col, border = F)
  }
  for (i in 1:nc)
  {
    for (j in 1:nc)
    {
      square(i, nc - j + 1, colour[101 + round(cor_mat[i,j] * 100)])
    }
  }
  invisible()
}

set.working.directory <-
function (keep = TRUE, working.directory)
{
  if (keep || missing(working.directory))
  {
    if (missing(working.directory))
    {
        working.directory <- wd.start <- "mcmc-run"
    } else
    {
      wd.start <- working.directory
    }

    kk <- 0
    working.directory <- paste(wd.start, "-", kk <- kk + 1, sep = "")
    while (file.exists(working.directory))
    {
      working.directory <- paste(wd.start, "-", kk <- kk + 1, sep = "")
    }
    wk.dir <- working.directory
    xx <- dir.create(wk.dir, showWarnings = FALSE)
    if (!xx)
    {
      stop(paste("\n\terror creating the directory ", inla.dir))
    } else
    {
      cat("\nCurrent working directory:", getwd(),"\n\n")
      cat("Model and results are stored in directory [",wk.dir, "]\n", sep = "")
    }
  } else
  {
    wk.dir <- tempfile()
    wk.dir <- gsub("\\\\", .Platform$file.sep, wk.dir)
    dir.create(wk.dir)
  }

  #
  # write a summary file giving the date and time etc of run - maybe even store required functions...
  #
  if (keep)
  {
    fp <- file(paste(wk.dir, .Platform$file.sep, "README.txt", sep=""), "w")
      fprintf(fp, "\nDirectory created on %s", date())
    close(fp)
  }

  return ( wk.dir )
}

