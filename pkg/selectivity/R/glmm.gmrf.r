
print.glmmgmrf <-
function (x, ...)
{
  printf("Generalized linear mixed model fit using GMRF\n")
  printf("Formula: %s\n", deparse(x$formula))
  printf("  Data: %s\n", "data...")
  printf("  Random effects:\n")
  printf("  ...table...\n")
  print(apply(x $ samp.kappa, 2, quantile, c(0.025, 0.5, 0.975)))
  printf("Number of obs: %d\n\n", x$n)
  printf("  Fixed effects:\n")
  printf("  ...table...\n")
}


plot.glmmgmrf <-
function (x, type = c("histogram", "line", "acf"), nburn = 0, ...)
{
  type <- match.arg( type )
  ntot <- length( x $ samp.kappa[,1] )
  stopifnot( (ntot - nburn) > 99 )
  n <- seq.int(nburn + 1, ntot)
  par(mfrow = c(2,2))
  if (type == "histogram")
  {
    apply( x $ samp.kappa[n,], 2, hist, nclass = 50)
  } else
  if (type == "line")
  {
    apply( x $ samp.kappa[n,], 2, plot.ts)
  } else
  if (type == "acf")
  {
    apply( x $ samp.kappa[n,], 2, acf)
  }
  
  invisible( NULL )
}


mvrnormQ.chol <-
function (mu, L)
{
  # /*
  #  * Sample from a GMRF with mean mu and precision t(U) %*% U = Q
  #  */
  p <- length(mu)
  if ( L @ Dim [1] != p)
      stop("incompatible arguments")
      
  z <- rnorm(p)
  x <- drop(mu) + solve(L, z, system = "Lt")[order(L @ perm)]
  ldet <- determinant(L, logarithm = TRUE)
  ldet <- as.numeric(ldet $ sign * ldet $ modulus)
  ldens <- 0.5 * ldet - 0.5 * sum( z^2 )

  list(x = as.vector(x), ldens = ldens)
}


glmm.gmrf <-
function (M, control, family = c("gaussian","binomial"), ...)
{
  family <- match.arg(family)
  if (family == "gaussian")
  {
    return ( glmm.gmrf.gaussian (M, control) )
  } else
  if (family == "binomial")
  {
    return ( glmm.gmrf.gaussian (M, control) )
  }
}

O <- function (n,m) Matrix(0, n, m, sparse = TRUE)
I <- function (n) Diagonal(n, rep(1, n))

glmm.gmrf.gaussian <-
function ( M, control )
{
  printf("\n\tGMRF 'Gaussian' method - working on it ;)\n\n\tFast method ...\n")

  wk.dir <- set.working.directory()


  printf("\nPreparing Model Matrices ... "); tm <- cpu()
  n <- M $ n
  X <- Matrix(M $ X, sparse = TRUE)
  Z <- Matrix(M $ Z, sparse = TRUE)
  nf <- M $ nf
  nr <- sum(M $ nr)
  np <- n + nf + nr
  y <- c(unname(M $ Y), rep(0, nf + nr))

  # the prior
  
  A <- rBind(
         cBind(     I(n),        -X,      -Z  ),
         cBind( O(nf, n),     I(nf), O(nf,nr) ),
         cBind( O(nr, n), O(nr, nf),    I(nr) )
            )
  A <- unname(A)
  tA <- t(A)

  Q.func.val <- M $ Q.func.val 
  Q.struct <- M $ Q.struct
  environment(Q.func.val) <- environment(Q.struct) <- environment(NULL)
  Q.struct <- Q.struct()
  
  # the likelihood
  H.struct <- Diagonal( np, c( rep(1, n), rep(0, nf + nr) ) )
  llik_fun <- M $ llik
  llik_args <- M $ llik_args

  printf("done!\t(%.2f secs)\n",cpu() - tm)

  # try a sampling scheme

  res <- if (M $ family $ family == "gaussian") 1 else 0

  ntheta <- sum( M $ st) + res # set in glmm maybe
  if (!is.null(control $ theta.start))
  {
    if (length(control $ theta.start) != ntheta) stop("theta.start is wrong length")
    theta <- control $ theta.start
  } else
  {
    theta <- rep(4, ntheta)
  }
  theta.tune <- double(ntheta)
  suppressWarnings(theta.tune[] <- control $ theta.tune)

  nrho <- M $ nrho
  if (!is.null(control $ rho.start))
  {
    if (length(control $ rho.start) != nrho) stop("rho.start is wrong length")
    rho <- control $ rho.start
  } else
  {
    rho <- rep(0, nrho)
  }
  rho.tune <- double(nrho)
  suppressWarnings(rho.tune[] <- control $ rho.tune)
  
  fp.mon <- if (control $ verbose) "" else
              file(paste(wk.dir, "monitoring.txt", sep = .Platform$file.sep), "w")
  
  t1 <- t2 <- theta
  r1 <- r2 <- rho
  # to sample from the conditional posterior:

  sample.x <-
  function ( theta, rho )
  {
    #tm <- tm2 <- cpu()
    # fast ~ .00 secs
    Q.struct @ x <- Q.func.val(1e3, rep(1e-3, nf), theta[(res+1):ntheta], rho)
    Q <- tA %*% Q.struct %*% A
    
    # approximate likelihoood: THIS BIT IS SLOW!!
    #first find the mode
    
    fn <- function(x) M $ fast_llik (x, n, llik_args, tau_y = theta[1])
    gr <- function(x) M $ fast_gr (x, n, llik_args, tau_y = theta[1])
    #tm <- cpu()
    x.optim <- optim(rep(0, n), fn, gr, method = "BFGS", control = list(fnscale = -1))
    #cat("\n\tmode optim\t:", cpu() - tm)
    #stopifnot(x.optim $ convergence == 0) # not actually nessisary as only want an approximation.
    
    app <- approx_2order_fast( x.optim $ par, n, M $ fast_llik, M $ fast_gr, M $ fast_ggr, llik_args, tau_y = theta[1] )
    
    b <- double(np)
    b[1:n] <- app $ b
    H.struct @ x[1:n] <- app $ c
    Qp <- Q + H.struct

    # slow ~ 50% of function spent here, ~0.05 secs with n = 1000, h = 10 but can;t be made
    # any faster => fastest speed is 20 iter per sec
    L <- Cholesky(Qp, LDL = FALSE)

    # fast ~ 0.00
    mu <- solve(L, b, system = "A")

    # fast bits ~ 0.00
    x <- mvrnormQ.chol(mu, L)
    llik <- M $ fast_llik (x $ x, n, llik_args, tau_y = theta[1])

    # slow ~ 50% of function spent here ~ 0.05 wiht n = 1000
    ldet <- determinant(Q, logarithm = TRUE)
    ldet <- as.numeric(ldet $ sign * ldet $ modulus)
    lprior <- 0.5 * ldet - 0.5 * t(x $ x) %*% (Q %*% x $ x)
    #cat("\n\tTotal\t\t:", cpu() - tm2, "\n")
    list (x = x $ x, llik = llik, lprior = lprior, lprop = x $ ldens)
  }

  keep1 <- sample.x(t1, r1)

  timeref <- cpu()
  eprob <- 0.0
  samp <- matrix(NA, control $ niter %/% control $ thin , np + ntheta + nrho)
  j <- 0
  for (i in seq.int(control $ niter))
  { 
    if (ntheta)
        t2 <- t1 * sapply(theta.tune, scale_proposal)
    if (nrho)
        r2 <- r1 + sapply(rho.tune, function(sig) rnorm(1, 0, sig)) # perhaps a constant cv jump would be better...
    keep2 <- sample.x(t2, r2)
    #print(keep2[-1])

    logRx <- keep2$llik + keep2$lprior - keep2$lprop - keep1$llik - keep1$lprior + keep1$lprop
    logRt <- sum(dgamma(t2, 0.001, 0.001, log = TRUE) - dgamma(t1, 0.001, 0.001, log = TRUE))
    lacc <- logRx + logRt
    acc <- exp( min(lacc, 0.0) )

    if (runif(1) <= acc)
    {
      t1 <- t2
      r1 <- r2
      keep1 <- keep2
    }
    if (!(i %% control $ thin))
    {
      samp[j <- j + 1, ] <- c( keep1 $ x, t1, r1 )
    }
    
    # write some summary output
    eprob <- eprob + acc
    if (!(i %% 10))
    {
      fprintf(fp.mon, "niter: %d of %d | lacc= %.3f  E(accept_prob)= %.3f iter/sec= %.3f kappa= ( ",
                      i, control$niter, lacc, eprob / i, i / (cpu() - timeref) )
      fprintf(fp.mon, "%.4f ", t1)
      if (nrho)
      {
        fprintf(fp.mon, ") rho= ( ")
        fprintf(fp.mon, "%.3f ", tanh(r1) )
      }
      fprintf(fp.mon, ")\n")
      if (control $ verbose) flush.console()
    }
  }
  if (!control $ verbose) close(fp.mon)

  # write results
  
  if (.__write)
  {
    cat("\nwriting results ...\n")
    fp <- file(paste(wk.dir,"mcmc.out",sep = .Platform$file.sep), "w")
    fprintf(fp, "%15s ", colnames(samp)); fprintf(fp, "\n")
    for (i in seq.int(control $ niter))
    {
      fprintf(fp, "%15.8f ", samp[i,])
      fprintf(fp, "\n")
    }
    close(fp)
  }
  
  # accumulate output
  cat("\naccumulating results ...\n")
  out <- list()
  out$samp.kappa <- samp[, 1:ntheta + np]
  out$samp.rho <- if (nrho) samp[, 1:nrho + np + ntheta] else NULL
  out$samp.fixef <- samp[, 1:nf + n]
  out$samp.ranef <- samp[, 1:nr + n + nf]
  out$acc.prob <- eprob / control $ niter
  out$formula <- M$formula
  out$n <- n

  class(out) <- "glmmgmrf"

  cat("\ndone ...\n\n")

  return ( out )
}
