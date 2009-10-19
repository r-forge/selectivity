
#
# simulate data
#

sim.example.data <-
function (n, seed, family = c("gaussian", "binomial"))
{
  family <- match.arg(family)
  
  if (!missing(seed))
      set.seed(seed)

  H <- 10

  beta0 <- 0
  beta1 <- 1
  beta <- matrix(c(beta0, beta1), nc=1)

  kappa_b0 <- 0.75
  kappa_b1 <- 1

  kappa_y <- 4

  b0 <- mvrnormQ(rep(0,H), diag(kappa_b0, H))$x
  b1 <- mvrnormQ(rep(0,H), diag(kappa_b1, H))$x
  b <- matrix(c(b0, b1), nc=1)

  dat <- data.frame(haul = sample(factor(1:H), n, replace=TRUE), len = runif(n, -20, 20))

  Z <- sapply(levels(dat$haul), function(i) as.numeric(dat$haul == i))
  Z <- cbind(Z, sweep(Z, 1, dat$len, "*"))

  X <- cbind(1, dat$len)

  nu <- X %*% beta + Z %*% b

  if (family == "gaussian")
  {
    dat$y <- mvrnormQ(nu, diag(kappa_y, n))$x
  } else
  if (family == "binomial")
  {
    dat$y <- as.numeric( runif(n) < exp(nu) / (1 + exp(nu)) )
  }
  assign("example.data", dat, "library:ColinsLib")

  invisible ( NULL )
}
