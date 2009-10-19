.__f <-
function (..., model = c("iid"), extraconstr)
{

  vars <- as.list(substitute(list(...)))[-1]
  d <- length(vars)
  if (d == 0)
  {
    stop("At least one variable in f() needs to be defined")
  }
  term <- deparse(vars[[1]], backtick = TRUE, width.cutoff = 500)
  term <- attr(terms(reformulate(term)), "term.labels")
  if (d == 1)
  {
    weights <- NULL
    weights.intercept <- NULL
  } else
  {
    weights <- deparse(vars[[2]], backtick = TRUE, width.cutoff = 500)
    wt <- terms(reformulate(weights))
    weights <- attr(wt, "term.labels")
    weights.intercept <- attr(wt, "intercept")
    if (weights.intercept == 1) weights <- paste("1 +", weights)
  }
  if (d > 2)
  {
    stop("Only two variables can be passed to f():  f(term, weights)")
  }

#  inla.is.model(model, stop.on.error = TRUE)
  model <- match.arg(model)

  constr <- FALSE
  diagonal <- if (constr) 1e-6 else 0

  if (missing(extraconstr))
  {
    extraconstr <- NULL
  } else
  {
    A <- extraconstr$A
    e <- extraconstr$e
    if (!is.matrix(A))
    {
      stop("A(extraconstraint) has to be a matrix")
    } else
    if (nrow(A) != length(e))
    {
      stop("Dimension of A and e do not correspond")
    }
  }

  #prop <- inla.model.properties(model, stop.on.error = TRUE)
  prop <- list(ntheta = 1)

  if (prop$ntheta)
  {
    fixed <- rep(0, prop$ntheta)
  } else
  {
    fixed <- NULL
  }

  ret <- list(term = term, weights = weights, weights.intercept = weights.intercept,
              model = model, diagonal = diagonal, fixed = fixed, constr = constr,
              label = term, extraconstr = extraconstr)

  return(ret)
}

interpret.glmm <-
function (formula)
{
  env <- environment(formula)
  tf <- terms.formula(formula, specials = "f")

  terms <- attr(tf, "term.labels")
  nt <- length(terms)
  if (attr(tf, "response") > 0)
  {
    response <- as.character(attr(tf, "variables")[2])
    fixf <- randf <- weightf <- paste(response, "~")
  } else
  {
    stop("\n\tno response variable specified")
  }

  if (length(attr(tf, "offset")) > 0)
  {
    i <- grep("^offset[(]", attr(tf, "variables"))
    offset <- as.character(attr(tf, "variables")[i])
    if (.__debug)
      cat("found offset\n")
  } else
  {
    offset <- NULL
    if (.__debug)
      cat("no offset\n")
  }

  rt <- attr(tf, "specials")$f
  vtab <- attr(tf, "factors")
  if (length(rt) > 0)
  {
    for (i in 1:length(rt))
    {
      ind <- (1:nt)[as.logical(vtab[rt[i], ])]
      rt[i] <- ind
    }
  }

  k <- ks <- kp <- 1
  len.rt <- length(rt)
  random.spec <- list()
  if (nt > 0)
  {
    for (i in 1:nt)
    {
      if (k <= len.rt && ((ks <= len.rt && rt[ks] == i)))
      {
        f.call <- gsub("^f\\(", ".__f(", terms[i])
        st <- eval(parse(text = f.call), envir = env)
        random.spec[[k]] <- st
        if (ks <= len.rt && rt[ks] == i)
          ks <- ks + 1
        else kt <- kt + 1
        k <- k + 1
      } else
      {
        if (kp > 1)
        {
          fixf <- paste(fixf, " + ", terms[i], sep = "")
        } else
        {
          fixf <- paste(fixf, terms[i], sep = "")
        }
        kp <- kp + 1
      }
    }
  }

  intercept <- (attr(tf, "intercept") == 1)
  n.weights <- 0
  if (length(random.spec) > 0)
  {
    for (i in 1:length(random.spec))
    {
      ff1 <- random.spec[[i]]$term
      if (!is.null(random.spec[[i]]$weights))
      {
        ww1 <- random.spec[[i]]$weights
        n.weights <- n.weights + 1
      }

      if (is.null(random.spec[[i]]$weights))
      {
        randf <- paste(randf, if (i==1) " " else " + ",  ff1, sep = "")
      } else
      {
        randf <- paste(randf, if (i==1) " " else " + ", paste(ff1, ":", ww1, sep = ""), sep = "")
        weightf <- paste(weightf, if (i==1) " " else " + ", ww1, sep = "")
      }
    }
  }

  if (len.rt > 0)
  {
    randf <- as.formula(randf, env)
  } else
  {
    randf <- NULL
  }

  if (intercept)
  {
    fixf <- paste(fixf, "+ 1")
    nt <- nt + 1
  } else
  {
    fixf <- paste(fixf, "- 1")
    nt <- nt + 1
  }

  if ((nt - len.rt) > 0)
  {
    fixf <- as.formula(fixf, env)
  } else
  {
    fixf <- NULL
  }

  if (n.weights > 0)
  {
    weightf <- as.formula(weightf, env)
  } else
  {
    weightf <- NULL
  }

  ret <- list(randf = randf, random.spec = random.spec, n.random = len.rt,
              fixf = fixf, n.fix = (nt - len.rt), weightf = weightf,
              n.weights = n.weights, offset = offset, response = response)

  return ( ret )
}

glmm.control <-
function (niter = 100, thin = 1, theta.tune = 2, rho.tune = 1, verbose = FALSE, theta.start = NULL, rho.start = NULL)
{
  list(niter = niter, thin = thin,
       theta.tune = theta.tune, rho.tune = rho.tune,
       verbose = verbose,
       theta.start = theta.start, rho.start = rho.start)
}

glmm <-
function (formula, data, family = c("gaussian","binomial"), fit = TRUE, M = NULL, method = "default", control = glmm.control(...), ...)
{
  if (nargs() == 0)
  {
    cat("\tUsage: glmm(formula, data, family, other.arguments...); see ?glmm\n")
    return( invisible(NULL) )
  }

  if (is.null(M))
  {
    gp <- interpret.glmm(formula)
    call <- deparse(match.call())
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1, m)]
    if (gp$n.fix > 0)
    {
      new.fix.formula <- as.formula(paste("y ~", gp$fixf[3]))
      gp$model.matrix <- model.matrix(new.fix.formula, data = model.frame(new.fix.formula, data))
      gp$n.fix <- dim(gp$model.matrix)[2]
    } else
    {
      gp$model.matrix <- NULL
    }

    quantiles <- c(0.025, 0.975)
    family <- match.arg(family)

    if (gp$n.fix > 0)
    {
      mf$formula <- gp$fixf
    } else
    if (gp$n.random > 0)
    {
      mf$formula <- gp$randf
    } else
    {
      mf$formula <- y ~ 1
    }
    mf$na.action <- na.pass
    mf[[1]] <- as.name("model.frame")

    if (gp$n.random > 0)
    {
      rf <- mf
      rf$scale <- rf$Ntrials <- rf$offset <- rf$E <- NULL
      rf$formula <- gp$randf
      rf <- eval(rf, parent.frame())
    } else
    {
      rf <- NULL
    }
    if (gp$n.weights > 0)
    {
      wf <- mf
      wf$scale <- wf$Ntrials <- wf$offset <- wf$E <- NULL
      wf$formula <- gp$weightf
      wf <- eval(wf, parent.frame())
    } else
    {
      wf <- NULL
    }

    mf <- eval(mf, parent.frame())
    tot.dat <- length(mf[, 1])
    ind <- seq.int(tot.dat)
    if (!is.null(gp$offset))
    {
      off <- 0
      for (i in seq(along=gp$offset))
      {
        off <- off + as.vector(eval(parse(text = gp$offset[i]), data))
      }
    } else
    {
      off <- NULL
    }

    nr <- gp$n.random
    n.weights <- 0
    j <- 0
    extra.fixed <- 0
    if (nr > 0)
    {
      #if (nr != (ncol(rf) - 1))
      #    stop("\n\tSOMETHING STRANGE in nr...")
      location <- covariate <- list()
      count.linear <- count.random <- 0
      for (i in 1:nr)
      {
        count.random <- count.random + 1
        xx <- rf[, gp$random.spec[[i]]$term ]
        if (is.factor(xx))
        {
          location[[i]] <- sort(unique(xx))
          cov <- match(xx, location[[i]])
          cov[is.na(cov)] <- -1
          covariate[[i]] <- cov
        } else
        {
          if (!is.null(gp$random.spec[[i]]$values))
          {
            location[[i]] <- sort(unique(gp$random.spec[[i]]$values))
            cov <- match(xx, location[[i]])
            cov[is.na(cov)] <- -1
            covariate[[i]] <- cov
          } else
          {
            location[[i]] <- sort(unique(xx))
            cov <- match(xx, location[[i]])
            cov[is.na(cov)] <- -1
            covariate[[i]] <- cov
          }
        }

        n <- length(location[[i]])

        if (!is.null(gp$random.spec[[i]]$extraconstr))
        {
          A <- gp$random.spec[[i]]$extraconstr$A
          e <- gp$random.spec[[i]]$extraconstr$e
          if (ncol(A) != n)
          {
            stop("\n\tNcol in matrix A(extraconstraint) does not correspont to the length of f")
          }
        }

        if (!is.null(gp$random.spec[[i]]$weights))
        {
          www <- wf[, n.weights + 2]
          if (sum(is.na(www)) != 0) www[is.na(www)] <- -1
          n.weights <- n.weights + 1
        }
      }
    }

    num.threads <- 2

    #
    # now we have done things the inla way, make variables for a glmm fit
    #
    M <- glmm.setup()
    M$num.threads <- num.threads
    M$gp <- gp
    M$family <- family
    M$mterms <- terms(gp$fixf)
    M$rterms <- terms(gp$randf)
    M$mf <- mf
    M$rf <- rf
    M$wf <- wf
    M$offset <- off
    M$formula <- formula
    environment(M$formula) <- environment(formula)
    M$covariate <- covariate
    M$location <- location
    M$quantiles <- quantiles
    M$X <- model.matrix(M$mterms, M$mf)
    M$Z <- model.matrix.rand(M$rterms, M$rf)
    M$n <- tot.dat
    M$index <- ind
  }
  if (!fit) return ( M )

  #
  # now we have all we need to run the MCMC...
  #

  ret <- estimate.glmm(M, method, control, family)

  return ( ret )
}

glmm.setup <-
function (...)
{
  return ( list() )
}

model.matrix.rand <-
function (object, data = environment(object))
{
  t <- if (missing(data))
      terms(object)
  else terms(object, data = data)
  if (is.null(attr(data, "terms")))
  {
    data <- model.frame(object, data, xlev = xlev)
  } else
  {
    reorder <- match(sapply(attr(t, "variables"), deparse, width.cutoff = 500)[-1L], names(data))
    if (any(is.na(reorder)))
        stop("model frame and formula mismatch in model.matrix()")
    if (!identical(reorder, seq_len(ncol(data))))
        data <- data[, reorder, drop = FALSE]
  }
  int <- attr(t, "response")
  if (length(data))
  {
    namD <- names(data)
    for (i in namD)
    {
      if (is.character(data[[i]]))
      {
        data[[i]] <- factor(data[[i]])
        warning(gettextf("variable '%s' converted to a factor", i), domain = NA)
      }
    }
  } else
  {
    data <- list(x = rep(0, nrow(data)))
  }

  ans <- .Internal(model.matrix(t, data))
  attr(ans, "assign") <- NULL
  ans <- ans[,-1]
  
  # 'fix' model matrix so that there are no constraints...
  
  # ... hmmm not sure how to do this yet - will do a very dirty fix in the main code ! :(

  return ( ans )
}


estimate.glmm <-
function ( Model.spec, method = c("gmrf","default"), control, family )
{
  #
  # so far emtpy...
  #
  method <- match.arg(method)
  out <-
    switch(method,
      gmrf = glmm.gmrf( Model.spec, control, family ),
      default = glmm.default ( Model.spec, control, family )
    )

  return ( out )
}


glmm.default <-
function ( M, control, family )
{
  printf("\n\tdefault method - not implemeted yet!\n")
  return ( NULL )
}
