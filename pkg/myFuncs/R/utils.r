
.onAttach <-
function (...)
{
  print.myFuncs.version()
  set.myFuncs.options()
  invisible (NULL)
}

print.myFuncs.version <-
function ()
{
  version <- library(help = myFuncs)$info[[1]]
  version <- version[pmatch("Version", version)]
  um <- strsplit(version, " ")[[1]]
  version <- um[nchar(um) > 0][2]
  printf("\nThis is colins package of useful functions version %s.\nFor an overview type `help(\"myFuncs-package\")'.\n\n", version)
}

set.myFuncs.options <-
function ()
{
  pars <- list(las = 1, bty = "l", ann = FALSE, cex.axis = 0.7,
               mgp = c(2, 0.5, 0), tck = -.01, col.axis = grey(0.2),
               fg = grey(0.7))
  options(myFuncs.par = pars)
}

inv.logit <- function(x) 0.5 * (tanh(.5 * x) + 1)

logit <- function(p) 2 * atanh( 2*p - 1 )

ci <-
function(obj)
{
  obj $ fit + 2 * (-1:1) * sqrt(obj $ var)
}

zero.if.null <-
function (x)
{
  if (is.null(x)) return (0) else return(x)
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

eps <-
function(power) {
  # Return eps^power, where eps is the smalles number such that 1+eps != eps.

  eps <- 1.0
  val <- 2.0
  while (val > 1.0) {
    eps <- eps/2.0
    val <- 1 + eps
  }
  eps <- eps * 2.0
  return(eps^power)
}

clearfile <-
function (filename)
{
  if (file.exists(filename))
  {
    cat("\nover-writing '", filename, "'\n", sep="")
    unlink(filename)
  }
}

newline <-
function (filename)
{
  fprintf(filename, "\n")
}

getchar <-
function (obj, err.text = "input", use.name = FALSE)
{
  if( !missing(obj) && !is.null(obj) )
  {
 	  ischar <- tryCatch(is.character(obj) && length(obj) == 1L, error = identity)
 	  if(inherits(ischar, "error")) ischar <- FALSE

    if (ischar && !use.name)
    {
 	    if ( is.name(y <- substitute(obj)) ) warning(sprintf("%s is a character vector so its contents have been used", as.character(y)))
    }

   	## if this was not a length-one character vector, try for the name.
   	if(!ischar || use.name)
    {
      if( !is.name(substitute(obj)) )
          stop(sprintf("'%s' should be a name or a length-one character vector", err.text))
      obj <- deparse(substitute(obj))
   	}
  } else
  {
    obj <- NULL
  }
  obj
}
