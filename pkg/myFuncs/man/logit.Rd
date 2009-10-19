\name{logit}
\Rdversion{1.1}

\alias{logit}
\alias{inv.logit}

\title{ Univariate Logistic and Inverse Logistic Transformation }
\description{
  Transforms a vector of numeric values from logistic space [0,1] to
  the real line (-Inf,Inf) ... and vice versa.
}
\usage{
  logit(p)
  inv.logit(x)
}
\arguments{
  \item{p}{a vector of probabilities}
  \item{x}{a vector probabilities on the logistic scale}
}
\details{
  Calculations are based on the \code{tanh} and \code{atanh} function.
}
\value{
  Either a vector of values on the real line for \code{logit} or a vector of
  probabilities for \code{inv.logit}
}

\seealso{
 \code{\link{atan}} \code{\link{atanh}}
}
\examples{
  (x <- logit(c(-Inf, 0.5, Inf)))
  inv.logit(x)
}
\keyword{ models }
