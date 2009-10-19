\name{gmt4}
\alias{gmt4}
\docType{data}
%# Example 7.2 page 228.
\title{ Bitterness of White Wines  }
\description{
   In a study on the bitterness of white wine it is of interest wheater treatments
   that can be controlleds during pressing the grapes influence the bitterness
   of wines.  The two factors considered are the temperature and the admission of
   contact with skin when pressing the grapes. 
}
\usage{data(gmt4)}
\format{
  A data frame with 72 observations on the following 5 variables.
  \describe{
    \item{lon}{a factor, temperature,  with levels \code{high} \code{low}}
    \item{lat}{a factor with levels \code{no} \code{yes}}
  }
}
\source{
     Ludwig Fahrmeir, Gerhard Tutz (1994):  Multivariate Statistical Modelling
   Based on Generalized Linear Models. Springer Series in Statistics. 
   Springer Verlag. New-York  Berlin Heidelberg
}
\examples{
data(gmt4)
str(gmt4)
summary(gmt4)
}
\keyword{datasets}
