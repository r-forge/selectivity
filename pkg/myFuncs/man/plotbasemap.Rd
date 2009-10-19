\name{plotbasemap}
\Rdversion{1.1}
\alias{plotbasemap}
\title{
  Plot a map of Europe
}
\description{
  \code{plotbasemap} is a function
}
\usage{
  plotbasemap(lon1 = -4.5, lon2 = 9, lat1 = 50.5, lat2 = 62.5, new = TRUE,
              landcolor = grey(0.8), seacolor = NA, landborder = landcolor, data)
}
\arguments{
  \item{lon1}{~~Describe \code{formula} here~~}
  \item{lon2}{~~Describe \code{formula} here~~}
  \item{lat1}{~~Describe \code{formula} here~~}
  \item{lat2}{~~Describe \code{formula} here~~}
  \item{new}{~~Describe \code{formula} here~~}
  \item{landcolor}{~~Describe \code{formula} here~~}
  \item{seacolor}{~~Describe \code{formula} here~~}
  \item{landborder}{~~Describe \code{formula} here~~}
  \item{data}{~~Describe \code{formula} here~~}
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  \code{NULL} is returned invisibly
}
\references{
 ~put references to the literature/web site here ~
}
\author{
  Colin Millar
}

\seealso{
 ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
  data(gmt4)
  plotbasemap(data = gmt4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
