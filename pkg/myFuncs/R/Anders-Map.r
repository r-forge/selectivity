
#lon1 <- -4.5
#lon2 <- 9
#lat1 <- 50.5
#lat2 <- 62.5
#new <- TRUE
#plot.it <- TRUE
#landcolor <- grey(0.8)
#seacolor <- "blue"
#landborder <- landcolor
#data <- NULL

plotbasemap <-
function (lon1 = -4.5, lon2 = 9, lat1 = 50.5, lat2 = 62.5, new = TRUE,
          landcolor = grey(0.8), seacolor = NA, landborder = landcolor,
          data = NULL)
{
  if (is.null(data))
  {
    data(gmt4)
    data <- gmt4
  }

  if (new) plot.new()

  xrange <- c(lon1, lon2)
  yrange <- c(lat1, lat2)

  aspect <- c(cos((mean(yrange) * pi)/180), 0.8)
  d <- 2 * c(diff(xrange), diff(yrange)) * aspect
  p <- par("fin") - drop(matrix(c(0, 1, 1, 0, 0, 1, 1, 0), 2) %*% par("mai"))
  p <- d * min(p/d)
  
  d <- (p / min(p/d) - d) / 2 / aspect
  realusr <- c(xrange, yrange) + rep(c(-1, 1), 2) * rep(d, c(2, 2))

  par(pin = p, usr = realusr)
  rect(lon1, lat1, lon2, lat2, col = seacolor)
#  if (xrange[1] < 0)
#  {
#    par(usr = realusr + c(360, 360, 0, 0))
#    polygon(data, border = landborder, col = landcolor)
#  }
#  if (xrange[2] > 360)
#  {
#    par(usr = realusr - c(360, 360, 0, 0))
#    polygon(data, border = landborder, col = landcolor)
#  }
#  par(usr = realusr)
  polygon(data, border = landborder, col = landcolor)
  rect(lon1, lat1, lon2, lat2, lwd = 1)
  
  invisible (NULL)
}

bp <-
function (x, y, v, scale = 3, ...)
{
  val <- sqrt(abs(v))
  points(x, y, cex = val * scale, col = "darkblue", pch=16, ...)
  points(x, y, cex = val * scale, col = "white", pch=1, lwd=1.5, ...)
  
  invisible (NULL)
}
