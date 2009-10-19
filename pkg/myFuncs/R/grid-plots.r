#-------------------------------------------------------------------------------
#
# function  : grid.L.box
#
# purpose   : draws an L shaped box to fill out axes lines
#
# arguments : none
#
#-------------------------------------------------------------------------------
grid.L.box <-
function (x = unit(0.5, "npc"), y = unit(0.5, "npc"),
    default.units = "npc", name = NULL,
    gp = gpar(), draw = TRUE, vp = NULL)
{
    lg <- linesGrob(x = unit(c(0, 0, 1), "npc"),
                    y = unit(c(1, 0, 0), "npc"),
                    default.units = "npc",
                    arrow = NULL, name = name,
                    gp=gp, vp = vp)
    if (draw)
        grid.draw(lg)
    invisible(lg)
}

single.plot <-
function (x, y, y.ci = NULL,
          ylim = c(NA, NA), main = "",
          shade = FALSE,
          cex = 1, lty = NULL,
          xfmt = "%3.2f", yfmt = "%3.2f", lwd)
{

  # if y.ci is not present then multiple lines may be drawn

  xscale <- range(x, finite=TRUE) + diff(range(x, finite=TRUE))*0.04*c(-1,1)
  yscale <- if (is.null(y.ci))
    {
      ifelse(is.na(ylim), range(y, finite=TRUE), ylim)
    } else
    {
      ifelse(is.na(ylim), range(y, y.ci, finite=TRUE), ylim)
    }

  grid.text(main,
            x = 0.05, y = .975,
            default.units = "npc",
            hjust = 0, vjust = 1,
            gp= gpar(cex = cex))

  pushViewport(
    viewport(
      x = .10, y = .10,
      w = .80,  h = .80,
      just = c("left", "bottom"),
    )
  )

  pushViewport(
    viewport(
      xscale = xscale,
      yscale = yscale,
      clip = FALSE
    )
  )

  pushViewport(
    viewport(
      xscale = xscale,
      yscale = yscale,
      clip = TRUE
    )
  )

  if (!is.null(y.ci))
  {
    if (shade)
    {
      grid.polygon(c(x, rev(x)),
                   c(y.ci[,1],rev(y.ci[,2])),
                   default.units = "native",
                   gp = gpar(fill = grey(0.8), col = grey(0.8)))
    } else
    {
      grid.segments(x, y.ci[,1],
                    x, y.ci[,2],
                    default.units="native",
                    gp = gpar(lwd = cex),
                    arrow = arrow(angle = 90,
                                  ends = "both",
                                  length = unit(0.007,"npc")))
    }
  }
  if (is.null(ncol(y)))
  {
    if (is.null(lty)) lty <- 1
    grid.lines(x, y, default.units = "native", gp = gpar(lwd=2*cex, lty=lty))
  } else
  {
    if (is.null(lty))
    {
      lty <- rep(1, length.out = ncol(y))
    } else
    {
      lty <- rep(lty, length.out = ncol(y))
    }
    for (i in 1:ncol(y))
    {
      grid.lines(x, y[,i], default.units = "native", gp = gpar(lwd=lwd[i]*cex, lty=lty[i]))
    }
  }
  upViewport(1)

  grid.axis(1, fmt = xfmt, xscale, cex = cex)
  grid.axis(2, fmt = yfmt, yscale, dp=2, cex = cex)
  grid.L.box(gp = gpar(cex = cex))
  upViewport(2)
}

grid.axis <-
function (side = 1, fmt = "", axis.scale, tck = -0.01, mgp = -0.017, cex = 1, dp = NULL, col=1)
{
  at <- grid.pretty(axis.scale)
  if (!is.null(dp)) at <- round(at, dp)
  if (side==1)
  {
    grid.text(sprintf(fmt, at),
              x=unit(at, "native"),
              y=unit(mgp, "npc"),
              just="top",
              gp = gpar(cex = 0.75 * cex, col=col))
    grid.segments(unit(at, "native"), unit(0, "npc"),
                  unit(at, "native"), unit(tck, "npc"),
                  gp = gpar(lwd = cex, col=col))
  } else if (side == 2)
  {
    grid.text(sprintf(fmt, at),
              y=unit(at, "native"),
              x=unit(mgp, "npc"),
              just="right",
              gp = gpar(cex = 0.75 * cex, col=col))
    grid.segments(unit(0, "npc"), unit(at, "native"),
                  unit(tck, "npc"), unit(at, "native"),
                  gp = gpar(lwd = cex, col=col))
  }
}

new.plot <-
function ()
{
  grid.newpage()
  pushViewport(
    viewport(name = "base",
             x = .025, y = .025,
             w = .95,  h = .95,
             just = c("left", "bottom")) )
}

end.plot <-
function ()
{
  upViewport(1)
}

four.plot <-
function (dat, ...)
{

  pushViewport(
    viewport( layout = grid.layout(2, 2))
  )

  for (i in 1:2)
  {
    for (j in 1:2)
    {
      pushViewport(
        viewport(
          layout.pos.row = i,
          layout.pos.col = j,
          clip = TRUE
        )
      )
      pargs <- c(dat[(i-1)*2 + j], list(...))
      do.call("single.plot", pargs[[1]])
      upViewport(1)
    }
  }
  upViewport(1)
}
