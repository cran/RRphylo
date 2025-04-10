% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colorbar.R
\name{colorbar}
\alias{colorbar}
\title{Draw colorbar on a plot}
\usage{
colorbar(colors,x,y=NULL,direction="vertical",
  height=1,width=1,border="black",lwd=2,lty=1,
  labs=NULL,labs.pos=NULL,title=NULL,title.pos=NULL,
  ticks=TRUE,tck.pos=NULL,tck.length=1,xpd=FALSE,...)
}
\arguments{
\item{colors}{vector of colors.}

\item{x, y}{the x and y coordinates where the bottom left corner of the bar is
positioned. Keywords as in \code{\link[graphics]{legend}} are allowed.}

\item{direction}{either \code{"vertical"} or \code{"horizontal"}.}

\item{height}{a number indicating the amount by which the height of the bar
should be scaled relative to the default.}

\item{width}{a number indicating the amount by which the width of the bar
should be scaled relative to the default.}

\item{border}{color of the border around the bar. Set \code{NA} to suppress
border drawing.}

\item{lwd}{border line width.}

\item{lty}{border line type.}

\item{labs}{the vector of labels to place next to the bar.}

\item{labs.pos}{either \code{"left"}/\code{"right"} for
\code{direction="vertical"} or \code{"top"}/\code{"bottom"} for
\code{direction="horizontal"}. Default settings are \code{"right"} and
\code{"bottom"}.}

\item{title}{the title to be placed next to the bar.}

\item{title.pos}{either on the \code{"top"} or at the \code{"bottom"} of the
bar. Default setting is \code{"top"}.}

\item{ticks}{logical indicating whether ticks should be drawn next to each
label.}

\item{tck.pos}{indicates whether ticks should be plotter \code{"in"}side or
\code{"out"}side the bar border.}

\item{tck.length}{tick lengths}

\item{xpd}{a value of the \code{\link[graphics]{par}} \code{xpd}.}

\item{...}{further arguments passed to the functions \code{text} (for labels
and title) and \code{segments}. All these arguments must be hooked to the
element they refer to by indicating: \code{labs.}\* for labels,
\code{title.}\* for title, and \code{tck.}\* for ticks. See example for
further details.}
}
\description{
The function adds a color bar to the current plot.
}
\examples{

rainbow(30)->cols
replicate(4,paste(sample(letters,4),collapse=""))->labs

plot(rnorm(20),rnorm(20))
colorbar(cols,"topleft")

plot(rnorm(20),rnorm(20))
colorbar(cols,"topright",
         height=1.2,width=1.2,lwd=2,
         labs=labs,labs.pos="left",labs.cex=1.3,labs.adj=1,
         title="Colorbar!",title.cex=1.4,title.font=2,title.adj=c(0,0),
         tck.pos="out",tck.lwd=2,xpd=TRUE)
}
\author{
Silvia Castiglione
}
