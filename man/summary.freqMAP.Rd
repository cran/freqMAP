\name{summary.freqMAP}
\alias{summary.freqMAP}
\title{ Summarizing Frequency MAP Objects }
\description{
  Summary method for freqMAP objects.
}
\usage{
\method{summary}{freqMAP}(object,...)
}
\arguments{
  \item{object}{ An object of class \code{freqMAP}, usually a result of a call
    to \code{\link{freqMAP}}}
  \item{...}{Not used at this time.} 
}
\value{
  \code{summary} returns the category moving average table of the
  object, element \code{[["cat.ma"]]}
}
\examples{
  #Make two sets of 2-category frequency data, y1 & y2, which both vary as
  #a function of a continuous variable x
  x <- runif(2000,min=-2,max=2)
  y1 <- c("a","b")[1+rbinom(n=length(x),size=rep(1,length(x)),prob=pnorm(x/2))]
  y2 <- c("a","b")[1+rbinom(n=length(x),size=rep(1,length(x)),prob=pnorm(x/5))]

  #Create the frequency MAP objects for y1 and y2
  fp1 <- freqMAP(data.frame(x=x,y=y1,stringsAsFactors=FALSE),
                  x=seq(-2,2,by=.2),x.label="x",hw=.2)
  fp2 <- freqMAP(data.frame(x=x,y=y2,stringsAsFactors=FALSE),
                  x=seq(-2,2,by=.2),x.label="x",hw=.2)

  #Examine the frequency MAP objects
  summary(fp1)
  print(fp2)

  #Three example plots
  plot(fp1,ylim=matrix(c(0,1),nrow=2,ncol=2,byrow=TRUE))
  plot(fp1,fp2,type="freq",legend=c("y1","y2"),show.p.value.legend=TRUE)
  plot(fp1,fp2,type="or")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ attribute }
