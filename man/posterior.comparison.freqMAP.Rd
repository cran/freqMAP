\name{posterior.comparison.freqMAP}
\alias{posterior.comparison.freqMAP}
\title{ Compare Two Frequency Moving Average Plot (MAP) Objects }
\description{
  This function compares the posterior distributions on two frequency
  MAPs from two unmatched populations.  Posterior summaries are
  calculated of the population frequency differences and
  log odds ratios.
}
\usage{
posterior.comparison.freqMAP(group1, group2)
}
\arguments{
  \item{group1}{ A \code{freqMAP} object }
  \item{group2}{ A \code{freqMAP} object }
}
\details{

  Posterior summaries are calculated by simply comparing the posterior
  samples inside the two \code{freqMAP} objects. The probability that
  group1 frequency is greater than group2 frequency is calculated, the
  log odds ratio mean and 95\% central posterior interval (CPI) is
  also calculated for each pair of categories.
  
  For any given category, the posterior distributions on the frequencies
  in the two populations are assumed independent. Namely the category
  data in the two populations is assumed to be unmatched.

  The two frequency MAPs \code{group1} and \code{group2} must have the
  same values in their elements
  \code{[["cat.names"]]},\code{[["hw"]]}, and
  \code{[["x.label"]]}. Further, the values of the continuous covariate in
  \code{...$cat.ma[,...$x.label]} must be the same in both groups.
  
  
}
\value{
  Dataframe with the following columns:
  \item{ x }{ The values of the continuous covariate. }
  \item{ n1 }{ The number of observations in group 1 falling into the
    bucket centered at \code{x[i]}.}
  \item{ n2 }{ The number of observations in group 2 falling into the
    bucket centered at \code{x[i]}.}
  \item{ *.gr1.gt.gr2 }{ (One column for each category) The posterior
    probability that group 1 frequency for this category is greater than
    the frequency in group 2. }
  \item{ *.*.lor.mean }{ (One column for each possible pair of
    categories) The posterior mean on the log odds ratio for each pair
    of categories. }
  \item{ *.*.lor.lpi & *.*.lor.cpi }{(One column for each possible pair of
    categories) The posterior lower (lpi) and upper bounds (upi) on the
    CPI on the log odds ratio for each pair of categories. } 
}
\author{ Colin McCulloch <colin.mcculloch@themccullochgroup.com> }
\seealso{ \code{\link{freqMAP}}, \code{\link{plot.freqMAP}} }
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

  #Compare the posterior distributions on the two frequency MAPs
  pc <- posterior.comparison.freqMAP(group1=fp1,group2=fp2)

  #Three example plots
  plot(fp1,ylim=matrix(c(0,1),nrow=length(fp1$cat.names),ncol=2,byrow=TRUE))
  plot(fp1,fp2,type="freq",legend=c("y1","y2"),show.p.value.legend=TRUE)
  plot(fp1,fp2,type="or")

}
\keyword{ htest }
\keyword{ models }
