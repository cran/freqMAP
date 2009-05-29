\name{plot.freqMAP}
\alias{plot.freqMAP}
\title{ Plotting Frequency Moving Average Plot (MAP) Objects }
\description{
  This function can plot a single frequency MAP object, can
  overlay two frequency MAP object plots, and can plot odds ratios (OR)
  comparing two frequency MAP objects.
}
\usage{
\method{plot}{freqMAP}(x, y = NULL, xlim = NULL, ylim = NULL, legend = NULL, show.p.value.legend = FALSE,
              type = "freq", p.value.bar.alpha = c(0.05, 0.01), cex = 1, layout.matrix = NULL,...)
}
\arguments{
  \item{x}{ A \code{freqMAP} object, usually the result of a call to
    \code{\link{freqMAP} }}
  \item{y}{ An optional \code{freqMAP} object, usually the result of a call to
    \code{\link{freqMAP} } }
  \item{xlim}{Limits on the x-axis. Must be a vector of length
    2. The same limits will be used for all category or OR plots. If not
  specified, then the function will calculate reasonable limits.}
  \item{ylim}{Limits on the y-axis. If \code{type="freq"},
    then \code{ylim} must be a matrix with 2 columns and as many rows as the length
    of \code{x$cat.names}. In this case, row i will be used as the y limits
    on the frequency plot of category i. If \code{type="or"}, then ylim
    must be a vector of length 2 specifying the common y limits for all
    odds ratio plots. If ylim is not specified, then the function will
    calculate reasonable limits which will be common on all plots.} 
  \item{legend}{A 2-vector of strings specifying legend text
    for the \code{x} and \code{y} \code{freqMAP} objects when
    \code{type="freq"}. This argument is ignored when \code{type="or"}
    or if \code{y} is not specified.  }
  \item{show.p.value.legend}{When set to \code{TRUE}, a legend is
    generated for the two levels of Bayesian p-values plotted. This
    argument is ignored when \code{y} is not specified.}
  \item{type}{ The type of plot to generate. Must be either
    \code{"freq"} or \code{"or"}. See \bold{Details} below. }
  \item{p.value.bar.alpha}{ A 2-vector of p-value cutoffs for
    the p-value bar at the bottom of plots. This argument is ignored if
    \code{y} is not specified. See \bold{Details} below. }
  \item{cex}{The size of all plotting characters, legend text, axis
    annotation and axis labeling. }
  \item{layout.matrix}{The matrix passed to \code{layout()} prior to
    plotting. When this argument is not used, the default behavior is to
    show all plots on the same figure stacked vertically.}
  \item{...}{Not used at this time.}
}
\details{
  
  Three types of plots are possible. First, if \code{y} is not
  specified, then a frequency MAP will be plotted for each category
  in the \code{x} object. The mean is plotted surrounded by 95\% central
  posterior intervals (CPI). Second, if \code{y} is specified and
  \code{type="freq"}, then a frequency MAP will be plotted for each
  category in the \code{x} and \code{y} objects, where the two frequency
  MAPs are overlaid on each category plot. Third, if \code{y} is
  specified and \code{type="or"}, then a plot is made for all odds
  ratios comparing each pair of categories. The posterior mean OR and
  95\% CPI on the OR are plotted. The title of each OR plot gives the
  direction of the OR. For example, if there are two categories "a" and
  "b", and the plot title says "b" vs. "a", this
  indicates the OR plotted is (t\_xb t\_ya) / (t\_xa t\_yb), where t\_xj and
  t\_yj indicate the true frequency of category j in the x and y
  MAPs, respectively.

  If \code{y} is specified, then  \code{x} and \code{y} must have the
  same values in their elements \code{[["cat.names"]]},
  \code{[["cat.short"]]}, and \code{[["hw"]]}. The
  \code{[["x.label"]]} element can be different. (This element
  indicates the name of the column in the \code{[["cat.ma"]]} element
  which contains the x values when plotting.)

  When \code{y} is specified, a gray bar is added to the bottom of the
  plot when the two frequency plots are "significantly"
  different. When \code{type="freq"}, significance is defined as the
  posterior Pr(category frequency in MAP \code{x} > category frequency
  in MAP \code{y}) > alpha/2
  or < (1-alpha/2). When \code{type="or"}, significance is defined as
  the posterior Pr( OR(category 1 vs. category 2) > 1 ) > alpha/2
  or < (1-alpha/2). The value of alpha is defined at two levels by the
  \code{p.value.bar.alpha} argument and plotted in light gray (higher
  alpha) and dark gray (lower alpha). This can be thought of as a
  Bayesian analog to a two-sided p-value. See
  \code{\link{posterior.comparison.freqMAP}} for further details.
  
  Specifying a legend or asking for a p-value legend, while not
  specifying \code{xlim}, will cause the x axis limits to be expanded
  20\% on the left side to give room for the legend. 
  
  Along with the example given below, see \code{\link{freqMAP-package}}
  for an example based on genotype data.
  
}
\author{ Colin McCulloch <colin.mcculloch@themccullochgroup.com> }

\references{
  Payami, H., Kay, D.M., Zabetian, C.P., Schellenberg, G.D., Factor,
  S.A., and McCulloch, C.C. (2009) "Visualizing Disease Associations: Graphic 
  Analysis of Frequency Distributions as a Function of Age Using Moving 
  Average Plots (MAP) with Application to Alzheimer's and Parkinson's Disease",
  Genetic Epidemiology
}

\seealso{ \code{\link{freqMAP}}, \code{\link{summary.freqMAP}},
  \code{\link{layout}}}

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
  pc <-  posterior.comparison.freqMAP(group1=fp1,group2=fp2)

  #Three example plots
  plot(fp1,ylim=matrix(c(0,1),nrow=2,ncol=2,byrow=TRUE))
  plot(fp1,fp2,type="freq",legend=c("y1","y2"),show.p.value.legend=TRUE)
  plot(fp1,fp2,type="or")

}
\keyword{ hplot }
