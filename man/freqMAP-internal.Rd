\name{freqMAP-internal}
\alias{cat.moving.average}
\alias{cat.freq.post}
\alias{draw.signif.post.boxes}
\alias{add.p.value.legend}
\alias{rdirichlet}
\title{Internal freqMAP objects}
\description{Internal freqMAP objects.}
\details{
  Although these methods are not hidden, they should not generally be
  called directly by the user.
  \describe{
  \item{cat.moving.average}{Calculate the moving average estimates of observed
     category frequencies.}
  \item{cat.freq.post}{ Calculates marginal posterior intervals on
     the true category frequencies. }
  \item{draw.signif.post.boxes}{ Draws the Bayesian p-value bar on the
     current plot.}
  \item{add.p.value.legend}{ Adds a legend for the Bayesian p-value
     bar. }
  \item{rdirichlet}{ Generates random samples from the Dirichlet distribution.}
  }
}
\keyword{internal}
