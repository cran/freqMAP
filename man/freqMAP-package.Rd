\name{freqMAP-package}
\alias{freqMAP-package}
\docType{package}
\title{
  Frequency Moving Average Plots (MAP) of Multinomial Data by a
  Continuous Covariate
}
\description{
  
 This package can be used to estimate a frequency moving average plots
 (MAP) from multinomial data and a continuous covariate. The frequency
 MAP is a moving average estimate of category frequencies, where
 frequency means and posterior bounds are estimated. Comparisons of two
 frequency MAPs as well as odds ratios can be plotted.
 
}
\details{
\tabular{ll}{
Package: \tab freqMAP\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2008-09-24\cr
License: \tab GPL (>= 2) \cr
}

This package allows one to visualize frequency data (e.g. genotype) as a
function of a continuous covariate (e.g. age). Two populations can be
visualized on the same plot and inference can be made about significant
differences in the frequencies and odds ratios comparing two unmatched
populations. This can be useful for exploring possible confounding
effects in a case control study. In this application, one frequency
MAP would be used for cases and another MAP for controls. Then
one would plot the frequency MAP of, say, genotype frequencies as a
function of age in both groups. Interactions between age and genotype
would show up as non-parallel frequency MAPs in the two groups. This
information can then be used for formulating interactions in the
subsequent logistic regression modeling.

At each of a pre-specified set of values of the continuous covariate,
the frequency data is modeled as Multinomial with an unknown true
category frequency vector. The prior on the true category frequency
vector is assumed to be uniform ( Dirichlet(1,...,1) ). Then samples are
generated from the Dirichlet posterior distribution on the true category
frequency vector. The central 95\% posterior interval on all true
category frequencies are estimated from the posterior samples.

The function \code{\link{freqMAP}} creates a \code{freqMAP} object and
samples from the posterior on true category frequencies. The object can
then be plotted using \code{\link{plot.freqMAP}}. The posterior
distributions from two unmatched populations contained in two
\code{freqMAP} objects can be compared using
\code{\link{posterior.comparison.freqMAP}}. \code{\link{plot.freqMAP}}
can also plot two populations at once as well as odds ratios comparing
the two populations.

The example below analyzes a typical dataset one would find in a genetic
association study, where genotypes are collected for a group of
probands and controls of varying ages. See \code{\link{genotypebyage}}
for more details on the dataset. The example code shows how to generate MAP
plots by genotype and by allele.

}
\author{
Colin McCulloch <colin.mcculloch@themccullochgroup.com>
}

\references{
  Payami, H., Kay, D.M., Zabetian, C.P., Schellenberg, G.D., Factor,
  S.A., and McCulloch, C.C. (2009) "Visualizing Disease Associations: Graphic 
  Analysis of Frequency Distributions as a Function of Age Using Moving 
  Average Plots (MAP) with Application to Alzheimer's and Parkinson's Disease",
  Genetic Epidemiology
}

\examples{

data("genotypebyage")

####An analysis by genotype

#Create a freqMAP object for probands and controls, respectively.
fm.cont <- freqMAP(genotypebyage[genotypebyage$subject.type=="control",c("age","genotype")],
                   x=seq(25,95,by=5),x.label="age",hw=5)
fm.prob <- freqMAP(genotypebyage[genotypebyage$subject.type=="proband",c("age","genotype")],
                   x=seq(25,95,by=5),x.label="age",hw=5)

#Plot the two MAPs overlaid. 
plot(fm.cont,fm.prob,legend=c("controls","probands"),show.p.value.legend=TRUE)


####An analysis by allele

#Reformat the dataset to stack the two alleles into one column. The
#resulting dataset has twice as many rows as genotypebyage
allelebyage <- rbind(cbind(genotypebyage[,c("age","subject.type")],
                           allele=genotypebyage$allele1,stringsAsFactors=FALSE),
                     cbind(genotypebyage[,c("age","subject.type")],
                           allele=genotypebyage$allele2,stringsAsFactors=FALSE))

#Create a freqMAP object for probands and controls, respectively.
fm.cont <- freqMAP(allelebyage[allelebyage$subject.type=="control",c("age","allele")],
                   x=seq(25,95,by=5),x.label="age",hw=5)
fm.prob <- freqMAP(allelebyage[allelebyage$subject.type=="proband",c("age","allele")],
                   x=seq(25,95,by=5),x.label="age",hw=5)

plot(fm.cont,fm.prob,legend=c("controls","probands"),show.p.value.legend=TRUE)


}

\keyword{ package }
\keyword{ htest }
\keyword{ models }
\keyword{ multivariate }
\keyword{ smooth }
