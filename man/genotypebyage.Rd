\name{genotypebyage}
\alias{genotypebyage}
\docType{data}
\title{ Simulated Genotype Data by Age for 4000 Probands and 4000 Controls  }
\description{
  A simple simulated dataset comprising 4000 probands and 4000 controls, where
  subject age is tabulated against their genotype. 
}
\usage{data(genotypebyage)}
\format{
  A data frame with 8000 observations on the following 5 variables.
  \describe{
    \item{\code{age}}{Subject age (a numeric vector)}
    \item{\code{subject.type}}{Subject type: proband or control (a character vector)}
    \item{\code{allele1}}{Subject's first allele: "a" or "b" (a character vector)}
    \item{\code{allele2}}{Subject's second allele: "a" or "b" (a character vector)}
    \item{\code{genotype}}{Subject's genotype, coded "allele1.allele2" (a character vector)}
  }
}
\details{
  
  This dataset is used for the examples in the main help page for the
  freqMAP: \code{\link{freqMAP-package}}. See there for usage
  details. In MAP plots as a function of age, the genotype frequency
  curves in this dataset demonstrate non-homogeneous frequency
  patterns. For instance, at ages above 70, the "b.b" genotype becomes
  more frequent in probands than in younger probands. To compensate, the
  "a.a" genotype is less frequent in older probands than younger
  ones. The "a.b" genotype in probands is more homogeneous over all
  ages. Control genotype frequencies are also relatively homogeneous
  over all ages.
  
}
\references{
  Payami, H., Kay, D.M., Zabetian, C.P., Schellenberg, G.D., Factor,
  S.A., and McCulloch, C.C. (2009) "Visualizing Disease Associations: Graphic 
  Analysis of Frequency Distributions as a Function of Age Using Moving 
  Average Plots (MAP) with Application to Alzheimer's and Parkinson's Disease",
  Genetic Epidemiology
}
\examples{
data(genotypebyage)
}
\keyword{datasets}
