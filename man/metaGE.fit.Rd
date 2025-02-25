% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metaGEfit.R
\name{metaGE.fit}
\alias{metaGE.fit}
\title{Meta-analysis procedure: Fixed or Random effect.}
\usage{
metaGE.fit(Data, MatCorr, Method, NA.omit = TRUE, DropZScores = FALSE)
}
\arguments{
\item{Data}{A dataset containing the estimated marker effect and its associated pvalue of each marker (in rows) in each environment (in columns), as obtained from \code{\link[=metaGE.collect]{metaGE.collect()}}.}

\item{MatCorr}{The inter-environments correlation matrix. Can be computed using \code{\link[=metaGE.cor]{metaGE.cor()}}.}

\item{Method}{A string specifying the method to be performed: either "\code{Fe}" or "\code{Re}".}

\item{NA.omit}{A boolean specifying whether the markers with some \code{NA} values should be removed. (\code{TRUE} by default)}

\item{DropZScores}{A boolean specifying whether the Zscores should be dropped from the dataset or not.(\code{FALSE} by default)}
}
\value{
The dataset Data with supplementary columns:
\itemize{
\item PVALUE: The pvalue of the MA test,
\item Mu: Estimate of the mean marker effect,
\item Tau: Estimate of the variance of the marker effect, for the Random model only,
\item the Zscores for each environment if \code{DropZScores = FALSE}.
}
}
\description{
Quantitative trait loci detection via Fixed or Random effect meta-analysis GWAS procedure.
}
\details{
Different tests may be performed:
\itemize{
\item Fixed Effect (Fe), to identify markers with a stable effect across environments.
\item Random Effect (Re), to identify markers whose effects may be unstable across environments.
}
}
\examples{
require(dplyr)
# Import the data
data("metaData")

# Compute the inter-environment correlation matrix
matCorr <- metaGE.cor(metaData, Threshold = 0.8)

# Fixed Effect
FeDF <- metaGE.fit(metaData, matCorr, Method = "Fe")
head(FeDF \%>\% select(CHR, POS, MARKER, Mu, Tau, PVALUE))

# Random Effect
ReDF <- metaGE.fit(metaData, matCorr, Method = "Re")
head(ReDF \%>\% select(CHR, POS, MARKER, Mu, Tau, PVALUE))
}
