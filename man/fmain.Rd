% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fmain.R
\name{fmain}
\alias{fmain}
\title{Run HSA}
\usage{
fmain(lsmap0, lscov0, outfile, Maxiter, submaxiter, lambda, Leapfrog, epslon,
  mkfix = 0, rho = 0, mk, initialS = NULL, coarsefit = TRUE,
  rmoutlier = FALSE, fitmode = 0)
}
\arguments{
\item{lsmap0}{To be documented.}

\item{lscov0}{To be documented.}

\item{outfile}{To be documented.}

\item{Maxiter}{To be documented.}

\item{submaxiter}{To be documented.}

\item{lambda}{To be documented.}

\item{Leapfrog}{To be documented.}

\item{epslon}{To be documented.}

\item{mkfix}{To be documented.}

\item{rho}{To be documented.}

\item{mk}{To be documented.}

\item{initialS}{To be documented.}

\item{coarsefit}{To be documented.}

\item{rmoutlier}{To be documented.}

\item{fitmode}{To be documented.}
}
\value{
A named list with elements \code{Ao}, \code{bo}, \code{invSigmao}, \code{Po},
and \code{betao}.
}
\description{
This \code{fmain()} function is an tidied-up and packaged version of the one
available in the R scripts provided by Zuo et al. (2016) available from
the authors' website \url{http://ouyanglab.jax.org/hsa/}.
For details on usage, please see the instructions provided by the authors
at their website.
}
\references{
Zou C, Zhang Y, Ouyang Z. (2016) HSA: integrating multi-track Hi-C data
for genome-scale reconstruction of 3D chromatin structure.
Genome Biology, 17:40.
\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4774023/}{PMC4774023}
}
\seealso{
\url{http://ouyanglab.jax.org/hsa/}
}
