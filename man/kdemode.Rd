\name{kdemode}
\alias{kdemode}

\title{Mode estimation via KDE}

\description{Return mode of a vector as calculated using KDE.}

\usage{
kdemode(data,type=c("diffusion","SJ","nrd0"),outplot=c(FALSE,TRUE))
}

\arguments{
  \item{data}{
    One-dimensional vector of data.
}
  \item{type}{
    Bandwidth selection: 1. "diffusion": Kernel Density Estimation Via Diffusion, Botev el al. 2010; 2. SJ: Sheater and Jones method; 3. nrd0: Silverman heuristic.
}
  \item{outplot}{
    If TRUE provides plot of the KDE and the mean, median and mode.
}
  \item{...}{
    Additional arguments can be passed to plotting functions. For example use main="" to replace the title.
}
}
\value{
\item{mode}{Estimated mode value.}
\item{xd}{Estimated values.}
\item{fd}{Estimated density.}
}
\references{
Kourentzes, Nikolaos, Devon K. Barrow, and Sven F. Crone. "Neural network ensemble operators for time series forecasting." Expert Systems with Applications 41.9 (2014): 4235-4244.

Botev, Z. I., J. F. Grotowski, and D. P. Kroese. "Kernel density estimation via diffusion." The Annals of Statistics 38.5 (2010): 2916-2957.
}
\author{
Nikolaos Kourentzes
}
\examples{
data <- rnorm(200,mean=0,sd=1)
kdemode(data,outplot=TRUE)
}