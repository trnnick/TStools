\name{nemenyi}
\alias{nemenyi}

\title{Friedman and Nemenyi tests}

\description{Friedman and Nemenyi tests for nonparametric comparisons between methods.}

\usage{
nemenyi(data, conf.int=0.95, sort=c(TRUE,FALSE), 
        plottype=c("vline","none","mcb","vmcb","line"), 
        console=c(FALSE,TRUE), silent=c(FALSE,TRUE), 
        pcol=NULL, title=NULL, select=NULL, labels=NULL)
}

\arguments{
  \item{data}{
    A matrix that includes measures of accuracy for several methods (in columns) for several time series (rows), of size n x k.
}
  \item{conf.int}{
    The width of the confidence interval. Default is 0.95.
}
  \item{sort}{
    If TRUE function sorts the final values of mean ranks. If plots are request this is forced TRUE.
}
  \item{plottype}{
    Type of plot to produce:
1. "none": No plot; 2. "mcb": MCB style plot; 3. "vmcb": Vertical MCB style plot; 4. "line": Line visualisation (ISF style). Numbers next to method names are the mean rank; 5."vline": Vertical line visualisation [Default].
}
  \item{pcol}{
    Override default colours with single colour. Default is to use blue for MCB plots and rainbow colours for line plots. If no evidence of significant differences is found, based on Friedman's test, then all colours become grey.
}
  \item{console}{
    If FALSE, no console output is given.
}
  \item{silent}{
    If TRUE, console is overriden to FALSE and plottype to "none". Use to have only numerical output.
}
  \item{title}{
    Override default title.
}
  \item{select}{
    Highlight selected model (column). Number 1 to k. Select NULL for no highlighting.
}
  \item{labels}{
    Optional labels for models. If NULL column names of 'data' will be used. If number of labels < k, then it is assumed as NULL. 
}
}
\value{
\item{means}{Mean rank of each method.}
\item{intervals}{Nemenyi intervals for each method.}
\item{fpval}{Friedman test p-value.}
\item{cd}{Critical distance for nemenyi.}
}
\references{
Hibbon M., Crone S. F., Kourentzes N., 2012. Statistical significance of forecasting methods: An empirical evaluation of the robustness and interpretability of the MCB, ANOM and Nemenyi tests, The 32nd Annual International Symposium on Forecasting, Boston.
\url{http://www.forecasters.org/proceedings12/HibonMicheleISF2012.pdf}
}
\author{
Ivan Svetunkov and Nikolaos Kourentzes
}
\examples{
N <- 50
M <- 4
data <- matrix( rnorm(N*M,mean=0,sd=1), N, M) 
data[,2] <- data[,2]+1
data[,3] <- data[,3]+0.7
data[,4] <- data[,4]+0.5
colnames(data) <- c("Method A","Method B","Method C - long name","Method D");
out <- nemenyi(data,conf.int=0.95,plottype="vline")
}
\keyword{ ~friedman.test }