\name{GWishart_BIPS_pairwise}
\alias{GWishart_BIPS_pairwise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BIPS algorithm for sampling from G-Wishart distribution 
}
\description{
BIPS algorithm for sampling from G-Wishart distribution 
}
\usage{
GWishart_BIPS_pairwise(bG, DG, adj, C, burnin, nmc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bG}{
d.f.
}
  \item{DG}{
location
}
  \item{adj}{
adjacency matrix 
}
  \item{C}{
initial precision matrix
}
  \item{burnin}{
number of MCMC burnins
}
  \item{nmc}{
number of saved samples
}
}
\details{
     BIPS algorithm for sampling from G-Wishart distribution with  density:
          \eqn{p(C) \propto |C|^{(bG-2)/2} exp(-1/2 tr(C DG))}
     where     
      (1)  bG : d.f.
      (2)  DG: location
      (3)  adj: adjacency matrix 
     C: initial precision matrix;

}
\value{
\item{C }{Samples from G-Wishart distribution.
The result is 3D array with \code{c(dim(C)[1],dim(C)[2], nmc )}.
}
\item{Sig }{Inverse of C
The result is 3D array with \code{c(dim(C)[1],dim(C)[2], nmc )}.
}
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\references{Wang and Li (2011) "Efficient Gaussian Graphical Model Determination without 
      Approximating Normalizing Constant of the G-Wishart distribution "
\url{http://www.stat.sc.edu/~wang345/RESEARCH/GWishart/GWishart.html}
}
\author{Hao Wang \email{}; 
        Sophia Zhengzi Li \email{}
}
\examples{}
\keyword{models}
\keyword{methods}