\name{GWishart_NOij_Gibbs}
\alias{GWishart_NOij_Gibbs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample from Gwishart distribution
}
\description{
Sample from Gwishart distribution
}
\usage{
GWishart_NOij_Gibbs(bG, DG, adj, C, i, j, edgeij, burnin, nmc)
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
initial partial covariance matrix
}
  \item{i}{
%%     ~~Describe \code{i} here~~
}
  \item{j}{
%%     ~~Describe \code{j} here~~
}
  \item{edgeij}{
%%     ~~Describe \code{edgeij} here~~
}
  \item{burnin}{
 number of MCMC burnins 
}
  \item{nmc}{
number of saved samples
}
}
\details{

     Sample C from Gwishart distribution with  density:
          \eqn{p(C) \propto |C|^{(bG-2)/2} exp(-1/2 tr(C DG))}
     where     
      (1)  bG : d.f.
      (2)  DG: location
      (3)  adj: adjacency matrix 
     C: initial partial covariance matrix;


}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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