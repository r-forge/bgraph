\name{GWishart_PAS_DMH}
\alias{GWishart_PAS_DMH}
\title{
Gaussian Non-decomposable graph determination using PAS algorithm
}
\description{
MCMC algorithm for Gaussian Non-decomposable graph determination using PAS algorithm.
}
\usage{
GWishart_PAS_DMH(b_prior, D_prior, n, S, C, beta, burnin, nmc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b_prior}{
d.f.
}
  \item{D_prior}{
location
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{S}{
%%     ~~Describe \code{S} here~~
}
  \item{C}{
initial partial covariance matrix
}
  \item{beta}{
%%     ~~Describe \code{beta} here~~
}
  \item{burnin}{
  number of MCMC burnin
}
  \item{nmc}{
  number of MCMC samples
}
}
\details{     
Sample C from Gwishart distribution with  density:
          \eqn{p(C) \propto |C|^{(bG-2)/2} exp(-1/2 tr(C DG))}
     where     
      (1)  bG : d.f.
      (2)  DG: location
     C: initial partial covariance matrix;

}
\value{
\item{C }{Samples from G-Wishart distribution.
The result is 3D array with \code{c(dim(C)[1],dim(C)[2], nmc )}.
}
\item{Sig }{Inverse of C
The result is 3D array with \code{c(dim(C)[1],dim(C)[2], nmc )}.
}
\item{adj}{ Sampled adjacency matrix.
The result is 3D array with \code{c(dim(C)[1],dim(C)[2], nmc )}.
}
}
\references{Wang and Li (2011) "Efficient Gaussian Graphical Model Determination without 
     Approximating Normalizing Constant of the G-Wishart distribution "
\url{http://www.stat.sc.edu/~wang345/RESEARCH/GWishart/GWishart.html}
}
\author{Hao Wang \email{}; 
        Sophia Zhengzi Li \email{}
        }
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
    p = 6; tedge = p*(p-1)/2; beta = 0.5; 
    indmx = matrix(1:(p^2),p,p); 
    b_prior = 3; 
    D_prior = diag(p); 
    n = 3*p; 
    b_post = b_prior+n;
    A = toeplitz(c(1,0.5,matrix(0,1,p-2))); 
    A[1,p] = 0.4; 
    A[p,1] = 0.4;
    S = solve(A)*n; 
    D_post = D_prior + S;
    adjTrue = 1*(abs(A)>0.001);  
    burnin  = 30; nmc = 100; C = diag(p); 
    resPAS_DMH <- GWishart_PAS_DMH(b_prior,D_prior,n,S,C,beta,burnin,nmc);
    C_save  <- resPAS_DMH[[1]]
    Sig_save<- resPAS_DMH[[2]]
    adj_save<- resPAS_DMH[[3]]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{methods}
