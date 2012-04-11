\name{GWishart_BIPS_maximumClique}
\alias{GWishart_BIPS_maximumClique}

\title{
Sample from Gwishart distribution
}
\description{
Sample from Gwishart distribution
}
\usage{
GWishart_BIPS_maximumClique(bG, DG, adj, C, burnin, nmc)
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
  \item{burnin}{
number of MCMC burnins
}
  \item{nmc}{
number of MCMC saved samples
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
\examples{
  
##########################################################
##### e.g. (2) dense graph      ###########################
##########################################################
p     <- 10; 
alpha <- 0.3; 
J0    <- 0.5*diag(p); 
B0    <- 1*(matrix(runif(p*p),p,p)<alpha);
for( i in 1:p ){ 
    for( j in 1:i ){ 
        if ( B0[i,j] && i!=j ){
            J0[i,j] <- 0.5
        }
    }
}; 
J0    <- J0 + t(J0);
tmp   <- eigen( J0 ); w <- tmp$values;  
delta <- ( p * min( w ) - max( w ) ) / ( 1 - p ); 
J0    <- J0 + delta * diag( p );
D     <- diag(p) + 100*solve(J0);
DG    <- D;
adj   <- 1*(abs(J0)>1e-4); 
bG    <- 103;
    
burnin = 100; nmc = 100;

###  Maximum clique algorithm  ################
C = diag(p); # Initial value
resBIPSmax = GWishart_BIPS_maximumClique(bG,DG,adj,C,burnin,nmc);
C_save_maxC  <-resBIPSmax[[1]]
Sig_save_maxC<-resBIPSmax[[2]]

}
\keyword{models}
\keyword{methods}