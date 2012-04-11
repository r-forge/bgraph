\name{log_GWishart_ud_const_mc}
\alias{log_GWishart_ud_const_mc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
log normalizing constant for a G-Wishart
}
\description{
log normalizing constant approximation using importance sampling for a G-Wishart(b,D,adj):
}
\usage{
log_GWishart_ud_const_mc(b, D, adj, N)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b}{
%%     ~~Describe \code{b} here~~
}
  \item{D}{
%%     ~~Describe \code{D} here~~
}
   \item{adj}{
adjacency matrix 
}
  \item{N}{
Monte Carlo sample size
}
}
\details{
log normalizing constant approximation using importance sampling for a G-Wishart(b,D,adj):

  \eqn{p(K) = I(b,D)^{-1} |K|^{(d-2)/2} exp(-trace(K D)/2)}

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