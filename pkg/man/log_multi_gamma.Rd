\name{log_multi_gamma}
\alias{log_multi_gamma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
multivariate gamma function at value n,p
}
\description{
multivariate gamma function at value n,p
}
\usage{
log_multi_gamma(p, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
}
\details{
  MULTI_GAMMA returns multivariate gamma function at value n,p 
   \deqn{\Gamma_p(n/2)= \pi^{p(p-1)/4}\Pi_{j=1}^p \Gamma\left[ (n+1-j)/2\right]}. 
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