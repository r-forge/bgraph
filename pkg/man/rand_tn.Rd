\name{rand_tn}
\alias{rand_tn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Truncated Normal
}
\description{
Sample from truncated normal \eqn{N(mu,sig^2)1_{lower<r<upper}}
}
\usage{
rand_tn(mu, sig, lower, upper)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mu}{
mean
}
  \item{sig}{
standard deviation
}
  \item{lower}{
lower cutoff
}
  \item{upper}{
upper cutoff
}
}
\details{
Sample from truncated normal \eqn{N(mu,sig^2)1_{lower<r<upper}}
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