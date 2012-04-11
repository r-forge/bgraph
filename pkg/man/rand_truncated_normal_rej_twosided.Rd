\name{rand_truncated_normal_rej_twosided}
\alias{rand_truncated_normal_rej_twosided}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Truncated Normal
}
\description{
Generate one sample from \eqn{N(mu,sig^2)1_{above > x > below} }
}
\usage{
rand_truncated_normal_rej_twosided(mu, sig, below, above)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mu}{
mean
}
  \item{sig}{
standard deviation
}
  \item{below}{
lower cutoff
}
  \item{above}{
higher cutoff
}
}
\details{
Generate one sample from \eqn{N(mu,sig^2)1_{above > x > below}}
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

\references{Proposition 2.3  C. P. Robert 'Simulation of truncated normal variables'
see also Wang and Li (2011) "Efficient Gaussian Graphical Model Determination without 
      Approximating Normalizing Constant of the G-Wishart distribution "
\url{http://www.stat.sc.edu/~wang345/RESEARCH/GWishart/GWishart.html}
}
\author{Hao Wang \email{}; 
        Sophia Zhengzi Li \email{}
}
\examples{}
\keyword{models}
\keyword{methods}