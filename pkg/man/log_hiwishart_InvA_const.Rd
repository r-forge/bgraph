\name{log_hiwishart_InvA_const}
\alias{log_hiwishart_InvA_const}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
normalizing constant for a HyperInv-Wishart
}
\description{
Generates the normalizing constant "cons."  for a HyperInv-Wishart
}
\usage{
log_hiwishart_InvA_const(G, df, S)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{
decomposed graph object
}
  \item{df}{
d.f.
}
  \item{S}{
covariance matrix
}
}
\details{
  Generates the normalizing constant "cons."  for a HyperInv-Wishart(G,df,S).
 Nonsingular  pdf is \eqn{p(clique) = cons. |K|^{-(df+2|p|)/2} exp(-trace(inv(K) S)/2) }
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