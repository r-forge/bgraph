\name{log_GWishart_NOij_pdf}
\alias{log_GWishart_NOij_pdf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Log density of G-Wishart upto the normalizing constant of 
}
\description{
Compute log p(Cc(i,j) ) upto the normalizing constant of G-Wishart
}
\usage{
log_GWishart_NOij_pdf(b_prior, D_prior, C, i, j, edgeij)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b_prior}{
%%     ~~Describe \code{b_prior} here~~
}
  \item{D_prior}{
%%     ~~Describe \code{D_prior} here~~
}
  \item{C}{
%%     ~~Describe \code{C} here~~
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
}
\details{
Compute \eqn{log p(C c(i,j) )} upto the normalizing constant of G-Wishart
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