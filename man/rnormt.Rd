\name{rnormt}
\alias{rnormt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
rnormt
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rnormt(mean, n = 1, dp, min, max)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mean}{
%%     ~~Describe \code{mean} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{dp}{
%%     ~~Describe \code{dp} here~~
}
  \item{min}{
%%     ~~Describe \code{min} here~~
}
  \item{max}{
%%     ~~Describe \code{max} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (mean, n = 1, dp, min, max) 
{
    vals <- (min - 0.5):(max + 0.5)
    p <- pnorm(vals, mean = mean, sd = dp)
    p2 <- diff(p)
    sample(min:max, n, prob = p2, replace = T)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
