\name{fert.t}
\alias{fert.t}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
fert.t
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fert.t1(cod.sp, n.propag, fun = mean)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cod.sp}{
%%     ~~Describe \code{cod.sp} here~~
}
  \item{n.propag}{
%%     ~~Describe \code{n.propag} here~~
}
  \item{fun}{
%%     ~~Describe \code{fun} here~~
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
function (cod.sp, n.propag, fun = mean) 
{
    especie <- unique(cod.sp[, 1])
    sp.level <- factor(cod.sp, levels = cod.sp)
    t.a <- function(x) {
        tapply(n.propag[, x], factor(cod.sp[, x], levels = especie), 
            fun)
    }
    res <- sapply(1:ncol(n.propag), t.a)
    colnames(res) <- colnames(n.propag)
    rownames(res) <- paste("sp", especie, sep = "")
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
