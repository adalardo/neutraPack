\name{simula.neutra.hub}
\alias{simula.neutra.hub}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
simula.neutra.hub
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simula.neutra.hub(S = 100, j = 10, ciclo = 10000, step = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
%%     ~~Describe \code{S} here~~
}
  \item{j}{
%%     ~~Describe \code{j} here~~
}
  \item{ciclo}{
%%     ~~Describe \code{ciclo} here~~
}
  \item{step}{
%%     ~~Describe \code{step} here~~
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
function (S = 100, j = 10, ciclo = 10000, step = 100) 
{
    J <- S * j
    ind.mat = matrix(nrow = J, ncol = 1 + ciclo/step)
    n.dead <- c()
    n.dead[1] <- 0
    ind.mat[, 1] <- rep(1:S, each = j)
    cod.sp <- ind.mat[, 1]
    p.death <- 1/J
    for (i in 2:(1 + ciclo/step)) {
        n.mortes <- 0
        for (a in 1:step) {
            morte = rbinom(J, 1, prob = p.death)
            D = sum(morte)
            n.mortes <- n.mortes + D
            if (D > 0) {
                nascer = which(morte == 1)
                novos <- sample(1:J, D, replace = T)
                cod.sp[nascer] <- cod.sp[novos]
            }
        }
        ind.mat[, i] <- cod.sp
        n.dead[i] <- n.mortes
    }
    tempo <- seq(0, ciclo, by = step)
    colnames(ind.mat) <- tempo
    names(n.dead) <- tempo
    resulta = list(tempo = tempo, sp.list = ind.mat, n.mortes = n.dead)
    return(resulta)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
