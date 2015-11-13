\name{simula.neutra.step}
\alias{simula.neutra.step}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
simula.neutra.step
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simula.neutra.step(S = 100, j = 10, X = 1000, dp = 0.1, ciclo = 1e+06, step = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
%%     ~~Describe \code{S} here~~
}
  \item{j}{
%%     ~~Describe \code{j} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{dp}{
%%     ~~Describe \code{dp} here~~
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
function (S = 100, j = 10, X = 1000, dp = 0.1, ciclo = 1e+06, 
    step = 100) 
{
    t0 = proc.time()[[3]]
    J <- S * j
    if (abs(X/J - round(X/J)) > .Machine$double.eps^0.5) {
        stop("\n\tO potencial reprodutivo (X) precisa ser multiplo do tamanho da comunidade (J). Tente novamente!\n\n")
    }
    ind.mat = matrix(nrow = J, ncol = 1 + ciclo/step)
    prop.mat = matrix(nrow = J, ncol = 1 + ciclo/step)
    dead.mat = matrix(nrow = J, ncol = 1 + ciclo/step)
    n.dead <- c()
    n.dead[1] <- 0
    ind.mat[, 1] <- rep(1:S, each = j)
    cod.sp <- ind.mat[, 1]
    dead.mat[, 1] <- 1/J
    p.death <- dead.mat[, 1]
    prop.mat[, 1] <- X/J
    n.propag <- prop.mat[, 1]
    for (i in 2:(1 + ciclo/step)) {
        n.mortes <- 0
        for (a in 1:step) {
            morte = rbinom(J, 1, prob = p.death)
            D = sum(morte)
            n.mortes <- n.mortes + D
            if (D > 0) {
                seed.bank <- rep(1:J, n.propag)
                nascer = which(morte == 1)
                mami = sample(seed.bank, D)
                papi <- c()
                for (w in 1:D) {
                  papi[w] <- sample(n.propag[cod.sp == cod.sp[mami[w]]], 
                    1)
                }
                medias.prop = (n.propag[mami] + papi)/2
                cod.sp[nascer] <- cod.sp[mami]
                n.propag[nascer] <- sapply(medias.prop, rnormt, 
                  dp = dp, min = 1, max = X)
                p.death[nascer] <- n.propag[nascer]/X
            }
        }
        ind.mat[, i] <- cod.sp
        dead.mat[, i] <- p.death
        prop.mat[, i] <- n.propag
        n.dead[i] <- n.mortes
    }
    tempo <- seq(0, ciclo, by = step)
    colnames(ind.mat) <- tempo
    colnames(dead.mat) <- tempo
    colnames(prop.mat) <- tempo
    names(n.dead) <- tempo
    resulta = list(tempo = tempo, sp.list = ind.mat, sementes = prop.mat, 
        prob.morte = dead.mat, n.mortes = n.dead)
    t1 = proc.time()[[3]]
    cat("\n\t tempo de processamento: ", round((t1 - t0)/60, 
        2), "\n")
    attributes(resulta)$start = list(especies = S, individuos = j, 
        nprop = X, sd = dp, ciclos = ciclo, passos = step)
    return(resulta)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
