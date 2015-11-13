\name{graf.abund}
\alias{graf.abund}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
graf.abund

}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
graf.abund(dados = resulta, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dados}{
%%     ~~Describe \code{dados} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (dados = resulta, ...) 
{
    info = attributes(dados)
    n.prop = info$start$nprop
    cv = info$start$cv
    tempo = dados$tempo
    ntempo = length(tempo)
    nspp = length(unique(dados$sp.list[, 1]))
    nmax = max(table(dados$sp.list[, dim(dados$sp.list)[2]]))
    plot(x = c(1, nspp), y = c(1, nmax), log = "y", type = "n", 
        ylab = "Abundance", xlab = "Rank order", cex.lab = 1.2, 
        cex.axis = 1.2, sub = paste("total seeds=", n.prop, "\t cv = ", 
            cv), ...)
    stempos = round(seq(2, (ntempo - 100), length.out = 100))
    colors = rainbow(length(stempos) * 10)
    ncol = length(colors)
    ncol.meio = round(ncol * 0.5)
    i = ncol.meio
    for (t in stempos) {
        i = i + 1
        points(sort(table(factor(dados$sp.list[, t], levels = 1:nspp)), 
            decreasing = TRUE), type = "l", lwd = 0.5, col = colors[i])
    }
    points(sort(table(factor(dados$sp.list[, 1], levels = 1:nspp)), 
        decreasing = TRUE), type = "l", col = "green", lwd = 2)
    points(sort(table(factor(dados$sp.list[, ncol.meio], levels = 1:nspp)), 
        decreasing = TRUE), type = "l", col = "blue", lwd = 2)
    points(sort(table(factor(dados$sp.list[, dim(dados$sp.list)[2]], 
        levels = 1:nspp)), decreasing = TRUE), type = "p", col = "red", 
        lwd = 2)
    legend("topright", lty = 1, col = c("green", "blue", "red"), 
        bty = "n", legend = c("start", "middle", "end"))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
