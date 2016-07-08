\name{simula.neutra}
\alias{simula.neutra}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
simula.neutra
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simula.neutra(S = 100, j = 10, xi0 = rep(seq(10,10,length.out = S),each=j), X = 10000, pi0 = NULL, dp = 500, dist.pos = NULL, dist.int = NULL, ciclo = 3e5, step = 100, trade.off = T, prop.neutro = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
%%     ~~Describe \code{S} here~~
}
  \item{j}{
%%     ~~Describe \code{j} here~~
}
  \item{xi0}{
%%     ~~Describe \code{xi0} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{pi0}{
%%     ~~Describe \code{pi0} here~~
}
  \item{dp}{
%%     ~~Describe \code{dp} here~~
}
  \item{dist.pos}{
%%     ~~Describe \code{dist.pos} here~~
}
  \item{dist.int}{
%%     ~~Describe \code{dist.int} here~~
}
  \item{ciclo}{
%%     ~~Describe \code{ciclo} here~~
}
  \item{step}{
%%     ~~Describe \code{step} here~~
}
  \item{trade.off}{
%%     ~~Describe \code{trade.off} here~~
}
  \item{prop.neutro}{
%%     ~~Describe \code{prop.neutro} here~~
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
function(S= 100, j=10, xi0=rep(seq(10,10,length.out = S),each=j), X=10000, pi0=NULL, dp=0.1, dist.pos=NULL, dist.int=NULL, ciclo=1e6, step=100, trade.off=T, prop.neutro=F)
{
  t0=proc.time()[[3]]
  J <- S*j
  if(max(xi0)>X)
  {
    stop("\n\tO número de propágulos produzidos por indivíduo a cada ciclo (xi0) não pode ser maior do que o número de propágulos máximo (X). Tente novamente!\n\n")
  }
  if(sum(dist.int>1)>0)
  {
    stop("\n\tA intensidade do distúrbio deve ser entre 0 e 1. Tente novamente!\n\n")
  } 
  if(sum(0<=dist.pos & dist.pos<=ciclo)<length(dist.pos))
  {
    stop("\n\tA posição dos eventos de distúrbio (dist.pos) precisa ser condizente com o número de ciclos a serem rodados (ciclo). Tente novamente!\n\n")
  } 
  if(sum(dist.pos==0 & dist.int>0)>0)
  {
    stop("\n\tAtenção! O primeiro ciclo das simulações é o ciclo 1. Atribua 0 ao dist.pos apenas quando não desejar implementar a ocorrência de distúrbios.\n\n")
  }
  if(step>ciclo)
  {
    stop("\n\tAtenção! O valor do argumento 'step' não pode ser maior que o de 'ciclo'.\n\n")
  }
  if(trade.off==T & prop.neutro==T)
  {
    stop("\n\tAtenção! Para que o caracter referente à produção de propágulos seja neutro, não pode haver trade-off entre este caracter e a probabilidade de morte, i.e., o argumento trade.off deve ser 'FALSE'.\n\n")
  }
  ind.mat=matrix(nrow=J,ncol=1+ciclo/step)
  prop.mat=matrix(nrow=J,ncol=1+ciclo/step)
  dead.mat=matrix(nrow=J,ncol=1+ciclo/step)
  n.dead <- 0
  n.dead.vetor<-c()
  n.dead.vetor[1] <- 0
  ind.mat[,1] <- rep(1:S,each=j)
  cod.sp <- ind.mat[,1]
  if(trade.off==T)
  {
    dead.mat[,1] <- xi0/X
  }
  if(trade.off==F)
  {
    dead.mat[,1] <- pi0
  }
  p.death <- dead.mat[,1]
  prop.mat[,1] <- xi0
  n.propag <- prop.mat[,1]
  prop.mat[,1]<-round(prop.mat[,1])
  sc=2
  for(i in 1:ciclo)
  {
    morte=rbinom(J, 1, prob=p.death)
    if(sum(dist.pos==i)>0)
    {
      vivos <- which(morte==0)
      nvivos <- length(vivos)
      if(length(dist.int)>1)
      {
        posdist <- which(dist.pos==i)
        ndist <- round(nvivos* dist.int[posdist])
      }
      if(length(dist.int)==1)
      {
        ndist <- round(nvivos* dist.int)
      }
      posmort <- sample(vivos, ndist)
      morte[posmort] <- 1
    }
    n.mortes <- sum(morte)
    n.dead <- n.dead + n.mortes
    if(n.mortes>0)
    {
      cod.ind<-1:J
      seed.bank <- rep(cod.ind,round(n.propag))
      nascer= which(morte==1)
      if(prop.neutro==F)
      {
        mami=sample(seed.bank, n.mortes)
      }
      if(prop.neutro==T)
      {
        mami=sample(cod.ind, n.mortes)
      }
      papi <- c()
      for(w in 1:n.mortes)
      {
        papi[w] <- sample(cod.ind[ cod.sp==cod.sp[mami[w]] ],1)
      }
      medias.prop=(n.propag[mami]+n.propag[papi])/2
      cod.sp[nascer]<-cod.sp[mami]
      n.propag[nascer] <- sapply(1,rtruncnorm,a=1, b=X , mean= medias.prop,sd=dp)
      if (trade.off==T)
      {
        p.death[nascer] <- n.propag[nascer]/X
      }
      if (trade.off==F)
      {
        medias.p.death <- (p.death[mami]+p.death[papi])/2
        p.death[nascer] <- sapply(1,rtruncnorm,a=1/200000, b=1, mean=medias.p.death,sd=dp)
      }
    }
    if(sum(i==seq(step,ciclo,step))==1)
    {
      ind.mat[,sc] <- cod.sp
      dead.mat[,sc] <- p.death
      prop.mat[,sc] <- round(n.propag)
      n.dead.vetor[sc] <- n.dead
      sc <- sc+1
    } 
  }
  tempo <- seq(0,ciclo,by=step)
  colnames(ind.mat) <- tempo
  colnames(dead.mat) <- tempo
  colnames(prop.mat) <- tempo
  names(n.dead.vetor) <- tempo
  resulta=list(tempo=tempo,sp.list=ind.mat,sementes=prop.mat,prob.morte=dead.mat,n.mortes.cumulativo=n.dead.vetor)
  t1=proc.time()[[3]]
  cat("Tempo de processamento: ", round((t1-t0)/60,2),"min\n")
  attributes(resulta)$start=list(especies=S, individuos=j, nprop=xi0, sd=dp, posicao_disturbios=dist.pos, intensidade_disturbios=dist.int, ciclos=ciclo, passos=step)
  return(resulta)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
