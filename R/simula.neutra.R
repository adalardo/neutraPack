##############################################################################
######### Evolucao de estrategias de vida a partir de cenario neutro #########
##############################################################################
############################# simula.neutra.step #############################
##############################################################################
############################ Listagem das versoes ############################
##############################################################################
# 4. Versao modificada por Luisa Novara e Alexandre Adalardo (fevereiro 2016)
# 3. Disturbio implementado por Alexandre Adalardo e Luisa Novara (novembro 2015) 
# 2. Versao modificada por Luisa Novara (2014)
# 1. Versao inicial de Alexandre Adalardo (outubro 2009) modificada por Paulo Inacio Prado (outubro 2009)
##############################################################################
########################## Detalhamento das versoes ##########################
##############################################################################
# 4. Troca o argumento X por xi0, para que X (que eh dado por xi0*J) seja inteiro
# 2. Troca o antigo argumento cv (coeficiente de variacao) pelo dp (desvio padrao) da distribuicao normal da herdabilidade de xi. Isso permite que a herdabilidade permaneça constante e evita o erro de gerar herdabilidade maior quando os valores de xi sao mais baixos (i.e., evita que os valores de xi da populacao fiquem "presos" em valores mais baixos e isso mascare os reais resultados das simulacoes)
##############################################################################
############################ Argumentos da funcao ############################
##############################################################################
#      S =  número de especies da comunidade
#      j = número inicial de individuos por especie
#      xi0 =  número de propagulos que o individuo i produz por ciclo no inicio da simulacao
#      dp =  desvio padrao da distribuicao normal da herdabilidade de xi
#      dist.pos = vetor com a identificacao dos ciclos em que ocorrem eventos de disturbio
#      dist.int = vetor com a intensidade de cada evento de disturbio, se ha soh um valor este eh aplicado para todos
#      ciclo = numero de ciclos da simulacao
#      step = intervalo de ciclos a cada qual sao salvos resultados da simulacao
##############################################################################
################################# Deduzidos ##################################
############################################################################## 
#      Dado que o esforco reprodutivo total eh o mesmo para todos os individuos (X), o tempo medio de vida de cada individuo eh E[ti] = X/xi
#      A cada intervalo, cada individuo tem uma probabilidade fixa pi de morrer. Portanto, a probabilidade de sobreviver t intervalos eh dada por uma distribuicao geometrica com suporte 0,8. A esperanca desta distribuicao e E[ti]=pi-1, portanto:
#            pi = xi/X
#      Dado s igual para todas as especies, e uma distribuicao gaussiana de caracteres continuos, o numero de propagulos por intervalo (xi) da prole de um individuo pode ser definido como uma variavel normal com parametros:
#            µ = xi
#            sd = dp
##############################################################################
############################## Inicio da funcao ##############################
##############################################################################
simula.neutra.step=function(S= 100, j=10, xi0=10, dp=0.1, dist.pos=NULL, dist.int=NULL, ciclo=1e6, step=100)
  {
  t0=proc.time()[[3]] ### Marca o inicio da contagem de tempo de processamento da funcao
  #############################################################################
  ########################### Argumentos deduzidos ############################
  #############################################################################
  J <- S*j ### Calcula o tamanho da comunidade (J)
  X <- xi0*J ### Calcula o numero total de propagulos produzidos por um individuo (X)
  #############################################################################
  ############################### Verificacoes ################################
  #############################################################################
  if(sum(0<=dist.pos & dist.pos<=ciclo)<length(dist.pos)) ### Verifica se ha disturbio programado para ciclo inexistente
  {
    stop("\n\tA posição dos eventos de distúrbio (dist.pos) precisa ser condizente com o número de ciclos a serem rodados (ciclo). Tente novamente!\n\n")
  } 
  if(sum(dist.pos==0 & dist.int>0)>0) ### Verifica se ha disturbio programado para o ciclo 0 (o primeiro ciclo eh o 1)
  {
    stop("\n\tAtenção! O primeiro ciclo das simulações é o ciclo 1. Atribua 0 ao dist.pos apenas quando não desejar implementar a ocorrência de distúrbios.\n\n")
  }
  #############################################################################
  ###################### Matrizes para guardar resultados #####################
  #############################################################################
  ind.mat=matrix(nrow=J,ncol=1+ciclo/step) ### Gera matriz da identidade (especie) de cada individuo por ciclo    
  prop.mat=matrix(nrow=J,ncol=1+ciclo/step) ### Gera matriz de propagulos produzidos por individuo em cada ciclo
  dead.mat=matrix(nrow=J,ncol=1+ciclo/step) ### Gera matriz de probabilidade de morte de cada individuo, por ciclo
  ### Vetor com numero de mortos por ciclo
  n.dead <- c()
  n.dead[1] <- 0
  #############################################################################
  ############################# Condicoes iniciais ############################
  #############################################################################
  ### Guarda a identidade inicial de cada individuo
  ind.mat[,1] <- rep(1:S,each=j)
  cod.sp <- ind.mat[,1] ### Transfere informacao para vetor temporario que se atualiza a cada ciclo, depois de ser copiado para uma coluna da matriz
  ### Guarda a probabilidade de morte de cada individuo. A probabilidade inicial de morte eh tomada de uma geometrica considerando que o numero de mortes esperado por ciclo eh a mesma para todos (p=1/J)
  dead.mat[,1] <- 1/J
  p.death <- dead.mat[,1] ### Transfere informacao para vetor temporario que se atualiza a cada ciclo, depois de ser copiado para uma coluna da matriz
  ### Guarda o numero de propagulos produzidos por ciclo de cada individuo
  prop.mat[,1] <- xi0
  n.propag <- prop.mat[,1] ### Transfere informacao para vetor temporario que se atualiza a cada ciclo, depois de ser copiado para uma coluna da matriz
  #############################################################################
  ############ Contador para salvar resultados a cada step ciclos #############
  #############################################################################
  sc=2 ### Contador que comeca na posicao 2, ja que na primeira posicao foram guardadas as condicoes iniciais
  #############################################################################
  ############################## Inicio do ciclo ##############################
  #############################################################################
  for(i in 1:ciclo)
    {
    #n.mortes <- 0
    morte=rbinom(J, 1, prob=p.death) ### Sorteio dos individuos que morrerao
    #########################################################################
    ######################### Inicio dos disturbios #########################
    #########################################################################
    if(sum(dist.pos==i)>0) ### Identifica se ocorre um evento de disturbio no ciclo atual
      {
      vivos <- which(morte==0) ### Identifica individuos sobreviventes
      nvivos <- length(vivos) ### Conta numero de individuos sobreviventes
      if(length(dist.int)>1) ### Identifica se os eventos apresentam diferentes intensidades
        {
        posdist <- which(dist.pos==i) ### Guarda qual o numero do evento do disturbio
        ndist <- round(nvivos* dist.int[posdist]) ### Calcula o numero de individuos mortos com o evento
        }
      if(length(dist.int)==1) ### Identifica se os eventos apresentam a mesma intensidade
        {
        ndist <- round(nvivos* dist.int) ### Calcula o numero de individuos mortos com o evento
        }
      posmort <- sample(vivos, ndist) ### Sorteia quais individuos serao mortos no evento de disturbio
      morte[posmort] <- 1 ### Marca os individuos que morreram como mortos (numero 1 no vetor 0/1)
      }
    #########################################################################
    ######################### Termino dos disturbios ########################
    #########################################################################
    #D=sum(morte)
    #n.mortes <- n.mortes+D
    n.mortes <- sum(morte) ### Grava o numero total de mortes no ciclo
    #########################################################################
    ######################## Substituicao de valores ########################
    #########################################################################
    if(n.mortes>0) ### Identifica se houve mortes no ciclo
      {
      seed.bank <- rep(1:J,round(n.propag)) ### Banco de propagulos: cada propagulo tem o codigo numerico do individuo. Como o fenotipo n.propag pode ter valores nao inteiros, arredondamos.
      nascer= which(morte==1) ### Armazena indices dos individuos que morreram
      mami=sample(seed.bank, n.mortes) ### Sorteia os propagulos que irao repor os mortos
      papi <- c() ### Cria vetor para armazenar o fenotipo do pai
      for(w in 1:n.mortes) ### Cria loop para sortear o pai entre os individuos da especie de cada propagulo-mae sorteado
        {
        papi[w] <- sample(n.propag[ cod.sp==cod.sp[mami[w]] ],1)
        }
      medias.prop=(n.propag[mami]+papi)/2 ### Calcula o valor esperado de propagulos produzidos por ciclo dos filhotes, representado pela media do numero medio de propagulos produzidos pelos parentais
      cod.sp[nascer]<-cod.sp[mami] ### Substitui codigos das especies dos mortos pelos codigos dos individuos novos
      n.propag[nascer] <- sapply(1,rtruncnorm,a=1, b=X , mean= medias.prop,sd=dp) ### Sorteia o numero de propagulos produzidos por ciclo pelos novos individuos de uma distribuicao normal discretizada e truncada entre 1 e X
      p.death[nascer] <- n.propag[nascer]/X ### Atualiza a matriz de probabilidades de morrer
      }
    #########################################################################
    ####################### Salvamento dos resultados #######################
    #########################################################################
    if(sum(i==seq(step,ciclo,step))==1) ### Verifica se o ciclo atual eh um dos que devem ser salvos
      {
      ind.mat[,sc] <- cod.sp ### Guarda na posicao sc a identificacao nova dos individuos
      dead.mat[,sc] <- p.death ### Guarda na posicao sc a probabilidade de morte nova dos individuos
      prop.mat[,sc] <- round(n.propag) ### Guarda na posicao sc o numero de propagulos produzidos por ciclo novo dos individuos
      n.dead[sc] <- n.mortes ### Guarda na posicao sc o numero de mortes do ciclo atual
      sc <- sc+1 ### Atualiza o contador que salva os resultados para o proximo ciclo a ser rodado
    ##cat(format(Sys.time(), "%d%b%Y_%H:%M"), "\t ciclo = ", i, "\n") # para avisar a cada ciclo! desligar se estiver usando Rcloud
      } 
  }
  #############################################################################
  ############################## Termino do ciclo #############################
  #############################################################################
  ########################### Organizacao do output ###########################
  #############################################################################
  tempo <- seq(0,ciclo,by=step) ### Cria vetor que dara nome as colunas das matrizes 
  colnames(ind.mat) <- tempo ### Nomeia colunas da matriz ind.mat
  colnames(dead.mat) <- tempo ### Nomeia colunas da matriz dead.mat
  colnames(prop.mat) <- tempo ### Nomeia colunas da matriz prop.mat
  names(n.dead) <- tempo ### Nomeia os elementos do vetor
  resulta=list(tempo=tempo,sp.list=ind.mat,sementes=prop.mat,prob.morte=dead.mat,n.mortes=n.dead)
  t1=proc.time()[[3]] ### Marca o termino da contagem de tempo de processamento da funcao
  cat("\n\t tempo de processamento: ", round((t1-t0)/60,2),"\n") ### Mostra o tempo de processamento no console
  attributes(resulta)$start=list(especies=S, individuos=j, nprop=X, sd=dp, posicao_disturbios=dist.pos, intensidade_disturbios=dist.int, ciclos=ciclo, passos=step) ### Inclui atributos no objeto resulta
  return(resulta) ### Retorna o objeto resulta
}
#############################################################################
############################# Termino da funcao #############################
#############################################################################

# Normal Discretizada Truncada
#Sorteia n numeros inteiros de uma aproximacao discreta da normal truncada em seus limites inferiores 
#e superiores dados sua media e coeficiente de variacao. 
##Funcao de sorteio de uma normal discretizada e truncada
# rnormt <- function(mean,n=1,dp,min,max)
# {
#   vals <- (min-0.5):(max+0.5)
#   p <- pnorm(vals,mean=mean,sd=dp)
#   p2 <- diff(p)
#   sample(min:max,n,prob=p2,replace=T)


#### em vez de usar a rnormt, usar o pacote truncnorm
