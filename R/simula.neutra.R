############################## Pacote necessario #############################
require(truncnorm)
##############################################################################

##############################################################################
######### Evolucao de estrategias de vida a partir de cenario neutro #########
################################# "GERALZONA" ################################
####### VERSAO COM TODAS AS POSSIBILIDADES TRANSFORMADAS EM ARGUMENTOS #######
##############################################################################
################################ simula.neutra ###############################
##############################################################################
############################ Argumentos da funcao ############################
##############################################################################
# S : numero de especies da comunidade (escalar)
# j: abundancia inicial das especies (escalar)
# xi0: numero de propagulos produzidos a cada ciclo pelos individuos (vetor de tamanho J)
# X: numero maximo de propagulos que qualquer individuo pode produzir a cada ciclo (escalar)
# pi0: probabilidade de morte/ciclo inicial dos individuos (vetor de tamanho J) - deve conter algum valor apenas quando trade.off=F
# dp: desvio padrao da(s) distribuicao(oes) normal(is) gerada(s) para o sorteio do numero de propagulos por ciclo (e da probabilidade de morte a cada ciclo) de individuos novos a partir da(s) media(s) do(s) caracter(es) dos parentais (escalar)
# dist.pos: identificacao dos ciclos em que ocorrem eventos de disturbio (vetor cujo tamanho corresponde ao numero de eventos de disturbio que ocorre ao longo da simulacao)
# dist.int: intensidade de cada evento de disturbio. Se ha soh um valor este eh aplicado para todos (vetor cujo tamanho corresponde ao numero de eventos de disturbio que ocorre ao longo da simulacao)
# ciclo: numero de ciclos da simulacao (escalar)
# step: intervalo de ciclos a cada qual sao salvos resultados da simulacao (escalar)
# trade.off: determina se o calculo da probabilidade de morte dos individuos novos ocorre a partir do valor de numero de propagulos produzidos por ciclo (xi0/X) ou se eh independente deste caracter (vetor logico T/F)
# prop.neutro: determina se o numero de propagulos produzidos por ciclo eh um dos componentes do fitness ou se eh um caracter neutro (vetor logico T/F) - deve ser TRUE apenas quando trade.off=F
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
simula.neutra=function(S= 100, j=10, xi0=rep(seq(10,10,length.out = S),each=j), X=10000, pi0=NULL, dp=0.1, dist.pos=NULL, dist.int=NULL, ciclo=1e6, step=100, trade.off=T, prop.neutro=F)
{
  t0=proc.time()[[3]] ### Marca o inicio da contagem de tempo de processamento da funcao
  cat("Inicio simulacao... Ciclos: ")
  #############################################################################
  ########################### Argumentos deduzidos ############################
  #############################################################################
  J <- S*j ### Calcula o tamanho da comunidade (J)
  #X <- min(xi0)*J ### Calcula o numero total de propagulos produzidos por um individuo (X).
  #############################################################################
  ############################### Verificacoes ################################
  #############################################################################
  if(max(xi0)>X) ### Verifica se o valor/algum dos valores de xi0 (numero de propagulos produzidos por um indivíduo por ciclo) eh maior do que o numero de propagulos produzidos ao longo da vida toda do individuo
  {
    stop("\n\tO número de propágulos produzidos por indivíduo a cada ciclo (xi0) não pode ser maior do que o número de propágulos máximo (X). Tente novamente!\n\n")
  }
  if(sum(dist.int>1)>0) ### Verifica se ha disturbio programado com valor de intensidade que nao faz sentido
  {
    stop("\n\tA intensidade do distúrbio deve ser entre 0 e 1. Tente novamente!\n\n")
  } 
  if(sum(0<=dist.pos & dist.pos<=ciclo)<length(dist.pos)) ### Verifica se ha disturbio programado para ciclo inexistente
  {
    stop("\n\tA posição dos eventos de distúrbio (dist.pos) precisa ser condizente com o número de ciclos a serem rodados (ciclo). Tente novamente!\n\n")
  } 
  if(sum(dist.pos==0 & dist.int>0)>0) ### Verifica se ha disturbio programado para o ciclo 0 (o primeiro ciclo eh o 1)
  {
    stop("\n\tAtenção! O primeiro ciclo das simulações é o ciclo 1. Atribua 0 ao dist.pos apenas quando não desejar implementar a ocorrência de distúrbios.\n\n")
  }
  if(step>ciclo) ### Verifica se o intervalo de ciclos entre os quais os resultados devem ser guardados esta contido no numero de ciclos total da simulacao
  {
    stop("\n\tAtenção! O valor do argumento 'step' não pode ser maior que o de 'ciclo'.\n\n")
  }
  if(trade.off==T & prop.neutro==T) ### Verifica se estah programado trade-off e neutralidade do caracter de producao de propagulos por ciclo
  {
    stop("\n\tAtenção! Para que o caracter referente à produção de propágulos seja neutro, não pode haver trade-off entre este caracter e a probabilidade de morte, i.e., o argumento trade.off deve ser 'FALSE'.\n\n")
  }
  #############################################################################
  ###################### Matrizes para guardar resultados #####################
  #############################################################################
  ind.mat=matrix(nrow=J,ncol=1+ciclo/step) ### Gera matriz da identidade (especie) de cada individuo por ciclo    
  prop.mat=matrix(nrow=J,ncol=1+ciclo/step) ### Gera matriz de propagulos produzidos por individuo em cada ciclo
  dead.mat=matrix(nrow=J,ncol=1+ciclo/step) ### Gera matriz de probabilidade de morte de cada individuo, por ciclo
  ### Vetor com numero de mortes cumulativo ate determinado ciclo
  n.dead <- 0
  n.dead.vetor<-c()
  n.dead.vetor[1] <- 0
  #############################################################################
  ############################# Condicoes iniciais ############################
  #############################################################################
  ### Guarda a identidade inicial de cada individuo
  ind.mat[,1] <- rep(1:S,each=j)
  cod.sp <- ind.mat[,1] ### Transfere informacao para vetor temporario que se atualiza a cada ciclo, depois de ser copiado para uma coluna da matriz
  ### Guarda a probabilidade de morte de cada individuo. A probabilidade inicial de morte eh tomada de uma geometrica considerando que o numero de mortes esperado por ciclo eh a mesma para todos (p=1/J)
  if(trade.off==T)
  {
    dead.mat[,1] <- xi0/X
  }
  if(trade.off==F)
  {
    dead.mat[,1] <- pi0
  }
  p.death <- dead.mat[,1] ### Transfere informacao para vetor temporario que se atualiza a cada ciclo, depois de ser copiado para uma coluna da matriz
  ### Guarda o numero de propagulos produzidos por ciclo de cada individuo
  prop.mat[,1] <- xi0 ### Transfere informacao para vetor temporario que se atualiza a cada ciclo, depois de ser copiado para uma coluna da matriz
  n.propag <- prop.mat[,1]
  prop.mat[,1]<-round(prop.mat[,1]) ### Arredonda valor de numero de propagulos produzidos por ciclo na matriz de resultados para que se retorne como saida da funcao o "fenotipo" dos individuos
  #############################################################################
  ############ Contador para salvar resultados a cada step ciclos #############
  #############################################################################
  sc=2 ### Contador que comeca na posicao 2, ja que na primeira posicao foram guardadas as condicoes iniciais
  #############################################################################
  ############################## Inicio do ciclo ##############################
  #############################################################################
  for(i in 1:ciclo)
  {
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
    n.mortes <- sum(morte) ### Grava o numero total de mortes no ciclo
    n.dead <- n.dead + n.mortes ### Atualiza n.dead com numero de mortes acumulado ate entao mais o numero de mortes do ciclo atual
    #########################################################################
    ######################## Substituicao de valores ########################
    #########################################################################
    if(n.mortes>0) ### Identifica se houve mortes no ciclo
    {
      cod.ind<-1:J
      seed.bank <- rep(cod.ind,round(n.propag)) ### Banco de propagulos: cada propagulo tem o codigo numerico do individuo. Como o fenotipo n.propag pode ter valores nao inteiros, arredondamos.
      nascer= which(morte==1) ### Armazena indices dos individuos que morreram
      if(prop.neutro==F)
      {
        mami=sample(seed.bank, n.mortes) ### Sorteia os propagulos que irao repor os mortos
      }
      if(prop.neutro==T)
      {
        mami=sample(cod.ind, n.mortes) ### Sorteia os individuos adultos que irao repor os mortos
      }
      papi <- c() ### Cria vetor para armazenar o fenotipo do pai
      for(w in 1:n.mortes) ### Cria loop para sortear o pai entre os individuos da especie de cada propagulo-mae sorteado
      {
        papi[w] <- sample(cod.ind[ cod.sp==cod.sp[mami[w]] ],1)
      }
      medias.prop=(n.propag[mami]+n.propag[papi])/2 ### Calcula o valor esperado de propagulos produzidos por ciclo dos filhotes, representado pela media do numero medio de propagulos produzidos pelos parentais
      cod.sp[nascer]<-cod.sp[mami] ### Substitui codigos das especies dos mortos pelos codigos dos individuos novos
      n.propag[nascer] <- sapply(1,rtruncnorm,a=(1-0.000001), b=(X+0.000001), mean= medias.prop,sd=dp) ### Sorteia o numero de propagulos produzidos por ciclo pelos novos individuos de uma distribuicao normal discretizada e truncada entre 1 e X e cuja media eh a media dos pais
      if (trade.off==T)
      {
        p.death[nascer] <- n.propag[nascer]/X ### Calcula a probabilidade de morte dos individuos novos a partir do trade-off
      }
      if (trade.off==F)
      {
        medias.p.death <- (p.death[mami]+p.death[papi])/2 ### Calcula a media da probabilidade de morte dos pais
        p.death[nascer] <- sapply(1,rtruncnorm,a=1/200000, b=1, mean=medias.p.death,sd=dp) ### Sorteia a probabilidade de morte dos novos individuos de uma distribuicao normal discretizada e truncada entre um numero muito pequeno e 1 e cuja media eh a media dos pais
      }
    }
    #########################################################################
    ####################### Salvamento dos resultados #######################
    #########################################################################
    if(sum(i==seq(step,ciclo,step))==1) ### Verifica se o ciclo atual eh um dos que devem ser salvos
    {
      ind.mat[,sc] <- cod.sp ### Guarda na posicao sc a identificacao nova dos individuos
      dead.mat[,sc] <- p.death ### Guarda na posicao sc a probabilidade de morte nova dos individuos
      prop.mat[,sc] <- round(n.propag) ### Guarda na posicao sc o numero de propagulos produzidos por ciclo novo dos individuos
      n.dead.vetor[sc] <- n.dead ### Guarda na posicao sc o numero de mortes acumulado ate o ciclo atual
      sc <- sc+1 ### Atualiza o contador que salva os resultados para o proximo ciclo a ser rodado
      ##cat(format(Sys.time(), "%d%b%Y_%H:%M"), "\t ciclo = ", i, "\n") # para avisar a cada ciclo! desligar se estiver usando Rcloud
      cat(i," ")
    } 
  }
  #############################################################################
  ############################## Termino do ciclo #############################
  #############################################################################
  cat("...Termino simulacao\n")
  #############################################################################
  ########################### Organizacao do output ###########################
  #############################################################################
  tempo <- seq(0,ciclo,by=step) ### Cria vetor que dara nome as colunas das matrizes 
  colnames(ind.mat) <- tempo ### Nomeia colunas da matriz ind.mat
  colnames(dead.mat) <- tempo ### Nomeia colunas da matriz dead.mat
  colnames(prop.mat) <- tempo ### Nomeia colunas da matriz prop.mat
  names(n.dead.vetor) <- tempo ### Nomeia os elementos do vetor
  resulta=list(tempo=tempo,sp.list=ind.mat,sementes=prop.mat,prob.morte=dead.mat,n.mortes.cumulativo=n.dead.vetor)
  t1=proc.time()[[3]] ### Marca o termino da contagem de tempo de processamento da funcao
  cat("Tempo de processamento: ", round((t1-t0)/60,2),"min\n") ### Mostra o tempo de processamento no console
  attributes(resulta)$start=list(especies=S, individuos=j, nprop=xi0, sd=dp, posicao_disturbios=dist.pos, intensidade_disturbios=dist.int, ciclos=ciclo, passos=step) ### Inclui atributos no objeto resulta
  return(resulta) ### Retorna o objeto resulta
}
#############################################################################
############################# Termino da funcao #############################
#############################################################################


##############################################################################
######### Evolucao de estrategias de vida a partir de cenario neutro #########
##############################################################################
############################# simula.neutra.trade ############################
##############################################################################
############################ Listagem das versoes ############################
##############################################################################
# 7. Versao modificada por Luisa Novara para incluir X como um argumento e alguns erros de parada (maio 2016)
# 6. Versao modificada por Luisa Novara para incluir spp com diferentes estrategias de vida iniciais (maio 2016)
# 5. Versao modificada por Luisa Novara para incluir retorno do numero de mortes cumulativo (maio 2016)
# 4. Versao modificada por Luisa Novara e Alexandre Adalardo (fevereiro 2016)
# 3. Disturbio implementado por Alexandre Adalardo e Luisa Novara (novembro 2015) 
# 2. Versao modificada por Luisa Novara (2014)
# 1. Versao inicial de Alexandre Adalardo (outubro 2009) modificada por Paulo Inacio Prado (outubro 2009)
##############################################################################
########################## Detalhamento das versoes ##########################
##############################################################################
# 7. Desliga a linha de comando que calcula X a partir de xi0J e inclui X como um argumento da funcao, para evitar que a estrategia de vida inicial (xi0) esteja sempre muito distante do maximo possivel (X)
# 6. Troca a forma do argumento xi0 de um valor para um vetor de tamanho J, com os valores das estrategias iniciais de cada individuo da comunidade.
# 4. Troca o argumento X por xi0, para que X (que eh dado por xi0*J) seja inteiro
# 2. Troca o antigo argumento cv (coeficiente de variacao) pelo dp (desvio padrao) da distribuicao normal da herdabilidade de xi. Isso permite que a herdabilidade permaneça constante e evita o erro de gerar herdabilidade maior quando os valores de xi sao mais baixos (i.e., evita que os valores de xi da populacao fiquem "presos" em valores mais baixos e isso mascare os reais resultados das simulacoes)
##############################################################################
############################## Inicio da funcao ##############################
##############################################################################
simula.neutra.trade=function(S= 100, j=10, xi0=rep(seq(10,10,length.out = S),each=j), X=10000, dp=0.1, dist.pos=NULL, dist.int=NULL, ciclo=1e6, step=100)
{
  t0=proc.time()[[3]] ### Marca o inicio da contagem de tempo de processamento da funcao
  cat("Inicio simulacao... Ciclos: ")
  #############################################################################
  ########################### Argumentos deduzidos ############################
  #############################################################################
  J <- S*j ### Calcula o tamanho da comunidade (J)
  #X <- min(xi0)*J ### Calcula o numero total de propagulos produzidos por um individuo (X).
  #############################################################################
  ############################### Verificacoes ################################
  #############################################################################
  if(max(xi0)>X) ### Verifica se o valor/algum dos valores de xi0 (numero de propagulos produzidos por um indivíduo por ciclo) eh maior do que o numero de propagulos produzidos ao longo da vida toda do individuo
  {
    stop("\n\tO número de propágulos produzidos por indivíduo a cada ciclo (xi0) não pode ser maior do que o número de propágulos máximo (X). Tente novamente!\n\n")
  }
  if(sum(dist.int>1)>0) ### Verifica se ha disturbio programado com valor de intensidade que nao faz sentido
  {
    stop("\n\tA intensidade do distúrbio deve ser entre 0 e 1. Tente novamente!\n\n")
  } 
  if(sum(0<=dist.pos & dist.pos<=ciclo)<length(dist.pos)) ### Verifica se ha disturbio programado para ciclo inexistente
  {
    stop("\n\tA posição dos eventos de distúrbio (dist.pos) precisa ser condizente com o número de ciclos a serem rodados (ciclo). Tente novamente!\n\n")
  } 
  if(sum(dist.pos==0 & dist.int>0)>0) ### Verifica se ha disturbio programado para o ciclo 0 (o primeiro ciclo eh o 1)
  {
    stop("\n\tAtenção! O primeiro ciclo das simulações é o ciclo 1. Atribua 0 ao dist.pos apenas quando não desejar implementar a ocorrência de distúrbios.\n\n")
  }
  if(step>ciclo) ### Verifica se o intervalo de ciclos entre os quais os resultados devem ser guardados esta contido no numero de ciclos total da simulacao
  {
    stop("\n\tAtenção! O valor do argumento 'step' não pode ser maior que o de 'ciclo'.\n\n")
  }
  #############################################################################
  ###################### Matrizes para guardar resultados #####################
  #############################################################################
  ind.mat=matrix(nrow=J,ncol=1+ciclo/step) ### Gera matriz da identidade (especie) de cada individuo por ciclo    
  prop.mat=matrix(nrow=J,ncol=1+ciclo/step) ### Gera matriz de propagulos produzidos por individuo em cada ciclo
  dead.mat=matrix(nrow=J,ncol=1+ciclo/step) ### Gera matriz de probabilidade de morte de cada individuo, por ciclo
  ### Vetor com numero de mortes cumulativo ate determinado ciclo
  n.dead <- 0
  n.dead.vetor<-c()
  n.dead.vetor[1] <- 0
  #############################################################################
  ############################# Condicoes iniciais ############################
  #############################################################################
  ### Guarda a identidade inicial de cada individuo
  ind.mat[,1] <- rep(1:S,each=j)
  cod.sp <- ind.mat[,1] ### Transfere informacao para vetor temporario que se atualiza a cada ciclo, depois de ser copiado para uma coluna da matriz
  ### Guarda a probabilidade de morte de cada individuo (calculada a partir do trade-off)
  dead.mat[,1] <- xi0/X
  p.death <- dead.mat[,1] ### Transfere informacao para vetor temporario que se atualiza a cada ciclo, depois de ser copiado para uma coluna da matriz
  ### Guarda o numero de propagulos produzidos por ciclo de cada individuo
  prop.mat[,1] <- xi0 ### Transfere informacao para vetor temporario que se atualiza a cada ciclo, depois de ser copiado para uma coluna da matriz
  n.propag <- prop.mat[,1]
  prop.mat[,1]<-round(prop.mat[,1]) ### Arredonda valor de numero de propagulos produzidos por ciclo na matriz de resultados para que se retorne como saida da funcao o "fenotipo" dos individuos
  #############################################################################
  ############ Contador para salvar resultados a cada step ciclos #############
  #############################################################################
  sc=2 ### Contador que comeca na posicao 2, ja que na primeira posicao foram guardadas as condicoes iniciais
  #############################################################################
  ############################## Inicio do ciclo ##############################
  #############################################################################
  for(i in 1:ciclo)
  {
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
    n.mortes <- sum(morte) ### Grava o numero total de mortes no ciclo
    n.dead <- n.dead + n.mortes ### Atualiza n.dead com numero de mortes acumulado ate entao mais o numero de mortes do ciclo atual
    #########################################################################
    ######################## Substituicao de valores ########################
    #########################################################################
    if(n.mortes>0) ### Identifica se houve mortes no ciclo
    {
      cod.ind<-1:J
      seed.bank <- rep(cod.ind,round(n.propag)) ### Banco de propagulos: cada propagulo tem o codigo numerico do individuo. Como o fenotipo n.propag pode ter valores nao inteiros, arredondamos.
      nascer= which(morte==1) ### Armazena indices dos individuos que morreram
      mami=sample(seed.bank, n.mortes) ### Sorteia os propagulos que irao repor os mortos
      papi <- c() ### Cria vetor para armazenar o fenotipo do pai
      for(w in 1:n.mortes) ### Cria loop para sortear o pai entre os individuos da especie de cada propagulo-mae sorteado
      {
        papi[w] <- sample(cod.ind[ cod.sp==cod.sp[mami[w]] ],1)
      }
      medias.prop=(n.propag[mami]+n.propag[papi])/2 ### Calcula o valor esperado de propagulos produzidos por ciclo dos filhotes, representado pela media do numero medio de propagulos produzidos pelos parentais
      cod.sp[nascer]<-cod.sp[mami] ### Substitui codigos das especies dos mortos pelos codigos dos individuos novos
      n.propag[nascer] <- sapply(1,rtruncnorm,a=(1-0.000001), b=(X+0.000001), mean= medias.prop,sd=dp) ### Sorteia o numero de propagulos produzidos por ciclo pelos novos individuos de uma distribuicao normal discretizada e truncada entre 1 e X
      p.death[nascer] <- n.propag[nascer]/(X+0.000001) ### Calcula a probabilidade de morte dos individuos novos a partir do trade-off
    }
    #########################################################################
    ####################### Salvamento dos resultados #######################
    #########################################################################
    if(sum(i==seq(step,ciclo,step))==1) ### Verifica se o ciclo atual eh um dos que devem ser salvos
    {
      ind.mat[,sc] <- cod.sp ### Guarda na posicao sc a identificacao nova dos individuos
      dead.mat[,sc] <- p.death ### Guarda na posicao sc a probabilidade de morte nova dos individuos
      prop.mat[,sc] <- round(n.propag) ### Guarda na posicao sc o numero de propagulos produzidos por ciclo novo dos individuos
      n.dead.vetor[sc] <- n.dead ### Guarda na posicao sc o numero de mortes acumulado ate o ciclo atual
      sc <- sc+1 ### Atualiza o contador que salva os resultados para o proximo ciclo a ser rodado
      ##cat(format(Sys.time(), "%d%b%Y_%H:%M"), "\t ciclo = ", i, "\n") # para avisar a cada ciclo! desligar se estiver usando Rcloud
      cat(i," ")
    } 
  }
  #############################################################################
  ############################## Termino do ciclo #############################
  #############################################################################
  cat("...Termino simulacao\n")
  #############################################################################
  ########################### Organizacao do output ###########################
  #############################################################################
  tempo <- seq(0,ciclo,by=step) ### Cria vetor que dara nome as colunas das matrizes 
  colnames(ind.mat) <- tempo ### Nomeia colunas da matriz ind.mat
  colnames(dead.mat) <- tempo ### Nomeia colunas da matriz dead.mat
  colnames(prop.mat) <- tempo ### Nomeia colunas da matriz prop.mat
  names(n.dead.vetor) <- tempo ### Nomeia os elementos do vetor
  resulta=list(
    tempo=tempo,
    sp.list=ind.mat,
    sementes=prop.mat,
    prob.morte=dead.mat,
    n.mortes.cumulativo=n.dead.vetor)
  t1=proc.time()[[3]] ### Marca o termino da contagem de tempo de processamento da funcao
  cat("Tempo de processamento: ", round((t1-t0)/60,2),"min\n") ### Mostra o tempo de processamento no console
  attributes(resulta)$start=list(especies=S, individuos=j, nprop=xi0, nprop_max=X, sd=dp, posicao_disturbios=dist.pos, intensidade_disturbios=dist.int, ciclos=ciclo, passos=step) ### Inclui atributos no objeto resulta
  return(resulta) ### Retorna o objeto resulta
}
#############################################################################
############################# Termino da funcao #############################
#############################################################################


##############################################################################
######### Evolucao de estrategias de vida a partir de cenario neutro #########
###################### SEM TRADE-OFF - DEMONIO DE DARWIN #####################
##############################################################################
############################## simula.neutra.dd ##############################
##############################################################################
############################ Listagem das versoes ############################
##############################################################################
# 9. Versao modificada por Luisa Novara para incluir pi0 como um argumento (maio 2016)
# 8. Versao modificada por Luisa Novara para incluir X como um argumento e alguns erros de parada (maio 2016)
# 7. Versao modificada por Luisa Novara para retirar o trade-off entre a fecundidade e a probabilidade de morte (maio 2016)
# 6. Versao modificada por Luisa Novara para incluir spp com diferentes estrategias de vida iniciais (maio 2016)
# 5. Versao modificada por Luisa Novara para incluir retorno do numero de mortes cumulativo (maio 2016)
# 4. Versao modificada por Luisa Novara e Alexandre Adalardo (fevereiro 2016)
# 3. Disturbio implementado por Alexandre Adalardo e Luisa Novara (novembro 2015) 
# 2. Versao modificada por Luisa Novara (2014)
# 1. Versao inicial de Alexandre Adalardo (outubro 2009) modificada por Paulo Inacio Prado (outubro 2009)
##############################################################################
########################## Detalhamento das versoes ##########################
##############################################################################
# 9. Retira o calculo da probabilidade de morte inicial de 1/J e substitui pela entrada de valores a escolha do usuario (argumento pi0)
# 8. Desliga a linha de comando que calcula X a partir de xi0J e inclui X como um argumento da funcao, para evitar que a estrategia de vida inicial (xi0) esteja sempre muito distante do maximo possivel (X)
# 7. Troca o calculo da probabilidade de morte de individuos novos de xi0/X (trade-off) para um sorteio de uma distribuicao normal cuja media eh a media dos parentais (assim como eh feito para o numero de propagulos produzidos por ciclo)
# 6. Troca a forma do argumento xi0 de um valor para um vetor de tamanho J, com os valores das estrategias iniciais de cada individuo da comunidade
# 4. Troca o argumento X por xi0, para que X (que eh dado por xi0*J) seja inteiro
# 2. Troca o antigo argumento cv (coeficiente de variacao) pelo dp (desvio padrao) da distribuicao normal da herdabilidade de xi. Isso permite que a herdabilidade permaneça constante e evita o erro de gerar herdabilidade maior quando os valores de xi sao mais baixos (i.e., evita que os valores de xi da populacao fiquem "presos" em valores mais baixos e isso mascare os reais resultados das simulacoes)
##############################################################################
############################## Inicio da funcao ##############################
##############################################################################
simula.neutra.dd=function(S= 100, j=10, xi0=rep(seq(10,10,length.out = S),each=j), X=10000, pi0=rep(seq(0.01,0.01,length.out = S),each=j), dp=0.1, dist.pos=NULL, dist.int=NULL, ciclo=1e6, step=100)
{
  t0=proc.time()[[3]] ### Marca o inicio da contagem de tempo de processamento da funcao
  #cat("Inicio simulacao... Ciclos: ")
  #############################################################################
  ########################### Argumentos deduzidos ############################
  #############################################################################
  J <- S*j ### Calcula o tamanho da comunidade (J)
  #X <- min(xi0)*J ### Calcula o numero total de propagulos produzidos por um individuo (X).
  #############################################################################
  ############################### Verificacoes ################################
  #############################################################################
  if(max(xi0)>X) ### Verifica se o valor/algum dos valores de xi0 (numero de propagulos produzidos por um indivíduo por ciclo) eh maior do que o numero de propagulos produzidos ao longo da vida toda do individuo
  {
    stop("\n\tO número de propágulos produzidos por indivíduo a cada ciclo (xi0) não pode ser maior do que o número de propágulos máximo (X). Tente novamente!\n\n")
  }
  if(sum(dist.int>1)>0) ### Verifica se ha disturbio programado com valor de intensidade que nao faz sentido
  {
    stop("\n\tA intensidade do distúrbio deve ser entre 0 e 1. Tente novamente!\n\n")
  } 
  if(sum(0<=dist.pos & dist.pos<=ciclo)<length(dist.pos)) ### Verifica se ha disturbio programado para ciclo inexistente
  {
    stop("\n\tA posição dos eventos de distúrbio (dist.pos) precisa ser condizente com o número de ciclos a serem rodados (ciclo). Tente novamente!\n\n")
  } 
  if(sum(dist.pos==0 & dist.int>0)>0) ### Verifica se ha disturbio programado para o ciclo 0 (o primeiro ciclo eh o 1)
  {
    stop("\n\tAtenção! O primeiro ciclo das simulações é o ciclo 1. Atribua 0 ao dist.pos apenas quando não desejar implementar a ocorrência de distúrbios.\n\n")
  }
  if(step>ciclo) ### Verifica se o intervalo de ciclos entre os quais os resultados devem ser guardados esta contido no numero de ciclos total da simulacao
  {
    stop("\n\tAtenção! O valor do argumento 'step' não pode ser maior que o de 'ciclo'.\n\n")
  }
  #############################################################################
  ###################### Matrizes para guardar resultados #####################
  #############################################################################
  ind.mat=matrix(nrow=J,ncol=1+ciclo/step) ### Gera matriz da identidade (especie) de cada individuo por ciclo    
  prop.mat=matrix(nrow=J,ncol=1+ciclo/step) ### Gera matriz de propagulos produzidos por individuo em cada ciclo
  dead.mat=matrix(nrow=J,ncol=1+ciclo/step) ### Gera matriz de probabilidade de morte de cada individuo, por ciclo
  ### Vetor com numero de mortes cumulativo ate determinado ciclo
  n.dead <- 0
  n.dead.vetor<-c()
  n.dead.vetor[1] <- 0
  #############################################################################
  ############################# Condicoes iniciais ############################
  #############################################################################
  ### Guarda a identidade inicial de cada individuo
  ind.mat[,1] <- rep(1:S,each=j)
  cod.sp <- ind.mat[,1] ### Transfere informacao para vetor temporario que se atualiza a cada ciclo, depois de ser copiado para uma coluna da matriz
  ### Guarda a probabilidade de morte de cada individuo. A probabilidade inicial de morte eh tomada de uma geometrica considerando que o numero de mortes esperado por ciclo eh a mesma para todos (p=1/J)
  dead.mat[,1] <- pi0
  p.death <- dead.mat[,1] ### Transfere informacao para vetor temporario que se atualiza a cada ciclo, depois de ser copiado para uma coluna da matriz
  ### Guarda o numero de propagulos produzidos por ciclo de cada individuo
  prop.mat[,1] <- xi0 ### Transfere informacao para vetor temporario que se atualiza a cada ciclo, depois de ser copiado para uma coluna da matriz
  n.propag <- prop.mat[,1]
  prop.mat[,1]<-round(prop.mat[,1]) ### Arredonda valor de numero de propagulos produzidos por ciclo na matriz de resultados para que se retorne como saida da funcao o "fenotipo" dos individuos
  #############################################################################
  ############ Contador para salvar resultados a cada step ciclos #############
  #############################################################################
  sc=2 ### Contador que comeca na posicao 2, ja que na primeira posicao foram guardadas as condicoes iniciais
  #############################################################################
  ############################## Inicio do ciclo ##############################
  #############################################################################
  for(i in 1:ciclo)
  {
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
    n.mortes <- sum(morte) ### Grava o numero total de mortes no ciclo
    n.dead <- n.dead + n.mortes ### Atualiza n.dead com numero de mortes acumulado ate entao mais o numero de mortes do ciclo atual
    #########################################################################
    ######################## Substituicao de valores ########################
    #########################################################################
    if(n.mortes>0) ### Identifica se houve mortes no ciclo
    {
      cod.ind<-1:J
      seed.bank <- rep(cod.ind,round(n.propag)) ### Banco de propagulos: cada propagulo tem o codigo numerico do individuo. Como o fenotipo n.propag pode ter valores nao inteiros, arredondamos.
      nascer= which(morte==1) ### Armazena indices dos individuos que morreram
      mami=sample(seed.bank, n.mortes) ### Sorteia os propagulos que irao repor os mortos
      papi <- c() ### Cria vetor para armazenar o fenotipo do pai
      for(w in 1:n.mortes) ### Cria loop para sortear o pai entre os individuos da especie de cada propagulo-mae sorteado
      {
        papi[w] <- sample(cod.ind[ cod.sp==cod.sp[mami[w]] ],1)
      }
      medias.prop=(n.propag[mami]+n.propag[papi])/2 ### Calcula o valor esperado de propagulos produzidos por ciclo dos filhotes, representado pela media do numero medio de propagulos produzidos pelos parentais
      cod.sp[nascer]<-cod.sp[mami] ### Substitui codigos das especies dos mortos pelos codigos dos individuos novos
      n.propag[nascer] <- sapply(1,rtruncnorm,a=(1-0.000001), b=(X+0.000001), mean= medias.prop,sd=dp) ### Sorteia o numero de propagulos produzidos por ciclo pelos novos individuos de uma distribuicao normal discretizada e truncada entre 1 e X e cuja media eh a media dos pais
      medias.p.death <- (p.death[mami]+p.death[papi])/2 ### Calcula a media da probabilidade de morte dos pais
      p.death[nascer] <- sapply(1,rtruncnorm,a=1/200000, b=1, mean=medias.p.death,sd=dp) ### Sorteia a probabilidade de morte dos novos individuos de uma distribuicao normal discretizada e truncada entre um numero muito pequeno e 1 e cuja media eh a media dos pais
    }
    #########################################################################
    ####################### Salvamento dos resultados #######################
    #########################################################################
    if(sum(i==seq(step,ciclo,step))==1) ### Verifica se o ciclo atual eh um dos que devem ser salvos
    {
      ind.mat[,sc] <- cod.sp ### Guarda na posicao sc a identificacao nova dos individuos
      dead.mat[,sc] <- p.death ### Guarda na posicao sc a probabilidade de morte nova dos individuos
      prop.mat[,sc] <- round(n.propag) ### Guarda na posicao sc o numero de propagulos produzidos por ciclo novo dos individuos
      n.dead.vetor[sc] <- n.dead ### Guarda na posicao sc o numero de mortes acumulado ate o ciclo atual
      sc <- sc+1 ### Atualiza o contador que salva os resultados para o proximo ciclo a ser rodado
      ##cat(format(Sys.time(), "%d%b%Y_%H:%M"), "\t ciclo = ", i, "\n") # para avisar a cada ciclo! desligar se estiver usando Rcloud
      #cat(i," ")
    } 
  }
  #############################################################################
  ############################## Termino do ciclo #############################
  #############################################################################
  #cat("...Termino simulacao\n")
  #############################################################################
  ########################### Organizacao do output ###########################
  #############################################################################
  tempo <- seq(0,ciclo,by=step) ### Cria vetor que dara nome as colunas das matrizes 
  colnames(ind.mat) <- tempo ### Nomeia colunas da matriz ind.mat
  colnames(dead.mat) <- tempo ### Nomeia colunas da matriz dead.mat
  colnames(prop.mat) <- tempo ### Nomeia colunas da matriz prop.mat
  names(n.dead.vetor) <- tempo ### Nomeia os elementos do vetor
  resulta=list(tempo=tempo,sp.list=ind.mat,sementes=prop.mat,prob.morte=dead.mat,n.mortes.cumulativo=n.dead.vetor)
  t1=proc.time()[[3]] ### Marca o termino da contagem de tempo de processamento da funcao
  cat("Tempo de processamento: ", round((t1-t0)/60,2),"min\n") ### Mostra o tempo de processamento no console
  attributes(resulta)$start=list(especies=S, individuos=j, nprop=xi0, sd=dp, posicao_disturbios=dist.pos, intensidade_disturbios=dist.int, ciclos=ciclo, passos=step) ### Inclui atributos no objeto resulta
  return(resulta) ### Retorna o objeto resulta
}
#############################################################################
############################# Termino da funcao #############################
#############################################################################


##############################################################################
######### Evolucao de estrategias de vida a partir de cenario neutro #########
###################### SEM TRADE-OFF - DEMONIO DE DARWIN #####################
######################### VERSAO COM CARACTER NEUTRO #########################
##############################################################################
############################## simula.neutra.cn ##############################
##############################################################################
############################ Listagem das versoes ############################
##############################################################################
# 10. Versao modificada por Luisa Novara para incluir pi0 como um argumento (maio 2016)
# 9. Versao modificada por Luisa Novara para incluir X como um argumento e alguns erros de parada (maio 2016)
# 8. Versao modificada por Luisa Novara para tornar a fecundidade (numero de propagulos produzidos por ciclo) um caracter não importante para o fitness (caracter neutro) (maio 2016)
# 7. Versao modificada por Luisa Novara para retirar o trade-off entre a fecundidade e a probabilidade de morte (maio 2016)
# 6. Versao modificada por Luisa Novara para incluir spp com diferentes estrategias de vida iniciais (maio 2016)
# 5. Versao modificada por Luisa Novara para incluir retorno do numero de mortes cumulativo (maio 2016)
# 4. Versao modificada por Luisa Novara e Alexandre Adalardo (fevereiro 2016)
# 3. Disturbio implementado por Alexandre Adalardo e Luisa Novara (novembro 2015) 
# 2. Versao modificada por Luisa Novara (2014)
# 1. Versao inicial de Alexandre Adalardo (outubro 2009) modificada por Paulo Inacio Prado (outubro 2009)
##############################################################################
########################## Detalhamento das versoes ##########################
##############################################################################
# 10. Retira o calculo da probabilidade de morte inicial de 1/J e substitui pela entrada de valores a escolha do usuario (argumento pi0)
# 9. Desliga a linha de comando que calcula X a partir de xi0J e inclui X como um argumento da funcao, para evitar que a estrategia de vida inicial (xi0) esteja sempre muito distante do maximo possivel (X)
# 8. Sorteia os individuos que irao repor os mortos a partir dos individuos adultos da comunidade, e nao a partir do banco de propagulos. Logo, o numero de propagulos que um individuo produz por ciclo nao interfere em sua persistencia na comunidade, e a fecundidade se torna um caracter "neutro" (porque, alem disso, nao existe o trade-off entre fecundidade e sobrevivencia). Isso eh interessante para avaliarmos como ocorre a evolucao de um carater neutro neste modelo
# 7. Troca o calculo da probabilidade de morte de individuos novos de xi0/X (trade-off) para um sorteio de uma distribuicao normal cuja media eh a media dos parentais (assim como eh feito para o numero de propagulos produzidos por ciclo)
# 6. Troca a forma do argumento xi0 de um valor para um vetor de tamanho J, com os valores das estrategias iniciais de cada individuo da comunidade
# 4. Troca o argumento X por xi0, para que X (que eh dado por xi0*J) seja inteiro
# 2. Troca o antigo argumento cv (coeficiente de variacao) pelo dp (desvio padrao) da distribuicao normal da herdabilidade de xi. Isso permite que a herdabilidade permaneça constante e evita o erro de gerar herdabilidade maior quando os valores de xi sao mais baixos (i.e., evita que os valores de xi da populacao fiquem "presos" em valores mais baixos e isso mascare os reais resultados das simulacoes)
##############################################################################
############################## Inicio da funcao ##############################
##############################################################################
simula.neutra.cn=function(S= 100, j=10, xi0=rep(seq(10,10,length.out = S),each=j), X=10000, pi0=rep(seq(0.01,0.01,length.out = S),each=j), dp=0.1, dist.pos=NULL, dist.int=NULL, ciclo=1e6, step=100)
{
  t0=proc.time()[[3]] ### Marca o inicio da contagem de tempo de processamento da funcao
  #cat("Inicio simulacao... Ciclos: ")
  #############################################################################
  ########################### Argumentos deduzidos ############################
  #############################################################################
  J <- S*j ### Calcula o tamanho da comunidade (J)
  #X <- min(xi0)*J ### Calcula o numero total de propagulos produzidos por um individuo (X).
  #############################################################################
  ############################### Verificacoes ################################
  #############################################################################
  if(max(xi0)>X) ### Verifica se o valor/algum dos valores de xi0 (numero de propagulos produzidos por um indivíduo por ciclo) eh maior do que o numero de propagulos produzidos ao longo da vida toda do individuo
  {
    stop("\n\tO número de propágulos produzidos por indivíduo a cada ciclo (xi0) não pode ser maior do que o número de propágulos máximo (X). Tente novamente!\n\n")
  }
  if(sum(dist.int>1)>0) ### Verifica se ha disturbio programado com valor de intensidade que nao faz sentido
  {
    stop("\n\tA intensidade do distúrbio deve ser entre 0 e 1. Tente novamente!\n\n")
  } 
  if(sum(0<=dist.pos & dist.pos<=ciclo)<length(dist.pos)) ### Verifica se ha disturbio programado para ciclo inexistente
  {
    stop("\n\tA posição dos eventos de distúrbio (dist.pos) precisa ser condizente com o número de ciclos a serem rodados (ciclo). Tente novamente!\n\n")
  } 
  if(sum(dist.pos==0 & dist.int>0)>0) ### Verifica se ha disturbio programado para o ciclo 0 (o primeiro ciclo eh o 1)
  {
    stop("\n\tAtenção! O primeiro ciclo das simulações é o ciclo 1. Atribua 0 ao dist.pos apenas quando não desejar implementar a ocorrência de distúrbios.\n\n")
  }
  if(step>ciclo) ### Verifica se o intervalo de ciclos entre os quais os resultados devem ser guardados esta contido no numero de ciclos total da simulacao
  {
    stop("\n\tAtenção! O valor do argumento 'step' não pode ser maior que o de 'ciclo'.\n\n")
  }
  #############################################################################
  ###################### Matrizes para guardar resultados #####################
  #############################################################################
  ind.mat=matrix(nrow=J,ncol=1+ciclo/step) ### Gera matriz da identidade (especie) de cada individuo por ciclo    
  prop.mat=matrix(nrow=J,ncol=1+ciclo/step) ### Gera matriz de propagulos produzidos por individuo em cada ciclo
  dead.mat=matrix(nrow=J,ncol=1+ciclo/step) ### Gera matriz de probabilidade de morte de cada individuo, por ciclo
  ### Vetor com numero de mortes cumulativo ate determinado ciclo
  n.dead <- 0
  n.dead.vetor<-c()
  n.dead.vetor[1] <- 0
  #############################################################################
  ############################# Condicoes iniciais ############################
  #############################################################################
  ### Guarda a identidade inicial de cada individuo
  ind.mat[,1] <- rep(1:S,each=j)
  cod.sp <- ind.mat[,1] ### Transfere informacao para vetor temporario que se atualiza a cada ciclo, depois de ser copiado para uma coluna da matriz
  ### Guarda a probabilidade de morte de cada individuo
  dead.mat[,1] <- pi0
  p.death <- dead.mat[,1] ### Transfere informacao para vetor temporario que se atualiza a cada ciclo, depois de ser copiado para uma coluna da matriz
  ### Guarda o numero de propagulos produzidos por ciclo de cada individuo
  prop.mat[,1] <- xi0 ### Transfere informacao para vetor temporario que se atualiza a cada ciclo, depois de ser copiado para uma coluna da matriz
  n.propag <- prop.mat[,1]
  prop.mat[,1]<-round(prop.mat[,1]) ### Arredonda valor de numero de propagulos produzidos por ciclo na matriz de resultados para que se retorne como saida da funcao o "fenotipo" dos individuos
  #############################################################################
  ############ Contador para salvar resultados a cada step ciclos #############
  #############################################################################
  sc=2 ### Contador que comeca na posicao 2, ja que na primeira posicao foram guardadas as condicoes iniciais
  #############################################################################
  ############################## Inicio do ciclo ##############################
  #############################################################################
  for(i in 1:ciclo)
  {
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
    n.mortes <- sum(morte) ### Grava o numero total de mortes no ciclo
    n.dead <- n.dead + n.mortes ### Atualiza n.dead com numero de mortes acumulado ate entao mais o numero de mortes do ciclo atual
    #########################################################################
    ######################## Substituicao de valores ########################
    #########################################################################
    if(n.mortes>0) ### Identifica se houve mortes no ciclo
    {
      cod.ind<-1:J
      seed.bank <- rep(cod.ind,round(n.propag)) ### Banco de propagulos: cada propagulo tem o codigo numerico do individuo. Como o fenotipo n.propag pode ter valores nao inteiros, arredondamos.
      nascer= which(morte==1) ### Armazena indices dos individuos que morreram
      mami=sample(cod.ind, n.mortes) ### Sorteia os individuos adultos que irao repor os mortos
      papi <- c() ### Cria vetor para armazenar o fenotipo do pai
      for(w in 1:n.mortes) ### Cria loop para sortear o pai entre os individuos da especie de cada propagulo-mae sorteado
      {
        papi[w] <- sample(cod.ind[ cod.sp==cod.sp[mami[w]] ],1)
      }
      medias.prop=(n.propag[mami]+n.propag[papi])/2 ### Calcula o valor esperado de propagulos produzidos por ciclo dos filhotes, representado pela media do numero medio de propagulos produzidos pelos parentais
      cod.sp[nascer]<-cod.sp[mami] ### Substitui codigos das especies dos mortos pelos codigos dos individuos novos
      n.propag[nascer] <- sapply(1,rtruncnorm,a=(1-0.000001), b=(X+0.000001), mean= medias.prop,sd=dp) ### Sorteia o numero de propagulos produzidos por ciclo pelos novos individuos de uma distribuicao normal discretizada e truncada entre 1 e X e cuja media eh a media dos pais
      medias.p.death <- (p.death[mami]+p.death[papi])/2 ### Calcula a media da probabilidade de morte dos pais
      p.death[nascer] <- sapply(1,rtruncnorm,a=1/200000, b=1, mean=medias.p.death,sd=dp) ### Sorteia a probabilidade de morte dos novos individuos de uma distribuicao normal discretizada e truncada entre um numero muito pequeno e 1 e cuja media eh a media dos pais
    }
    #########################################################################
    ####################### Salvamento dos resultados #######################
    #########################################################################
    if(sum(i==seq(step,ciclo,step))==1) ### Verifica se o ciclo atual eh um dos que devem ser salvos
    {
      ind.mat[,sc] <- cod.sp ### Guarda na posicao sc a identificacao nova dos individuos
      dead.mat[,sc] <- p.death ### Guarda na posicao sc a probabilidade de morte nova dos individuos
      prop.mat[,sc] <- round(n.propag) ### Guarda na posicao sc o numero de propagulos produzidos por ciclo novo dos individuos
      n.dead.vetor[sc] <- n.dead ### Guarda na posicao sc o numero de mortes acumulado ate o ciclo atual
      sc <- sc+1 ### Atualiza o contador que salva os resultados para o proximo ciclo a ser rodado
      ##cat(format(Sys.time(), "%d%b%Y_%H:%M"), "\t ciclo = ", i, "\n") # para avisar a cada ciclo! desligar se estiver usando Rcloud
      #cat(i," ")
    } 
  }
  #############################################################################
  ############################## Termino do ciclo #############################
  #############################################################################
  #cat("...Termino simulacao\n")
  #############################################################################
  ########################### Organizacao do output ###########################
  #############################################################################
  tempo <- seq(0,ciclo,by=step) ### Cria vetor que dara nome as colunas das matrizes 
  colnames(ind.mat) <- tempo ### Nomeia colunas da matriz ind.mat
  colnames(dead.mat) <- tempo ### Nomeia colunas da matriz dead.mat
  colnames(prop.mat) <- tempo ### Nomeia colunas da matriz prop.mat
  names(n.dead.vetor) <- tempo ### Nomeia os elementos do vetor
  resulta=list(tempo=tempo,sp.list=ind.mat,sementes=prop.mat,prob.morte=dead.mat,n.mortes.cumulativo=n.dead.vetor)
  t1=proc.time()[[3]] ### Marca o termino da contagem de tempo de processamento da funcao
  cat("Tempo de processamento: ", round((t1-t0)/60,2),"min\n") ### Mostra o tempo de processamento no console
  attributes(resulta)$start=list(especies=S, individuos=j, nprop=xi0, nprop_max=X, prob_mort=pi0, sd=dp, posicao_disturbios=dist.pos, intensidade_disturbios=dist.int, ciclos=ciclo, passos=step) ### Inclui atributos no objeto resulta
  return(resulta) ### Retorna o objeto resulta
}
#############################################################################
############################# Termino da funcao #############################
#############################################################################
