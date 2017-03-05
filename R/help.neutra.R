####################################################################
########## FUNCOES ACESSORIAS E GRAFICAS DO MODELO NEUTRO ##########
####################################################################

###### OBSERVACAO 1: AQUI, simulacao = RESULTADO (LISTA FINAL) DE DETERMINADA SIMULACAO COM simula.neutra.step

####################################################################
########### Grafico da variacao na SAD ao longo do tempo ###########
########################## AAO junho 2010 ##########################
####################################################################

graf.abund=function(dados=resulta,...)
{
  info=attributes(dados)
  n.prop=info$start$nprop
  cv=info$start$cv
  tempo=dados$tempo
  ntempo=length(tempo)
  nspp=length(unique(dados$sp.list[,1]))
  #nmax=max(unlist(dados))
  nmax=max(table(dados$sp.list[,dim(dados$sp.list)[2]]))
  plot(x=c(1,nspp),y=c(1, nmax),log="y", type="n", ylab= "Abundance", xlab="Rank order", cex.lab=1.2, cex.axis=1.2, sub= paste("total seeds=", n.prop, "\t cv = ", cv),... )
  stempos=round(seq(2, (ntempo-100),length.out=100 ))
  colors=rainbow(length(stempos)*10)
  ncol=length(colors)
  ncol.meio=round(ncol*0.5)
  i=ncol.meio
  for (t in stempos )
  {
    i=i+1
    points(sort(table(factor(dados$sp.list[,t], levels=1:nspp)), decreasing=TRUE),type="l", lwd=0.5, col=colors[i])
  }
  points(sort(table(factor(dados$sp.list[,1], levels=1:nspp)), decreasing=TRUE), type="l", col="green", lwd=2)
  points(sort(table(factor(dados$sp.list[,ncol.meio], levels=1:nspp)), decreasing=TRUE), type="l", col="blue", lwd=2)
  points(sort(table(factor(dados$sp.list[,dim(dados$sp.list)[2]], levels=1:nspp)), decreasing=TRUE), type="l", col="red", lwd=2)
  legend("topright", lty=1, col=c("green", "blue", "red"), bty="n", legend=c("start", "middle", "end") )
}

####################################################################
######################## Rank-abundance plot ####################### 
###  grafico simples de Whitaker (log abundâncias x ranque) ########
############# entrada: vetor de abundancia de especies ############# 
########################## AAO junho 2010 ##########################
####################################################################

w.plot <- function(x,...)
{
  plot(1:length(x),sort(x,decreasing=T),
       ylab="Abundance",log="y",...)
}

####################################################################
################## Box-plot e rank-abundance plot ################## explora graficamente a distribuicao da fertilidade. Plota box-plots da fertlidade dos individuos de cada especie, sobre um grafico de rank-abundancia da fertilidade de cada especie, num dado instante gravado da simulacoo (argumento tempo)
###### entrada: dados = objeto (lista resultante) da simulacao #####  
########################## AAO junho 2010 ##########################
####################################################################

box.w <- function(dados, tempo=NULL )
{
  if(is.null(tempo))
  {tempo=dim(dados$sp.list)[2]}
  #tempo,cod.sp,n.propag, info=attributes(dados)
  ast=attributes(dados)$start
  sp <- dados$sp.list[,tempo]
  prop <- dados$sementes[,tempo]
  ab <- tapply(sp,sp,length)
  par(mfrow=c(2,1), mar=c(2,4,4,3), mgp=c(2,1,0))
  w.plot(ab, xlab="", xaxt="n", main=paste("NEUTRALITY SIMULATION\n", paste(names(ast), ast, sep="=", collapse="; "), "\n(time = ",tempo, ")"), cex.main=0.8)
  par(mar=c(5,4,0,3))#, new=TRUE, oma=c(0,2,2,2))
  boxplot(prop~factor(sp,levels=names(ab)[order(ab,decreasing=T)]),xlab="Species Codes", ylab="Number of Seeds")
  par(mfrow=c(1,1))
}

####################################################################
############## Fertilidade media por especie + Grafico ############# calcula uma estatistica descritiva (o default eh a media) da fertilidade de cada especie, em cada instante de tempo. Fora da funcao, pode-se plotar a matriz resultante em um grafico
############# entrada: cod.sp = $sp.list da simulacao, #############
################# n.propag = $sementes da simulacao ################
########################## AAO junho 2010 ##########################
####################################################################
#### versao da fert.t sem loop -renomeada, original fert.t1
######################### AAO fevereiro 2012 #######################

fert.t <- function(cod.sp,n.propag,fun=mean)
{ 
  especie <- unique(cod.sp[,1])
  sp.level<-factor(cod.sp,levels=cod.sp)
  t.a<-function(x){tapply(n.propag[,x],factor(cod.sp[,x],levels=especie),fun)}
  res<-sapply(1:ncol(n.propag), t.a)
  #names(res)<-colnames(n.propag)
  colnames(res) <- colnames(n.propag)
  rownames(res) <- paste("sp",especie, sep="")
  return(res)
}

#################################################################### 
######################### Input: hipercubo #########################
############### entrada: input_hipercubo= data.frame ###############
# Jm= tam comunidades, ciclos= numero ciclos, steps= numero passos #
# Obs: diferente de simula_externa_input, funciona diretamente e apenas com simula.neutra.step, e nao com simula_externa_output 
######################### NOVARA, MAR 2016 #########################

simula_input<-function(input_hipercubo,Jm=5000,ciclos=100000,steps=100){
  if(dim(input_hipercubo)[2]<5){stop("Deve ser atribuído valor de entrada para todos os argumentos de fun. Tente novamente!")}
  # input_hipercubo deve ser um data.frame, em que cada coluna eh um parametro do modelo simula.neutra.step (logo, cada linha eh um conjunto de valores de parametros que sera usado em uma simulacao)
  # Jm eh o numero medio de individuos da comunidade, varia dependendo do arredondamento de j e posterior recalculo de J (dentro de simula.neutra.step)
  input_hipercubo$j<-round(Jm/input_hipercubo[,1])
  input_hipercubo$ciclo<-ciclos
  input_hipercubo$step<-steps
  input<-cbind(input_hipercubo[,1],input_hipercubo[,6],input_hipercubo[,2],input_hipercubo[,3],input_hipercubo[,4],input_hipercubo[,5],input_hipercubo[,7],input_hipercubo[,8])
  res<-list()
  for(i in 1:(dim(input)[1]))
  {
    res[[i]]<-simula.neutra.step(S=input[i,1],
                                 j=input[i,2],
                                 xi0=input[i,3],
                                 dp=input[i,4],
                                 dist.pos=seq(input[i,5],ciclos,input[i,5]),
                                 dist.int=input[i,6],
                                 ciclo=input[i,7],
                                 step=input[i,8])
  }
  return(res) 
}

#################################################################### 
######################### Output reduzido ##########################
######### entrada: resultado das simulacoes com simula.neutra.step (ou simula_input) #########
######### saida: media, variancia, assimetria e curtose do #########
######## numero de propagulos/ciclo (proxy para estrategia) ########
######################### NOVARA, MAR 2016 #########################

require(PerformanceAnalytics)

simula_output<-function(lista_simulacoes,nger=5000)
{
  x<-lista_simulacoes
  resultado<-matrix()
  resultado_lista<-list()
  for (i in 1:length(x))
  {
    tam_com <- attributes(x[[i]])$start[[1]]*attributes(x[[i]])$start[[2]]
    mortes_acumuladas <- x[[i]]$n.mortes.cumulativo
    n.mortes_prox_nger <- mortes_acumuladas - (nger*tam_com)
    ngeracao <- which(abs(n.mortes_prox_nger)==min(abs(n.mortes_prox_nger)))
    prop <- as.matrix(x[[i]]$sementes)
    sp <- as.matrix(x[[i]]$sp.list)
    abund<-as.vector(table(sp[,ngeracao]))
    riq.fin <- length(unique(sp[,ngeracao]))
    riq.inic <- length(unique(sp[,1]))
    med_geral <- mean(prop[,ngeracao])
    med_sp <- tapply(X=prop[,ngeracao],INDEX=sp[,ngeracao],FUN=mean)
    median_geral <- median(prop[,ngeracao])
    median_sp <- tapply(X=prop[,ngeracao],INDEX=sp[,ngeracao],FUN=median)
    moda_geral <- prop[,ngeracao][max(table(prop[,ngeracao]))]
    moda_sp <- tapply(X=prop[,ngeracao],INDEX=sp[,ngeracao],FUN=function(coisa){return(coisa[max(table(coisa))])})
    vari_geral <- var(prop[,ngeracao])
    vari_sp <- tapply(X=prop[,ngeracao],INDEX=sp[,ngeracao],FUN=var)
    assim_geral <- skewness(prop[,ngeracao])
    assim_sp <- tapply(X=prop[,ngeracao],INDEX=sp[,ngeracao],FUN=skewness)
    exc_curt_geral <- kurtosis(prop[,ngeracao])
    exc_curt_sp <- tapply(X=prop[,ngeracao],INDEX=sp[,ngeracao],FUN=kurtosis)
    resultado <- matrix(data=c(riq.inic,rep(NA,riq.fin),riq.fin,rep(NA,riq.fin),tam_com,abund,med_geral,med_sp,median_geral,median_sp,moda_geral,moda_sp,vari_geral,vari_sp,assim_geral,assim_sp,exc_curt_geral,exc_curt_sp),ncol=9,dimnames=list(c("geral",paste("sp",(unique(sp[,ngeracao])[order(unique(sp[,ngeracao]))]),sep="")),c("riq_inic","riq_fin","abundancia","media","mediana","moda","variancia","assimetria","excesso_curtose")))
    resultado_lista[[i]]<-resultado
  }
  return(resultado_lista)
}