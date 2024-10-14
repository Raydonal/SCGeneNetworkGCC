#################################################################################################
######          Simulacao do coeficiente de correlacao da normal bivariada                 ######
######            Comparando com o coeficiente de correlacao de Chinchilli                 ######
######                                                                                     ######
#################################################################################################

#### Caso 3 : N(-0.5,0.5;1,1;rho)

#### Escolhendo Gerador #####
rm(list=ls(all=TRUE)) 
RNGkind("Mersenne-Twister")

library(mvtnorm) ## gerar numeros normais e t
library(BMS) ## funcao hypergeometric

###### Variaveis
NOBS_vetor=c(10,50,100,250,500) ## Tamanho das amostras 10,50,100,250,500
NREP=5000
mean=matrix(c(-0.5,0.5),ncol=2)
sigma1=1
sigma2=1
rho_vetor=c(0,0.3,0.9)

## Matrizes de estimativas
EMV=matrix(0,1,NREP) ## Estimativas de maxima verossimilhança
media_rho_mc=matrix(0,3,5) ## matriz de estimativas
var_rho_mc=matrix(0,3,5) ## matriz de variancias
vies_rho_mc=matrix(0,3,5) ## matriz do Vies do coeficiente de correlacao
eqm_rho_mc=matrix(0,3,5) ## matriz do Erro quadratico medio
raiz_eqm_rho_mc=matrix(0,3,5) ## Raiz quadrada do erro quadrático médio
ptm=proc.time() ## Tempo

## Diretorio da função de chinchilli
setwd("/home/cleber/Dropbox/Dissertacao-Cleber/Simulacao/R/Simulacao 17⁄11⁄2013/Simulacao 18⁄11⁄2013/Caso3/")

## Matrizes das estimativas de Chinchilli
Rho_Gama=matrix(0,3,5)
Gamma=c(0,0.5,1) ## Valores de Gamma
media_rho_ch=matrix(0,3,5)
var_rho_ch=matrix(0,3,5)
vies_rho_ch=matrix(0,3,5)
eqm_rho_ch=matrix(0,3,5)
raiz_eqm_rho_ch=matrix(0,3,5)
dyn.load("somac.so")

## Matrizes de estimativas de Spearman
ESP=matrix(0,1,NREP) ## Estimativas de maxima verossimilhança
media_rho_sp=matrix(0,3,5) ## matriz de estimativas
var_rho_sp=matrix(0,3,5) ## matriz de variancias
vies_rho_sp=matrix(0,3,5) ## matriz do Vies do coeficiente de correlacao
eqm_rho_sp=matrix(0,3,5) ## matriz do Erro quadratico medio
raiz_eqm_rho_sp=matrix(0,3,5) ## Raiz quadrada do erro quadrático médio

## Matriz da distribuicao de Chinchilli
EMVc=matrix(0,3,NREP)

## Funcao de log-verossimilhanca do coeficiente de correlacao
log_vero=function(x){
  rho_hat=x[1]
  c=0
  (-1)*(c-NOBS*(log(sqrt(sigma1))+log(sqrt(sigma2))+1/2*log(1-rho_hat^2)+mean[,1]^2/(2*(1-rho_hat^2)*sigma1)
                +mean[,2]^2/(2*(1-rho_hat^2)*sigma2)-(rho_hat*mean[,1]*mean[,2])/((1-rho_hat^2)*sigma1*sigma2))
        -sum(amostra_biv[,1]^2)/(2*(1-rho_hat^2)*sigma1)-sum(amostra_biv[,2]^2)/(2*(1-rho_hat^2)*sigma2)
        +(mean[,1]*sigma2-mean[,2]*sigma1*rho_hat)/(sigma1*sqrt(sigma2)*(1-rho_hat^2))*sum(amostra_biv[,1])
        +(mean[,2]*sigma1-mean[,1]*sigma2*rho_hat)/(sigma2*sqrt(sigma1)*(1-rho_hat^2))*sum(amostra_biv[,2])
        +rho_hat/(sqrt(sigma1*sigma2)*(1-rho_hat^2))*sum(amostra_biv[,1]*amostra_biv[,2]))
}


######### Funcao de Chinchilli #########
chinc=function(x,y,n,gama){
  dens=.C("chinc",as.double(x),as.double(y),as.integer(n),as.double(gama),resultado=as.double(1))
  dens[["resultado"]]
}

######### Funcao distribuicao de Chinchilli

hat_rho_gama=function(gama,x){
  rho_hat=x[1]
  2*pi^(-1/2)*(gamma(1/2*gama+1))^2/(gamma(gama+1/2))*rho_hat*f21hyper(1/2*(1-gama),1/2*(1-gama),3/2,rho_hat^2)
}

##### Inicio do Monte Carlo ######
for(m in 1:3){
  rho=rho_vetor[m]
  cov=rho*sqrt(sigma1*sigma2) # definido no artigo 
  matrix_cov=matrix(c(sigma1,cov,
                      cov,sigma2),ncol=2) # matriz de variancia e covariancia
  for (j in 1:3){
    gama=Gamma[j]
    rho_gama=hat_rho_gama(gama,rho)
    for(k in 1:5){
      NOBS=NOBS_vetor[k]
      
      set.seed(2) ## Semente
      
      for(i in 1:NREP){
        amostra_biv=rmvnorm(NOBS,mean,matrix_cov)
        ##### Estimacao sem o uso do gradiente
        ir=optim(0.1,log_vero,NULL,method="BFGS")
        EMV[1,i]=ir$par
        EMVc[1,i]=hat_rho_gama(gama,EMV[1,i]) # estimativas de maxima verossimilhanca na distribuicao
        
        cor_chinchilli=chinc(amostra_biv[,1],amostra_biv[,2],NOBS,gama)
        EMVc[2,i]=cor_chinchilli
        
        ESP[1,i]= cor(amostra_biv[,1],amostra_biv[,2],method=c("spearman"))
        EMVc[3,i]=hat_rho_gama(gama,ESP[1,i])
      }
      
      if (k==1) ENREP1=round(EMVc,4) ## Estimativas para o número de replicas com  n=10
      if (k==2) ENREP2=round(EMVc,4) ## Estimativas para o número de replicas com  n=50
      if (k==3) ENREP3=round(EMVc,4) ## Estimativas para o número de replicas com  n=100
      if (k==4) ENREP4=round(EMVc,4) ## Estimativas para o número de replicas com  n=250
      if (k==5) ENREP5=round(EMVc,4) ## Estimativas para o número de replicas com  n=500       
      
      #### Media das estimativas
      media_rho_mc[j,k]=mean(EMVc[1,]) # estimativas de maxima verossimilhança
      media_rho_ch[j,k]=mean(EMVc[2,]) # estimativas de chinchilli
      media_rho_sp[j,k]=mean(EMVc[3,]) # estimativas de Spearman
      #### Variancia das estimativas 
      var_rho_mc[j,k]= var(EMVc[1,]) #sum(EMVc[1,]-rho_gama)^2/(NOBS)
      var_rho_ch[j,k]= var(EMVc[2,]) #sum(EMVc[2,]-rho_gama)^2/(NOBS)
      var_rho_sp[j,k]= var(EMVc[3,]) #sum(EMVc[3,]-rho_gama)^2/(NOBS)
      #### Vies do coeficiente de correlacao
      vies_rho_mc[j,k]=mean(EMVc[1, ])-rho_gama
      vies_rho_ch[j,k]=mean(EMVc[2, ])-rho_gama
      vies_rho_sp[j,k]=mean(EMVc[3, ])-rho_gama
      #### Erro quadratico medio
      eqm_rho_mc[j,k]=vies_rho_mc[j,k]^2 + var_rho_mc[j,k]
      eqm_rho_ch[j,k]=vies_rho_ch[j,k]^2 + var_rho_ch[j,k]
      eqm_rho_sp[j,k]=vies_rho_sp[j,k]^2 + var_rho_sp[j,k]
      #### Raiz quadrada do erro quadrático médio
      raiz_eqm_rho_mc[j,k]=round(sqrt(eqm_rho_mc[j,k]),4)
      raiz_eqm_rho_ch[j,k]=round(sqrt(eqm_rho_ch[j,k]),4)
      raiz_eqm_rho_sp[j,k]=round(sqrt(eqm_rho_sp[j,k]),4)
      
      ####matriz com os valores de rho_gamma
      Rho_Gama[j,]=rho_gama
    }
    if (j==1) dados1=cbind(ENREP1,ENREP2,ENREP3,ENREP4,ENREP5) ## Matrix com os varios n para gama=0
    if (j==2) dados2=cbind(ENREP1,ENREP2,ENREP3,ENREP4,ENREP5) ## Matrix com os varios n para gama=0.5
    if (j==3) dados3=cbind(ENREP1,ENREP2,ENREP3,ENREP4,ENREP5) ## Matrix com os varios n para gama=1
    
  }
  #############  Impressao dos resultados  #################
  
  if (m==1){
    
    save.image("Rho0Mistura.Rdata")
    ## Saindo tabela
    write.table(round(rbind(Rho_Gama,media_rho_mc,media_rho_ch,media_rho_sp),4),file="EstimativaSemPerturbacaoRho=0.csv"
                ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
    write.table(round(rbind(vies_rho_mc,vies_rho_ch,vies_rho_sp),4),file="ViésSemPerturbacaoRho=0.csv"
                ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
    write.table(round(rbind(var_rho_mc,var_rho_ch,var_rho_sp),4),file="VarianciaSemPerturbacaoRho=0.csv"
                ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
    write.table(round(rbind(raiz_eqm_rho_mc,raiz_eqm_rho_ch,raiz_eqm_rho_sp),4),file="ReqmSemPerturbacaoRho=0.csv"
                ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))}
  
  if (m==2){
    
    save.image("Rho03Mistura.Rdata")
    ## Saindo tabela
    write.table(round(rbind(Rho_Gama,media_rho_mc,media_rho_ch,media_rho_sp),4),file="EstimativaSemPerturbacaoRho=03.csv"
                ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
    write.table(round(rbind(vies_rho_mc,vies_rho_ch,vies_rho_sp),4),file="ViésSemPerturbacaoRho=03.csv"
                ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
    write.table(round(rbind(var_rho_mc,var_rho_ch,var_rho_sp),4),file="VarianciaSemPerturbacaoRho=03.csv"
                ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
    write.table(round(rbind(raiz_eqm_rho_mc,raiz_eqm_rho_ch,raiz_eqm_rho_sp),4),file="ReqmSemPerturbacaoRho=03.csv"
                ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))}
  
  if (m==3){
    
    save.image("Rho09Mistura.Rdata")
    ## Saindo tabela
    write.table(round(rbind(Rho_Gama,media_rho_mc,media_rho_ch,media_rho_sp),4),file="EstimativaSemPerturbacaoRho=09.csv"
                ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
    write.table(round(rbind(vies_rho_mc,vies_rho_ch,vies_rho_sp),4),file="ViésSemPerturbacaoRho=09.csv"
                ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
    write.table(round(rbind(var_rho_mc,var_rho_ch,var_rho_sp),4),file="VarianciaSemPerturbacaoRho=09.csv"
                ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
    write.table(round(rbind(raiz_eqm_rho_mc,raiz_eqm_rho_ch,raiz_eqm_rho_sp),4),file="ReqmSemPerturbacaoRho=09.csv"
                ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))}
  
  
  ## Tempo
  proc.time()-ptm
}

#### Gerando Gráficos Boxplot
load("~/Dropbox/Dissertacao-Cleber/Simulacao/R/Simulacao 17⁄11⁄2013/Simulacao 18⁄11⁄2013/Caso3/Rho09Mistura.Rdata")
setwd("/home/cleber/Dropbox/Dissertacao-Cleber/Simulacao/Graficos/colorido/Simulacao 17⁄11⁄2013/Simulacao 21⁄11⁄2013/Caso3/Rho09/")

library(reshape2)
library(ggplot2)

## Resultados para Gamma=0
dados1_trans=t(dados1)
dados1_rep=cbind(dados1_trans,rep(c(10,50,100,250,500),each=NREP))
dados1_rep=as.data.frame(dados1_rep)
names(dados1_rep)=c("r1","r2","r3","n")

dat.m <- melt(dados1_rep,id.vars='n', measure.vars=c('r1','r2','r3'))
dat.m$n=factor(dat.m$n)
labl <- list(expression(hat(rho)[gamma]),expression(tilde(rho)[gamma]),expression(bar(rho)[S][gamma]))

ggplot(dat.m) +  theme_bw() + xlab("Tamanho da amostra")+ ylab("Valores")+
  geom_boxplot(aes(x=n, y=value, fill=variable),position = position_dodge(width = .9))+geom_abline(intercept=Rho_Gama[1,1],slope=0)+
  guides(fill = guide_legend(keywidth = 2, keyheight = 2,title="Coeficientes"))+scale_fill_grey(start=0.5,end=1,labels=labl)


## Salvando os graficos em eps e pdf 
dev.copy2pdf(file="Grafico1(Boxplot).pdf")
dev.copy2eps(file="Grafico1(Boxplot).eps")

### Resultados para Gammma=0.5
dados2_trans=t(dados2)
dados2_rep=cbind(dados2_trans,rep(c(10,50,100,250,500),each=NREP))
dados2_rep=as.data.frame(dados2_rep)
names(dados2_rep)=c("r1","r2","r3","n")

dat2.m <- melt(dados2_rep,id.vars='n', measure.vars=c('r1','r2','r3'))
dat2.m$n=factor(dat2.m$n)
labl <- list(expression(hat(rho)[gamma]),expression(tilde(rho)[gamma]),expression(bar(rho)[S][gamma]))

ggplot(dat2.m) +  theme_bw() + xlab("Tamanho da amostra")+ ylab("Valores")+
  geom_boxplot(aes(x=n, y=value, fill=variable),position = position_dodge(width = .9))+geom_abline(intercept=Rho_Gama[2,1],slope=0)+
  guides(fill = guide_legend(keywidth = 2, keyheight = 2,title="Coeficientes"))+scale_fill_grey(start=0.5,end=1,labels=labl)
## Salvando os graficos em eps e pdf 
dev.copy2pdf(file="Grafico2(Boxplot).pdf")
dev.copy2eps(file="Grafico2(Boxplot).eps")

### Resultados para Gama=1
dados3_trans=t(dados3)
dados3_rep=cbind(dados3_trans,rep(c(10,50,100,250,500),each=NREP))
dados3_rep=as.data.frame(dados3_rep)
names(dados3_rep)=c("r1","r2","r3","n")

dat3.m <- melt(dados3_rep,id.vars='n', measure.vars=c('r1','r2','r3'))
dat3.m$n=factor(dat3.m$n)
labl <- list(expression(hat(rho)[gamma]),expression(tilde(rho)[gamma]),expression(bar(rho)[S][gamma]))

ggplot(dat3.m) +  theme_bw() + xlab("Tamanho da amostra")+ ylab("Valores")+
  geom_boxplot(aes(x=n, y=value, fill=variable),position = position_dodge(width = .9))+geom_abline(intercept=Rho_Gama[3,1],slope=0)+
  guides(fill = guide_legend(keywidth = 2, keyheight = 2,title="Coeficientes"))+scale_fill_grey(start=0.5,end=1,labels=labl)

## Salvando os graficos em eps e pdf 
dev.copy2pdf(file="Grafico3(Boxplot).pdf")
dev.copy2eps(file="Grafico3(Boxplot).eps")
