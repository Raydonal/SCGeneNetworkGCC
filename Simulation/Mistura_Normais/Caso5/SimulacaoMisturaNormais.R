#################################################################################################
######          Simulacao do coeficiente de correlacao da normal bivariada                 ######
######            Comparando com o coeficiente de correlacao de Chinchilli                 ######
######                          Mistura de Normais                                         ######
######      N(0,0;1,1;0.9)           /       N(0,0;1,1;0)      ######
#################################################################################################

############### Caso 5 0.2 x N_1 + 0.8 x N_2

#### Escolhendo Gerador #####
rm(list=ls(all=TRUE)) 
RNGkind("Mersenne-Twister")

library(mvtnorm) ## gerar numeros normais e t
library(BMS) ## funcao hypergeometric

###### Variaveis
NOBS_vetor=c(10,50,100,250,500) ## Tamanho das amostras 10,50,100,250,500
NREP=5000
p=0.8

##### Definindo a normal padrão 
mean=matrix(c(0,0),ncol=2)
sigma1=sigma2=1
rho1=0.9

##### Definindo a outra normal variando média e variância
mean2=matrix(c(0,0),ncol=2)
sigma3=1
sigma4=1
rho2=0

## Matrizes de estimativas
EMV=matrix(0,1,NREP) ## Estimativas de maxima verossimilhança
media_rho_mc=matrix(0,4,5) ## matriz de estimativas
var_rho_mc=matrix(0,4,5) ## matriz de variancias
vies_rho_mc=matrix(0,4,5) ## matriz do Vies do coeficiente de correlacao
eqm_rho_mc=matrix(0,4,5) ## matriz do Erro quadratico medio
raiz_eqm_rho_mc=matrix(0,4,5) ## Raiz quadrada do erro quadrático médio
ptm=proc.time() ## Tempo

## Diretorio da função de chinchilli
setwd("~/Dropbox/Dissertacao-Cleber/Simulacao/R/NovaSimulacaoMistura/MisturaNormais/Caso1/")

## Matrizes das estimativas de Chinchilli
Rho_Gama=matrix(0,4,5)
Gamma=c(0,0.35,0.65,1) ## Valores de Gamma
media_rho_ch=matrix(0,4,5)
var_rho_ch=matrix(0,4,5)
vies_rho_ch=matrix(0,4,5)
eqm_rho_ch=matrix(0,4,5)
raiz_eqm_rho_ch=matrix(0,4,5)
dyn.load("somac.so")

## Matrizes de estimativas de Spearman
ESP=matrix(0,1,NREP) ## Estimativas de maxima verossimilhança
media_rho_sp=matrix(0,4,5) ## matriz de estimativas
var_rho_sp=matrix(0,4,5) ## matriz de variancias
vies_rho_sp=matrix(0,4,5) ## matriz do Vies do coeficiente de correlacao
eqm_rho_sp=matrix(0,4,5) ## matriz do Erro quadratico medio
raiz_eqm_rho_sp=matrix(0,4,5) ## Raiz quadrada do erro quadrático médio

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
hat_rho_gama1=function(gama,x){
  rho_hat=x[1]
  2*pi^(-1/2)*(gamma(1/2*gama+1))^2/(gamma(gama+1/2))*rho_hat*f21hyper(1/2*(1-gama),1/2*(1-gama),3/2,rho_hat^2)
}

## Funcao para calcular rho_gama quando temos mistura de normais
rho_gama2=function(gama,Rho1,Rho2){
  2*pi^(-1/2)*(gamma(1/2*gama+1))^2/(gamma(gama+1/2))*((1-p)*Rho1*(1-Rho1^2)^(gama+1/2)*f21hyper(1/2*gama+1,1/2*gama+1,3/2,Rho1^2)+
                                                         p*Rho2*(1-Rho2^2)^(gama+1/2)*f21hyper(1/2*gama+1,1/2*gama+1,3/2,Rho2^2))
}


#### Definindo a primeira normal 
cov=rho1*sqrt(sigma1*sigma2) # definido no artigo 
matrix_cov1=matrix(c(sigma1,cov,
                     cov,sigma2),ncol=2) 

#### Definindo a outra normal 
cov2=rho2*sqrt(sigma3*sigma4) # definido no artigo 
matrix_cov2=matrix(c(sigma3,cov2,
                     cov2,sigma4),ncol=2)

##### Inicio do Monte Carlo ######

for (j in 1:4){
  gama=Gamma[j]
  rho_gama=rho_gama2(gama,rho1,rho2)
  for(k in 1:5){
    NOBS=NOBS_vetor[k]
    
    set.seed(2) ## Semente
    
    for(i in 1:NREP){
      amostra_normal1=rmvnorm(NOBS,mean,matrix_cov1)
      amostra_normal2=rmvnorm(NOBS,mean2,matrix_cov2)
      amostra_biv=(1-p)*amostra_normal1+p*amostra_normal2
      
      ##### Estimacao sem o uso do gradiente
      ir=optim(0.1,log_vero,NULL,method="BFGS")
      EMV[1,i]=ir$par
      EMVc[1,i]=hat_rho_gama1(gama,EMV[1,i]) # estimativas de maxima verossimilhanca na distribuicao
      
      cor_chinchilli=chinc(amostra_biv[,1],amostra_biv[,2],NOBS,gama)
      EMVc[2,i]=cor_chinchilli
      
      ESP[1,i]= cor(amostra_biv[,1],amostra_biv[,2],method=c("spearman"))
      EMVc[3,i]=hat_rho_gama1(gama,ESP[1,i])
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
    var_rho_mc[j,k]=var(EMVc[1,])
    var_rho_ch[j,k]=var(EMVc[2,])
    var_rho_sp[j,k]=var(EMVc[3,])
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
  if (j==2) dados2=cbind(ENREP1,ENREP2,ENREP3,ENREP4,ENREP5) ## Matrix com os varios n para gama=0.35
  if (j==3) dados3=cbind(ENREP1,ENREP2,ENREP3,ENREP4,ENREP5) ## Matrix com os varios n para gama=0.65
  if (j==4) dados4=cbind(ENREP1,ENREP2,ENREP3,ENREP4,ENREP5) ## Matrix com os varios n para gama=1
}
## Tempo
proc.time()-ptm

setwd("/home/cleber/Dropbox/Dissertacao-Cleber/Simulacao/R/Simulacao 17⁄11⁄2013/MisturaNormais/Caso5/")

#############  Impressao dos resultados  #################

## Saindo tabela
write.table(round(rbind(Rho_Gama,media_rho_mc,media_rho_ch,media_rho_sp),4),file="EstimativaMisturaCaso5.csv"
            ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
write.table(round(rbind(vies_rho_mc,vies_rho_ch,vies_rho_sp),4),file="ViésMisturaCaso5.csv"
            ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
write.table(round(rbind(var_rho_mc,var_rho_ch,var_rho_sp),4),file="VarianciaMisturaCaso5.csv"
            ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
write.table(round(rbind(raiz_eqm_rho_mc,raiz_eqm_rho_ch,raiz_eqm_rho_sp),4),file="ReqmMisturaCaso5.csv"
            ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))


### Gerando os graficos boxplot
setwd("/home/cleber/Dropbox/Dissertacao-Cleber/Simulacao/Graficos/colorido/Simulacao 17⁄11⁄2013/MisturaNormais/Caso5/")  
save.image("RhoMistura.Rdata")
load("~/Dropbox/Dissertacao-Cleber/Simulacao/Graficos/colorido/Simulacao 17⁄11⁄2013/MisturaNormais/Caso5/RhoMistura.Rdata")

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
hline <- data.frame(yint = Rho_Gama[1,1],lt = 'Valor Real')
labl2 <- list(expression(rho[gamma]))

ggplot(dat.m) +  theme_bw() + xlab("Tamanho da amostra")+ ylab("Valores")+
  geom_boxplot(aes(x=n, y=value, fill=variable),position = position_dodge(width = .9))+
  guides(fill = guide_legend(title="Coeficientes"))+scale_fill_grey(start=0.5,end=1,labels=labl)+
  theme(legend.text = element_text(size = 20),legend.title = element_text(size=15),axis.title=element_text(size=20),axis.text=element_text(size=20))+
  geom_hline(data = hline,aes(yintercept=yint,linetype = lt),color = "Black",show_guide = TRUE) + 
  scale_colour_discrete(guide = "none") + 
  scale_linetype_manual(name = 'Verdadeiro',values = 1,guide = "legend",labels=labl2)
## Salvando os graficos em eps e pdf 
dev.copy2pdf(file="Grafico1(Boxplot).pdf")
dev.copy2eps(file="Grafico1(Boxplot).eps")

### Resultados para Gammma=0.35
dados2_trans=t(dados2)
dados2_rep=cbind(dados2_trans,rep(c(10,50,100,250,500),each=NREP))
dados2_rep=as.data.frame(dados2_rep)
names(dados2_rep)=c("r1","r2","r3","n")

dat2.m <- melt(dados2_rep,id.vars='n', measure.vars=c('r1','r2','r3'))
dat2.m$n=factor(dat2.m$n)
labl <- list(expression(hat(rho)[gamma]),expression(tilde(rho)[gamma]),expression(bar(rho)[S][gamma]))
hline <- data.frame(yint = Rho_Gama[2,1],lt = 'Valor Real')
labl2 <- list(expression(rho[gamma]))

ggplot(dat2.m) +  theme_bw() + xlab("Tamanho da amostra")+ ylab("Valores")+
  geom_boxplot(aes(x=n, y=value, fill=variable),position = position_dodge(width = .9))+
  guides(fill = guide_legend(title="Coeficientes"))+scale_fill_grey(start=0.5,end=1,labels=labl)+
  theme(legend.text = element_text(size = 20),legend.title = element_text(size=15),axis.title=element_text(size=20),axis.text=element_text(size=20))+
  geom_hline(data = hline,aes(yintercept=yint,linetype = lt),color = "Black",show_guide = TRUE) + 
  scale_colour_discrete(guide = "none") + 
  scale_linetype_manual(name = 'Verdadeiro',values = 1,guide = "legend",labels=labl2)
## Salvando os graficos em eps e pdf 
dev.copy2pdf(file="Grafico2(Boxplot).pdf")
dev.copy2eps(file="Grafico2(Boxplot).eps")

### Resultados para Gama=0.65
dados3_trans=t(dados3)
dados3_rep=cbind(dados3_trans,rep(c(10,50,100,250,500),each=NREP))
dados3_rep=as.data.frame(dados3_rep)
names(dados3_rep)=c("r1","r2","r3","n")

dat3.m <- melt(dados3_rep,id.vars='n', measure.vars=c('r1','r2','r3'))
dat3.m$n=factor(dat3.m$n)
labl <- list(expression(hat(rho)[gamma]),expression(tilde(rho)[gamma]),expression(bar(rho)[S][gamma]))
hline <- data.frame(yint = Rho_Gama[3,1],lt = 'Valor Real')
labl2 <- list(expression(rho[gamma]))

ggplot(dat3.m) +  theme_bw() + xlab("Tamanho da amostra")+ ylab("Valores")+
  geom_boxplot(aes(x=n, y=value, fill=variable),position = position_dodge(width = .9))+
  guides(fill = guide_legend(title="Coeficientes"))+scale_fill_grey(start=0.5,end=1,labels=labl)+
  theme(legend.text = element_text(size = 20),legend.title = element_text(size=15),axis.title=element_text(size=20),axis.text=element_text(size=20))+
  geom_hline(data = hline,aes(yintercept=yint,linetype = lt),color = "Black",show_guide = TRUE) + 
  scale_colour_discrete(guide = "none") + 
  scale_linetype_manual(name = 'Verdadeiro',values = 1,guide = "legend",labels=labl2)
## Salvando os graficos em eps e pdf 
dev.copy2pdf(file="Grafico3(Boxplot).pdf")
dev.copy2eps(file="Grafico3(Boxplot).eps")

### Resultados para Gama=1
dados4_trans=t(dados4)
dados4_rep=cbind(dados4_trans,rep(c(10,50,100,250,500),each=NREP))
dados4_rep=as.data.frame(dados4_rep)
names(dados4_rep)=c("r1","r2","r3","n")

dat4.m <- melt(dados4_rep,id.vars='n', measure.vars=c('r1','r2','r3'))
dat4.m$n=factor(dat4.m$n)
labl <- list(expression(hat(rho)[gamma]),expression(tilde(rho)[gamma]),expression(bar(rho)[S][gamma]))
hline <- data.frame(yint = Rho_Gama[4,1],lt = 'Valor Real')
labl2 <- list(expression(rho[gamma]))

ggplot(dat4.m) +  theme_bw() + xlab("Tamanho da amostra")+ ylab("Valores")+
  geom_boxplot(aes(x=n, y=value, fill=variable),position = position_dodge(width = .9))+
  guides(fill = guide_legend(title="Coeficientes"))+scale_fill_grey(start=0.5,end=1,labels=labl)+
  theme(legend.text = element_text(size = 20),legend.title = element_text(size=15),axis.title=element_text(size=20),axis.text=element_text(size=20))+
  geom_hline(data = hline,aes(yintercept=yint,linetype = lt),color = "Black",show_guide = TRUE) + 
  scale_colour_discrete(guide = "none") + 
  scale_linetype_manual(name = 'Verdadeiro',values = 1,guide = "legend",labels=labl2)
## Salvando os graficos em eps e pdf 
dev.copy2pdf(file="Grafico4(Boxplot).pdf")
dev.copy2eps(file="Grafico4(Boxplot).eps")
