#################################################################################################
######          Simulacao do coeficiente de correlacao da normal bivariada                 ######
######            Comparando com o coeficiente de correlacao de Chinchilli                 ######
######                                                                                     ######
#################################################################################################

  #### Escolhendo Gerador #####
rm(list=ls(all=TRUE)) 
RNGkind("Mersenne-Twister")

library(mvtnorm) ## gerar numeros normais e t
library(BMS) ## funcao hypergeometric

  ###### Variaveis
NOBS_vetor=c(10,50,100,250,500) ## Tamanho das amostras 10,50,100,250,500
NREP=5000
mean=matrix(c(0,0),ncol=2)
sigma1=1
sigma2=1
rho_vetor=c(0.1,0.3,0.5,0.7,0.9)
Porcentagem=c(10,30,50,70)

## Matrizes de estimativas
EMV=matrix(0,1,NREP) ## Estimativas de maxima verossimilhança
ptm=proc.time() ## Tempo

## Diretorio da função de chinchilli
# setwd("~/Dropbox/Dissertacao-Cleber/Simulacao/R/NovaSimluacao/")

## Matrizes das estimativas de Chinchilli
Rho_Gama=matrix(0,5,1)
Gamma=c(0,0.1,0.5,0.9,1) ## Valores de Gamma
dyn.load("somac.so")

## Matriz de estimativas de Spearman
ESP=matrix(0,1,NREP) ## Estimativas de maxima verossimilhança

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

###### Funcao para porcentagem de perturbacao
pertur=function(porcentagem){
  valor=porcentagem/100
}

   ##### Inicio do Monte Carlo ######
 for(m in 1:5){
  rho=rho_vetor[m]  ## mudar valores de rho retirando do rho_vetor
  cov=rho*sqrt(sigma1*sigma2) # definido no artigo 
  matrix_cov=matrix(c(sigma1,cov,
                      cov ,sigma2),ncol=2) # matriz de variancia e covariancia
  for (j in 1:5){
   gama=Gamma[j]
   rho_gama=hat_rho_gama(gama,rho)
  for(k in 1:5){
   NOBS=NOBS_vetor[k]
    
	 set.seed(2) ## Semente

   for(n in 1:4){
     porcentagem=Porcentagem[n] # niveis de porcentagem
     por_perturbacao=pertur(porcentagem) #porcentagem
     
     for(i in 1:NREP){
       a=runif(NOBS)
       amostra_per=rmvt(NOBS,sigma=matrix_cov,df=3) ## gerando numeros t
       amostra_biv=rmvnorm(NOBS,mean,matrix_cov) ## gerando numeros normais
       amostra_biv_per=matrix(c(ifelse(a>por_perturbacao,amostra_biv[,1],amostra_per[,1])
                                ,ifelse(a>por_perturbacao,amostra_biv[,2],amostra_per[,2])),ncol=2) ## perturbacao
       
       ##### Estimacao sem o uso do gradiente
       ir=optim(0.1,log_vero,NULL,method="BFGS")
       EMV[1,i]=ir$par
       EMVc[1,i]=hat_rho_gama(gama,EMV[1,i]) # estimativas de maxima verossimilhanca na distribuicao
  
       cor_chinchilli=chinc(amostra_biv[,1],amostra_biv[,2],NOBS,gama)
       EMVc[2,i]=cor_chinchilli
       
       ESP[1,i]= cor(amostra_biv[,1],amostra_biv[,2],method=c("spearman"))
       EMVc[3,i]=hat_rho_gama(gama,ESP[1,i])
                  }
   
   ####matriz com os valores de rho_gamma
   Rho_Gama[j,1]=rho_gama
      
   if (k==1){ 
     if (n==1){
       ENREP1010=round(EMVc,4) ## Estimativas para o número de replicas com  n=10 e 10 % perturbacao
          }
     if (n==2){
       ENREP1030=round(EMVc,4) ## Estimativas para o número de replicas com  n=10 e 30 % perturbacao
           }
     if (n==3){
        ENREP1050=round(EMVc,4) ## Estimativas para o número de replicas com  n=10 e 50 % perturbacao
     }
     if (n==4){
        ENREP1070=round(EMVc,4) ## Estimativas para o número de replicas com  n=10 e 70 % perturbacao
       }
   }         
     if (k==2){
       if (n==1){
          ENREP5010=round(EMVc,4) ## Estimativas para o número de replicas com  n=50 e 10 % perturbacao
         }
       if (n==2){
          ENREP5030=round(EMVc,4) ## Estimativas para o número de replicas com  n=50 e 30 % perturbacao
         }
       if (n==3){
         ENREP5050=round(EMVc,4) ## Estimativas para o número de replicas com  n=50 e 50 % perturbacao
       }
       if (n==4){
         ENREP5070=round(EMVc,4) ## Estimativas para o número de replicas com  n=50 e 70 % perturbacao
       } 
     }
      if (k==3){
        if (n==1){
          ENREP10010=round(EMVc,4) ## Estimativas para o número de replicas com  n=100 e 10 % perturbacao
        }
        if (n==2){
          ENREP10030=round(EMVc,4) ## Estimativas para o número de replicas com  n=100 e 30 % perturbacao
        }
        if (n==3){
          ENREP10050=round(EMVc,4) ## Estimativas para o número de replicas com  n=100 e 50 % perturbacao
        }
        if (n==4){
          ENREP10070=round(EMVc,4) ## Estimativas para o número de replicas com  n=100 e 70 % perturbacao
        }
      }
       if (k==4){
         if (n==1){
           ENREP25010=round(EMVc,4) ## Estimativas para o número de replicas com  n=250 e 10 % perturbacao
         }
         if (n==2){
           ENREP25030=round(EMVc,4) ## Estimativas para o número de replicas com  n=250 e 30 % perturbacao
         }
         if (n==3){
           ENREP25050=round(EMVc,4) ## Estimativas para o número de replicas com  n=250 e 50 % perturbacao
         }
         if (n==4){
           ENREP25070=round(EMVc,4) ## Estimativas para o número de replicas com  n=250 e 70 % perturbacao
         }
       }
        if (k==5){
          if (n==1){
            ENREP50010=round(EMVc,4) ## Estimativas para o número de replicas com  n=500 e 10 % perturbacao
          }
          if (n==2){
            ENREP50030=round(EMVc,4) ## Estimativas para o número de replicas com  n=500 e 30 % perturbacao
          }
          if (n==3){
            ENREP50050=round(EMVc,4) ## Estimativas para o número de replicas com  n=500 e 50 % perturbacao
          }
          if (n==4){
            ENREP50070=round(EMVc,4) ## Estimativas para o número de replicas com  n=500 e 70 % perturbacao
          }
        }
     }
  }
     if (j==1){
       dados10gamma0=cbind(ENREP1010,ENREP5010,ENREP10010,ENREP25010,ENREP50010) ## Matrix com os varios n para gama=0 e 10% pertur
       dados30gamma0=cbind(ENREP1030,ENREP5030,ENREP10030,ENREP25030,ENREP50030) ## Matrix com os varios n para gama=0 e 10% pertur
       dados50gamma0=cbind(ENREP1050,ENREP5050,ENREP10050,ENREP25050,ENREP50050) ## Matrix com os varios n para gama=0 e 10% pertur
       dados70gamma0=cbind(ENREP1070,ENREP5070,ENREP10070,ENREP25070,ENREP50070) ## Matrix com os varios n para gama=0 e 10% pertur
         }
     if (j==2){
       dados10gamma01=cbind(ENREP1010,ENREP5010,ENREP10010,ENREP25010,ENREP50010) ## Matrix com os varios n para gama=01 e 10% pertur
       dados30gamma01=cbind(ENREP1030,ENREP5030,ENREP10030,ENREP25030,ENREP50030) ## Matrix com os varios n para gama=01 e 10% pertur
       dados50gamma01=cbind(ENREP1050,ENREP5050,ENREP10050,ENREP25050,ENREP50050) ## Matrix com os varios n para gama=01 e 10% pertur
       dados70gamma01=cbind(ENREP1070,ENREP5070,ENREP10070,ENREP25070,ENREP50070) ## Matrix com os varios n para gama=01 e 10% pertur
     }
      if (j==3){
        dados10gamma05=cbind(ENREP1010,ENREP5010,ENREP10010,ENREP25010,ENREP50010) ## Matrix com os varios n para gama=05 e 10% pertur
        dados30gamma05=cbind(ENREP1030,ENREP5030,ENREP10030,ENREP25030,ENREP50030) ## Matrix com os varios n para gama=05 e 10% pertur
        dados50gamma05=cbind(ENREP1050,ENREP5050,ENREP10050,ENREP25050,ENREP50050) ## Matrix com os varios n para gama=05 e 10% pertur
        dados70gamma05=cbind(ENREP1070,ENREP5070,ENREP10070,ENREP25070,ENREP50070) ## Matrix com os varios n para gama=05 e 10% pertur
      }
       if (j==4){
         dados10gamma09=cbind(ENREP1010,ENREP5010,ENREP10010,ENREP25010,ENREP50010) ## Matrix com os varios n para gama=09 e 10% pertur
         dados30gamma09=cbind(ENREP1030,ENREP5030,ENREP10030,ENREP25030,ENREP50030) ## Matrix com os varios n para gama=09 e 10% pertur
         dados50gamma09=cbind(ENREP1050,ENREP5050,ENREP10050,ENREP25050,ENREP50050) ## Matrix com os varios n para gama=09 e 10% pertur
         dados70gamma09=cbind(ENREP1070,ENREP5070,ENREP10070,ENREP25070,ENREP50070) ## Matrix com os varios n para gama=09 e 10% pertur
       }
        if (j==5){
          dados10gamma1=cbind(ENREP1010,ENREP5010,ENREP10010,ENREP25010,ENREP50010) ## Matrix com os varios n para gama=1 e 10% pertur
          dados30gamma1=cbind(ENREP1030,ENREP5030,ENREP10030,ENREP25030,ENREP50030) ## Matrix com os varios n para gama=1 e 10% pertur
          dados50gamma1=cbind(ENREP1050,ENREP5050,ENREP10050,ENREP25050,ENREP50050) ## Matrix com os varios n para gama=1 e 10% pertur
          dados70gamma1=cbind(ENREP1070,ENREP5070,ENREP10070,ENREP25070,ENREP50070) ## Matrix com os varios n para gama=1 e 10% pertur
        }
 
   }
proc.time()-ptm
     if (m==1){ 
# setwd("~/Dropbox/Dissertacao-Cleber/Simulacao/Graficos/colorido/NovaSimulacaoPerturbada/Boxplot/Rho01/")
save.image("Rho01.Rdata")

#### Gerando Gráficos

library(reshape2)
library(ggplot2)
##########################################################
## Resultados para Gamma=0
##########################################################

dados1_trans=t(dados10gamma0)
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
### Salvando os graficos em eps e pdf 
dev.copy2pdf(file="Grafico1(BoxplotGama0-10).pdf")
dev.copy2eps(file="Grafico1(BoxplotGama0-10).eps")

dados1_trans=t(dados30gamma0)
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
### Salvando os graficos em eps e pdf 
dev.copy2pdf(file="Grafico1(BoxplotGama0-30).pdf")
dev.copy2eps(file="Grafico1(BoxplotGama0-30).eps")

dados1_trans=t(dados50gamma0)
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
### Salvando os graficos em eps e pdf 
dev.copy2pdf(file="Grafico1(BoxplotGama0-50).pdf")
dev.copy2eps(file="Grafico1(BoxplotGama0-50).eps")

#################################################################
### Resultados para Gama=0.5
#################################################################

dados3_trans=t(dados10gamma05)
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
### Salvando os graficos em eps e pdf 
dev.copy2pdf(file="Grafico3(BoxplotGama05-10).pdf")
dev.copy2eps(file="Grafico3(BoxplotGama05-10).eps")

dados3_trans=t(dados30gamma05)
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
### Salvando os graficos em eps e pdf 
dev.copy2pdf(file="Grafico3(BoxplotGama05-30).pdf")
dev.copy2eps(file="Grafico3(BoxplotGama05-30).eps")

dados3_trans=t(dados50gamma05)
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
### Salvando os graficos em eps e pdf 
dev.copy2pdf(file="Grafico3(BoxplotGama05-50).pdf")
dev.copy2eps(file="Grafico3(BoxplotGama05-50).eps")

##########################################################################
### Resultados para Gama=1
###########################################################################
dados5_trans=t(dados10gamma1)
dados5_rep=cbind(dados5_trans,rep(c(10,50,100,250,500),each=NREP))
dados5_rep=as.data.frame(dados5_rep)
names(dados5_rep)=c("r1","r2","r3","n")

dat5.m <- melt(dados5_rep,id.vars='n', measure.vars=c('r1','r2','r3'))
dat5.m$n=factor(dat5.m$n)
labl <- list(expression(hat(rho)[gamma]),expression(tilde(rho)[gamma]),expression(bar(rho)[S][gamma]))
hline <- data.frame(yint = Rho_Gama[5,1],lt = 'Valor Real')
labl2 <- list(expression(rho[gamma]))

ggplot(dat5.m) +  theme_bw() + xlab("Tamanho da amostra")+ ylab("Valores")+
  geom_boxplot(aes(x=n, y=value, fill=variable),position = position_dodge(width = .9))+
  guides(fill = guide_legend(title="Coeficientes"))+scale_fill_grey(start=0.5,end=1,labels=labl)+
  theme(legend.text = element_text(size = 20),legend.title = element_text(size=15),axis.title=element_text(size=20),axis.text=element_text(size=20))+
  geom_hline(data = hline,aes(yintercept=yint,linetype = lt),color = "Black",show_guide = TRUE) + 
  scale_colour_discrete(guide = "none") + 
  scale_linetype_manual(name = 'Verdadeiro',values = 1,guide = "legend",labels=labl2)
### Salvando os graficos em eps e pdf 
dev.copy2pdf(file="Grafico5(BoxplotGama1-10).pdf")
dev.copy2eps(file="Grafico5(BoxplotGama1-10).eps")

dados5_trans=t(dados30gamma1)
dados5_rep=cbind(dados5_trans,rep(c(10,50,100,250,500),each=NREP))
dados5_rep=as.data.frame(dados5_rep)
names(dados5_rep)=c("r1","r2","r3","n")

dat5.m <- melt(dados5_rep,id.vars='n', measure.vars=c('r1','r2','r3'))
dat5.m$n=factor(dat5.m$n)
labl <- list(expression(hat(rho)[gamma]),expression(tilde(rho)[gamma]),expression(bar(rho)[S][gamma]))
hline <- data.frame(yint = Rho_Gama[5,1],lt = 'Valor Real')
labl2 <- list(expression(rho[gamma]))

ggplot(dat5.m) +  theme_bw() + xlab("Tamanho da amostra")+ ylab("Valores")+
  geom_boxplot(aes(x=n, y=value, fill=variable),position = position_dodge(width = .9))+
  guides(fill = guide_legend(title="Coeficientes"))+scale_fill_grey(start=0.5,end=1,labels=labl)+
  theme(legend.text = element_text(size = 20),legend.title = element_text(size=15),axis.title=element_text(size=20),axis.text=element_text(size=20))+
  geom_hline(data = hline,aes(yintercept=yint,linetype = lt),color = "Black",show_guide = TRUE) + 
  scale_colour_discrete(guide = "none") + 
  scale_linetype_manual(name = 'Verdadeiro',values = 1,guide = "legend",labels=labl2)
### Salvando os graficos em eps e pdf 
dev.copy2pdf(file="Grafico5(BoxplotGama1-30).pdf")
dev.copy2eps(file="Grafico5(BoxplotGama1-30).eps")

dados5_trans=t(dados50gamma1)
dados5_rep=cbind(dados5_trans,rep(c(10,50,100,250,500),each=NREP))
dados5_rep=as.data.frame(dados5_rep)
names(dados5_rep)=c("r1","r2","r3","n")

dat5.m <- melt(dados5_rep,id.vars='n', measure.vars=c('r1','r2','r3'))
dat5.m$n=factor(dat5.m$n)
labl <- list(expression(hat(rho)[gamma]),expression(tilde(rho)[gamma]),expression(bar(rho)[S][gamma]))
hline <- data.frame(yint = Rho_Gama[5,1],lt = 'Valor Real')
labl2 <- list(expression(rho[gamma]))

ggplot(dat5.m) +  theme_bw() + xlab("Tamanho da amostra")+ ylab("Valores")+
  geom_boxplot(aes(x=n, y=value, fill=variable),position = position_dodge(width = .9))+
  guides(fill = guide_legend(title="Coeficientes"))+scale_fill_grey(start=0.5,end=1,labels=labl)+
  theme(legend.text = element_text(size = 20),legend.title = element_text(size=15),axis.title=element_text(size=20),axis.text=element_text(size=20))+
  geom_hline(data = hline,aes(yintercept=yint,linetype = lt),color = "Black",show_guide = TRUE) + 
  scale_colour_discrete(guide = "none") + 
  scale_linetype_manual(name = 'Verdadeiro',values = 1,guide = "legend",labels=labl2)
### Salvando os graficos em eps e pdf 
dev.copy2pdf(file="Grafico5(BoxplotGama1-50).pdf")
dev.copy2eps(file="Grafico5(BoxplotGama1-50).eps")
     }
  
  if (m==3){ 
    setwd("~/Dropbox/Dissertacao-Cleber/Simulacao/Graficos/colorido/NovaSimulacaoPerturbada/Boxplot/Rho05/")
    save.image("Rho05.Rdata")
    
    #### Gerando Gráficos
    
    library(reshape2)
    library(ggplot2)
    ##########################################################
    ## Resultados para Gamma=0
    ##########################################################
    
    dados1_trans=t(dados10gamma0)
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
    ### Salvando os graficos em eps e pdf 
    dev.copy2pdf(file="Grafico1(BoxplotGama0-10).pdf")
    dev.copy2eps(file="Grafico1(BoxplotGama0-10).eps")
    
    dados1_trans=t(dados30gamma0)
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
    ### Salvando os graficos em eps e pdf 
    dev.copy2pdf(file="Grafico1(BoxplotGama0-30).pdf")
    dev.copy2eps(file="Grafico1(BoxplotGama0-30).eps")
    
    dados1_trans=t(dados50gamma0)
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
    ### Salvando os graficos em eps e pdf 
    dev.copy2pdf(file="Grafico1(BoxplotGama0-50).pdf")
    dev.copy2eps(file="Grafico1(BoxplotGama0-50).eps")
    
    #################################################################
    ### Resultados para Gama=0.5
    #################################################################
    
    dados3_trans=t(dados10gamma05)
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
    ### Salvando os graficos em eps e pdf 
    dev.copy2pdf(file="Grafico3(BoxplotGama05-10).pdf")
    dev.copy2eps(file="Grafico3(BoxplotGama05-10).eps")
    
    dados3_trans=t(dados30gamma05)
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
    ### Salvando os graficos em eps e pdf 
    dev.copy2pdf(file="Grafico3(BoxplotGama05-30).pdf")
    dev.copy2eps(file="Grafico3(BoxplotGama05-30).eps")
    
    dados3_trans=t(dados50gamma05)
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
    ### Salvando os graficos em eps e pdf 
    dev.copy2pdf(file="Grafico3(BoxplotGama05-50).pdf")
    dev.copy2eps(file="Grafico3(BoxplotGama05-50).eps")
    
    ##########################################################################
    ### Resultados para Gama=1
    ###########################################################################
    dados5_trans=t(dados10gamma1)
    dados5_rep=cbind(dados5_trans,rep(c(10,50,100,250,500),each=NREP))
    dados5_rep=as.data.frame(dados5_rep)
    names(dados5_rep)=c("r1","r2","r3","n")
    
    dat5.m <- melt(dados5_rep,id.vars='n', measure.vars=c('r1','r2','r3'))
    dat5.m$n=factor(dat5.m$n)
    labl <- list(expression(hat(rho)[gamma]),expression(tilde(rho)[gamma]),expression(bar(rho)[S][gamma]))
    hline <- data.frame(yint = Rho_Gama[5,1],lt = 'Valor Real')
    labl2 <- list(expression(rho[gamma]))
    
    ggplot(dat5.m) +  theme_bw() + xlab("Tamanho da amostra")+ ylab("Valores")+
      geom_boxplot(aes(x=n, y=value, fill=variable),position = position_dodge(width = .9))+
      guides(fill = guide_legend(title="Coeficientes"))+scale_fill_grey(start=0.5,end=1,labels=labl)+
      theme(legend.text = element_text(size = 20),legend.title = element_text(size=15),axis.title=element_text(size=20),axis.text=element_text(size=20))+
      geom_hline(data = hline,aes(yintercept=yint,linetype = lt),color = "Black",show_guide = TRUE) + 
      scale_colour_discrete(guide = "none") + 
      scale_linetype_manual(name = 'Verdadeiro',values = 1,guide = "legend",labels=labl2)
    ### Salvando os graficos em eps e pdf 
    dev.copy2pdf(file="Grafico5(BoxplotGama1-10).pdf")
    dev.copy2eps(file="Grafico5(BoxplotGama1-10).eps")
    
    dados5_trans=t(dados30gamma1)
    dados5_rep=cbind(dados5_trans,rep(c(10,50,100,250,500),each=NREP))
    dados5_rep=as.data.frame(dados5_rep)
    names(dados5_rep)=c("r1","r2","r3","n")
    
    dat5.m <- melt(dados5_rep,id.vars='n', measure.vars=c('r1','r2','r3'))
    dat5.m$n=factor(dat5.m$n)
    labl <- list(expression(hat(rho)[gamma]),expression(tilde(rho)[gamma]),expression(bar(rho)[S][gamma]))
    hline <- data.frame(yint = Rho_Gama[5,1],lt = 'Valor Real')
    labl2 <- list(expression(rho[gamma]))
    
    ggplot(dat5.m) +  theme_bw() + xlab("Tamanho da amostra")+ ylab("Valores")+
      geom_boxplot(aes(x=n, y=value, fill=variable),position = position_dodge(width = .9))+
      guides(fill = guide_legend(title="Coeficientes"))+scale_fill_grey(start=0.5,end=1,labels=labl)+
      theme(legend.text = element_text(size = 20),legend.title = element_text(size=15),axis.title=element_text(size=20),axis.text=element_text(size=20))+
      geom_hline(data = hline,aes(yintercept=yint,linetype = lt),color = "Black",show_guide = TRUE) + 
      scale_colour_discrete(guide = "none") + 
      scale_linetype_manual(name = 'Verdadeiro',values = 1,guide = "legend",labels=labl2)
    ### Salvando os graficos em eps e pdf 
    dev.copy2pdf(file="Grafico5(BoxplotGama1-30).pdf")
    dev.copy2eps(file="Grafico5(BoxplotGama1-30).eps")
    
    dados5_trans=t(dados50gamma1)
    dados5_rep=cbind(dados5_trans,rep(c(10,50,100,250,500),each=NREP))
    dados5_rep=as.data.frame(dados5_rep)
    names(dados5_rep)=c("r1","r2","r3","n")
    
    dat5.m <- melt(dados5_rep,id.vars='n', measure.vars=c('r1','r2','r3'))
    dat5.m$n=factor(dat5.m$n)
    labl <- list(expression(hat(rho)[gamma]),expression(tilde(rho)[gamma]),expression(bar(rho)[S][gamma]))
    hline <- data.frame(yint = Rho_Gama[5,1],lt = 'Valor Real')
    labl2 <- list(expression(rho[gamma]))
    
    ggplot(dat5.m) +  theme_bw() + xlab("Tamanho da amostra")+ ylab("Valores")+
      geom_boxplot(aes(x=n, y=value, fill=variable),position = position_dodge(width = .9))+
      guides(fill = guide_legend(title="Coeficientes"))+scale_fill_grey(start=0.5,end=1,labels=labl)+
      theme(legend.text = element_text(size = 20),legend.title = element_text(size=15),axis.title=element_text(size=20),axis.text=element_text(size=20))+
      geom_hline(data = hline,aes(yintercept=yint,linetype = lt),color = "Black",show_guide = TRUE) + 
      scale_colour_discrete(guide = "none") + 
      scale_linetype_manual(name = 'Verdadeiro',values = 1,guide = "legend",labels=labl2)
    ### Salvando os graficos em eps e pdf 
    dev.copy2pdf(file="Grafico5(BoxplotGama1-50).pdf")
    dev.copy2eps(file="Grafico5(BoxplotGama1-50).eps")
  }
  
  if (m==5){ 
    setwd("~/Dropbox/Dissertacao-Cleber/Simulacao/Graficos/colorido/NovaSimulacaoPerturbada/Boxplot/Rho09/")
    save.image("Rho09.Rdata")
    
    #### Gerando Gráficos
    
    library(reshape2)
    library(ggplot2)
    ##########################################################
    ## Resultados para Gamma=0
    ##########################################################
    
    dados1_trans=t(dados10gamma0)
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
    ### Salvando os graficos em eps e pdf 
    dev.copy2pdf(file="Grafico1(BoxplotGama0-10).pdf")
    dev.copy2eps(file="Grafico1(BoxplotGama0-10).eps")
    
    dados1_trans=t(dados30gamma0)
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
    ### Salvando os graficos em eps e pdf 
    dev.copy2pdf(file="Grafico1(BoxplotGama0-30).pdf")
    dev.copy2eps(file="Grafico1(BoxplotGama0-30).eps")
    
    dados1_trans=t(dados50gamma0)
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
    ### Salvando os graficos em eps e pdf 
    dev.copy2pdf(file="Grafico1(BoxplotGama0-50).pdf")
    dev.copy2eps(file="Grafico1(BoxplotGama0-50).eps")
    
    #################################################################
    ### Resultados para Gama=0.5
    #################################################################
    
    dados3_trans=t(dados10gamma05)
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
    ### Salvando os graficos em eps e pdf 
    dev.copy2pdf(file="Grafico3(BoxplotGama05-10).pdf")
    dev.copy2eps(file="Grafico3(BoxplotGama05-10).eps")
    
    dados3_trans=t(dados30gamma05)
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
    ### Salvando os graficos em eps e pdf 
    dev.copy2pdf(file="Grafico3(BoxplotGama05-30).pdf")
    dev.copy2eps(file="Grafico3(BoxplotGama05-30).eps")
    
    dados3_trans=t(dados50gamma05)
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
    ### Salvando os graficos em eps e pdf 
    dev.copy2pdf(file="Grafico3(BoxplotGama05-50).pdf")
    dev.copy2eps(file="Grafico3(BoxplotGama05-50).eps")
    
    ##########################################################################
    ### Resultados para Gama=1
    ###########################################################################
    dados5_trans=t(dados10gamma1)
    dados5_rep=cbind(dados5_trans,rep(c(10,50,100,250,500),each=NREP))
    dados5_rep=as.data.frame(dados5_rep)
    names(dados5_rep)=c("r1","r2","r3","n")
    
    dat5.m <- melt(dados5_rep,id.vars='n', measure.vars=c('r1','r2','r3'))
    dat5.m$n=factor(dat5.m$n)
    labl <- list(expression(hat(rho)[gamma]),expression(tilde(rho)[gamma]),expression(bar(rho)[S][gamma]))
    hline <- data.frame(yint = Rho_Gama[5,1],lt = 'Valor Real')
    labl2 <- list(expression(rho[gamma]))
    
    ggplot(dat5.m) +  theme_bw() + xlab("Tamanho da amostra")+ ylab("Valores")+
      geom_boxplot(aes(x=n, y=value, fill=variable),position = position_dodge(width = .9))+
      guides(fill = guide_legend(title="Coeficientes"))+scale_fill_grey(start=0.5,end=1,labels=labl)+
      theme(legend.text = element_text(size = 20),legend.title = element_text(size=15),axis.title=element_text(size=20),axis.text=element_text(size=20))+
      geom_hline(data = hline,aes(yintercept=yint,linetype = lt),color = "Black",show_guide = TRUE) + 
      scale_colour_discrete(guide = "none") + 
      scale_linetype_manual(name = 'Verdadeiro',values = 1,guide = "legend",labels=labl2)
    ### Salvando os graficos em eps e pdf 
    dev.copy2pdf(file="Grafico5(BoxplotGama1-10).pdf")
    dev.copy2eps(file="Grafico5(BoxplotGama1-10).eps")
    
    dados5_trans=t(dados30gamma1)
    dados5_rep=cbind(dados5_trans,rep(c(10,50,100,250,500),each=NREP))
    dados5_rep=as.data.frame(dados5_rep)
    names(dados5_rep)=c("r1","r2","r3","n")
    
    dat5.m <- melt(dados5_rep,id.vars='n', measure.vars=c('r1','r2','r3'))
    dat5.m$n=factor(dat5.m$n)
    labl <- list(expression(hat(rho)[gamma]),expression(tilde(rho)[gamma]),expression(bar(rho)[S][gamma]))
    hline <- data.frame(yint = Rho_Gama[5,1],lt = 'Valor Real')
    labl2 <- list(expression(rho[gamma]))
    
    ggplot(dat5.m) +  theme_bw() + xlab("Tamanho da amostra")+ ylab("Valores")+
      geom_boxplot(aes(x=n, y=value, fill=variable),position = position_dodge(width = .9))+
      guides(fill = guide_legend(title="Coeficientes"))+scale_fill_grey(start=0.5,end=1,labels=labl)+
      theme(legend.text = element_text(size = 20),legend.title = element_text(size=15),axis.title=element_text(size=20),axis.text=element_text(size=20))+
      geom_hline(data = hline,aes(yintercept=yint,linetype = lt),color = "Black",show_guide = TRUE) + 
      scale_colour_discrete(guide = "none") + 
      scale_linetype_manual(name = 'Verdadeiro',values = 1,guide = "legend",labels=labl2)
    ### Salvando os graficos em eps e pdf 
    dev.copy2pdf(file="Grafico5(BoxplotGama1-30).pdf")
    dev.copy2eps(file="Grafico5(BoxplotGama1-30).eps")
    
    dados5_trans=t(dados50gamma1)
    dados5_rep=cbind(dados5_trans,rep(c(10,50,100,250,500),each=NREP))
    dados5_rep=as.data.frame(dados5_rep)
    names(dados5_rep)=c("r1","r2","r3","n")
    
    dat5.m <- melt(dados5_rep,id.vars='n', measure.vars=c('r1','r2','r3'))
    dat5.m$n=factor(dat5.m$n)
    labl <- list(expression(hat(rho)[gamma]),expression(tilde(rho)[gamma]),expression(bar(rho)[S][gamma]))
    hline <- data.frame(yint = Rho_Gama[5,1],lt = 'Valor Real')
    labl2 <- list(expression(rho[gamma]))
    
    ggplot(dat5.m) +  theme_bw() + xlab("Tamanho da amostra")+ ylab("Valores")+
      geom_boxplot(aes(x=n, y=value, fill=variable),position = position_dodge(width = .9))+
      guides(fill = guide_legend(title="Coeficientes"))+scale_fill_grey(start=0.5,end=1,labels=labl)+
      theme(legend.text = element_text(size = 20),legend.title = element_text(size=15),axis.title=element_text(size=20),axis.text=element_text(size=20))+
      geom_hline(data = hline,aes(yintercept=yint,linetype = lt),color = "Black",show_guide = TRUE) + 
      scale_colour_discrete(guide = "none") + 
      scale_linetype_manual(name = 'Verdadeiro',values = 1,guide = "legend",labels=labl2)
    ### Salvando os graficos em eps e pdf 
    dev.copy2pdf(file="Grafico5(BoxplotGama1-50).pdf")
    dev.copy2eps(file="Grafico5(BoxplotGama1-50).eps")
  }
}

