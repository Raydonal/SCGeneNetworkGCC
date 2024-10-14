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

## Matrizes de estimativas
EMV=matrix(0,1,NREP) ## Estimativas de maxima verossimilhança
media_rho_mc=matrix(0,4*5,5) ## matriz de estimativas
var_rho_mc=matrix(0,4*5,5) ## matriz de variancias
vies_rho_mc=matrix(0,4*5,5) ## matriz do Vies do coeficiente de correlacao
eqm_rho_mc=matrix(0,4*5,5) ## matriz do Erro quadratico medio
raiz_eqm_rho_mc=matrix(0,4*5,5) ## Raiz quadrada do erro quadrático médio
ptm=proc.time() ## Tempo

# ## Diretorio da função de chinchilli
# setwd("../../Simulacao/ ~/Dropbox/Dissertacao-Cleber/Simulacao/R/NovaSimulacaoPerturbada/")

## Matrizes das estimativas de Chinchilli
ECH=matrix(0,5,NREP) ## Estimativas de chinchilli
Rho_Gama=matrix(0,5,5)
Gamma=c(0,0.1,0.5,0.9,1) ## Valores de Gamma
media_rho_ch=matrix(0,4*5,5)
var_rho_ch=matrix(0,4*5,5)
vies_rho_ch=matrix(0,4*5,5)
eqm_rho_ch=matrix(0,4*5,5)
raiz_eqm_rho_ch=matrix(0,4*5,5)
dyn.load("somac.so")

## Matrizes de estimativas de Spearman
ESP=matrix(0,1,NREP) ## Estimativas de maxima verossimilhança
media_rho_sp=matrix(0,4*5,5) ## matriz de estimativas
var_rho_sp=matrix(0,4*5,5) ## matriz de variancias
vies_rho_sp=matrix(0,4*5,5) ## matriz do Vies do coeficiente de correlacao
eqm_rho_sp=matrix(0,4*5,5) ## matriz do Erro quadratico medio
raiz_eqm_rho_sp=matrix(0,4*5,5) ## Raiz quadrada do erro quadrático médio

## Matriz da distribuicao de Chinchilli
EMVc=matrix(0,3,NREP)

## Porcentagem de Perturbacao
Porcentagem=c(10,30,50,70)

   ## Funcao de log-verossimilhanca do coeficiente de correlacao
 log_vero=function(x){
  rho_hat=x[1]
  c=0
  (-1)*(c-NOBS*(log(sqrt(sigma1))+log(sqrt(sigma2))+1/2*log(1-rho_hat^2)+mean[,1]^2/(2*(1-rho_hat^2)*sigma1)
                +mean[,2]^2/(2*(1-rho_hat^2)*sigma2)-(rho_hat*mean[,1]*mean[,2])/((1-rho_hat^2)*sigma1*sigma2))
        -sum(amostra_biv_per[,1]^2)/(2*(1-rho_hat^2)*sigma1)-sum(amostra_biv_per[,2]^2)/(2*(1-rho_hat^2)*sigma2)
        +(mean[,1]*sigma2-mean[,2]*sigma1*rho_hat)/(sigma1*sqrt(sigma2)*(1-rho_hat^2))*sum(amostra_biv_per[,1])
        +(mean[,2]*sigma1-mean[,1]*sigma2*rho_hat)/(sigma2*sqrt(sigma1)*(1-rho_hat^2))*sum(amostra_biv_per[,2])
        +rho_hat/(sqrt(sigma1*sigma2)*(1-rho_hat^2))*sum(amostra_biv_per[,1]*amostra_biv_per[,2]))
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
 rho=rho_vetor[m]
 cov=rho*sqrt(sigma1*sigma2) # definido no artigo 
 matrix_cov=matrix(c(sigma1,cov,
                     cov,sigma2),ncol=2) # matriz de variancia e covariancia
 
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
  
       cor_chinchilli=chinc(amostra_biv_per[,1],amostra_biv_per[,2],NOBS,gama)
       EMVc[2,i]=cor_chinchilli
       
       ESP[1,i]= cor(amostra_biv_per[,1],amostra_biv_per[,2],method=c("spearman"))
       EMVc[3,i]=hat_rho_gama(gama,ESP[1,i])
                  }
     
  if(j==1){
   #### Media das estimativas
   media_rho_mc[n,k]=mean(EMVc[1,]) # estimativas de maxima verossimilhança
   media_rho_ch[n,k]=mean(EMVc[2,]) # estimativas de chinchilli
   media_rho_sp[n,k]=mean(EMVc[3,]) # estimativas de Spearman
   #### Variancia das estimativas 
   var_rho_mc[n,k]=var(EMVc[1,])
   var_rho_ch[n,k]=var(EMVc[2,])
   var_rho_sp[n,k]=var(EMVc[3,])
   #### Vies do coeficiente de correlacao
   vies_rho_mc[n,k]=mean(EMVc[1, ])-rho_gama
   vies_rho_ch[n,k]=mean(EMVc[2, ])-rho_gama
   vies_rho_sp[n,k]=mean(EMVc[3, ])-rho_gama
   #### Erro quadratico medio
   eqm_rho_mc[n,k]=vies_rho_mc[n,k]^2 + var_rho_mc[n,k]
   eqm_rho_ch[n,k]=vies_rho_ch[n,k]^2 + var_rho_ch[n,k]
   eqm_rho_sp[n,k]=vies_rho_sp[n,k]^2 + var_rho_sp[n,k]
   #### Raiz quadrada do erro quadrático médio
   raiz_eqm_rho_mc[n,k]=round(sqrt(eqm_rho_mc[n,k]),4)
   raiz_eqm_rho_ch[n,k]=round(sqrt(eqm_rho_ch[n,k]),4)
   raiz_eqm_rho_sp[n,k]=round(sqrt(eqm_rho_sp[n,k]),4)}
    
    if(j==2){
      #### Media das estimativas
      media_rho_mc[4+n,k]=mean(EMVc[1,]) # estimativas de maxima verossimilhança
      media_rho_ch[4+n,k]=mean(EMVc[2,]) # estimativas de chinchilli
      media_rho_sp[4+n,k]=mean(EMVc[3,]) # estimativas de Spearman
      #### Variancia das estimativas 
      var_rho_mc[4+n,k]=var(EMVc[1,])
      var_rho_ch[4+n,k]=var(EMVc[2,])
      var_rho_sp[4+n,k]=var(EMVc[3,])
      #### Vies do coeficiente de correlacao
      vies_rho_mc[4+n,k]=mean(EMVc[1, ])-rho_gama
      vies_rho_ch[4+n,k]=mean(EMVc[2, ])-rho_gama
      vies_rho_sp[4+n,k]=mean(EMVc[3, ])-rho_gama
      #### Erro quadratico medio
      eqm_rho_mc[4+n,k]=vies_rho_mc[4+n,k]^2 + var_rho_mc[4+n,k]
      eqm_rho_ch[4+n,k]=vies_rho_ch[4+n,k]^2 + var_rho_ch[4+n,k]
      eqm_rho_sp[4+n,k]=vies_rho_sp[4+n,k]^2 + var_rho_sp[4+n,k]
      #### Raiz quadrada do erro quadrático médio
      raiz_eqm_rho_mc[4+n,k]=round(sqrt(eqm_rho_mc[4+n,k]),4)
      raiz_eqm_rho_ch[4+n,k]=round(sqrt(eqm_rho_ch[4+n,k]),4)
      raiz_eqm_rho_sp[4+n,k]=round(sqrt(eqm_rho_sp[4+n,k]),4)}
    
    if(j==3){
      #### Media das estimativas
      media_rho_mc[8+n,k]=mean(EMVc[1,]) # estimativas de maxima verossimilhança
      media_rho_ch[8+n,k]=mean(EMVc[2,]) # estimativas de chinchilli
      media_rho_sp[8+n,k]=mean(EMVc[3,]) # estimativas de Spearman
      #### Variancia das estimativas 
      var_rho_mc[8+n,k]=var(EMVc[1,])
      var_rho_ch[8+n,k]=var(EMVc[2,])
      var_rho_sp[8+n,k]=var(EMVc[3,])
      #### Vies do coeficiente de correlacao
      vies_rho_mc[8+n,k]=mean(EMVc[1, ])-rho_gama
      vies_rho_ch[8+n,k]=mean(EMVc[2, ])-rho_gama
      vies_rho_sp[8+n,k]=mean(EMVc[3, ])-rho_gama
      #### Erro quadratico medio
      eqm_rho_mc[8+n,k]=vies_rho_mc[8+n,k]^2 + var_rho_mc[8+n,k]
      eqm_rho_ch[8+n,k]=vies_rho_ch[8+n,k]^2 + var_rho_ch[8+n,k]
      eqm_rho_sp[8+n,k]=vies_rho_sp[8+n,k]^2 + var_rho_sp[8+n,k]
      #### Raiz quadrada do erro quadrático médio
      raiz_eqm_rho_mc[8+n,k]=round(sqrt(eqm_rho_mc[8+n,k]),4)
      raiz_eqm_rho_ch[8+n,k]=round(sqrt(eqm_rho_ch[8+n,k]),4)
      raiz_eqm_rho_sp[8+n,k]=round(sqrt(eqm_rho_sp[8+n,k]),4)}
    
    if(j==4){
      #### Media das estimativas
      media_rho_mc[12+n,k]=mean(EMVc[1,]) # estimativas de maxima verossimilhança
      media_rho_ch[12+n,k]=mean(EMVc[2,]) # estimativas de chinchilli
      media_rho_sp[12+n,k]=mean(EMVc[3,]) # estimativas de Spearman
      #### Variancia das estimativas 
      var_rho_mc[12+n,k]=var(EMVc[1,])
      var_rho_ch[12+n,k]=var(EMVc[2,])
      var_rho_sp[12+n,k]=var(EMVc[3,])
      #### Vies do coeficiente de correlacao
      vies_rho_mc[12+n,k]=mean(EMVc[1, ])-rho_gama
      vies_rho_ch[12+n,k]=mean(EMVc[2, ])-rho_gama
      vies_rho_sp[12+n,k]=mean(EMVc[3, ])-rho_gama
      #### Erro quadratico medio
      eqm_rho_mc[12+n,k]=vies_rho_mc[12+n,k]^2 + var_rho_mc[12+n,k]
      eqm_rho_ch[12+n,k]=vies_rho_ch[12+n,k]^2 + var_rho_ch[12+n,k]
      eqm_rho_sp[12+n,k]=vies_rho_sp[12+n,k]^2 + var_rho_sp[12+n,k]
      #### Raiz quadrada do erro quadrático médio
      raiz_eqm_rho_mc[12+n,k]=round(sqrt(eqm_rho_mc[12+n,k]),4)
      raiz_eqm_rho_ch[12+n,k]=round(sqrt(eqm_rho_ch[12+n,k]),4)
      raiz_eqm_rho_sp[12+n,k]=round(sqrt(eqm_rho_sp[12+n,k]),4)}
    
    if(j==5){
      #### Media das estimativas
      media_rho_mc[16+n,k]=mean(EMVc[1,]) # estimativas de maxima verossimilhança
      media_rho_ch[16+n,k]=mean(EMVc[2,]) # estimativas de chinchilli
      media_rho_sp[16+n,k]=mean(EMVc[3,]) # estimativas de Spearman
      #### Variancia das estimativas 
      var_rho_mc[16+n,k]=var(EMVc[1,])
      var_rho_ch[16+n,k]=var(EMVc[2,])
      var_rho_sp[16+n,k]=var(EMVc[3,])
      #### Vies do coeficiente de correlacao
      vies_rho_mc[16+n,k]=mean(EMVc[1, ])-rho_gama
      vies_rho_ch[16+n,k]=mean(EMVc[2, ])-rho_gama
      vies_rho_sp[16+n,k]=mean(EMVc[3, ])-rho_gama
      #### Erro quadratico medio
      eqm_rho_mc[16+n,k]=vies_rho_mc[16+n,k]^2 + var_rho_mc[16+n,k]
      eqm_rho_ch[16+n,k]=vies_rho_ch[16+n,k]^2 + var_rho_ch[16+n,k]
      eqm_rho_sp[16+n,k]=vies_rho_sp[16+n,k]^2 + var_rho_sp[16+n,k]
      #### Raiz quadrada do erro quadrático médio
      raiz_eqm_rho_mc[16+n,k]=round(sqrt(eqm_rho_mc[16+n,k]),4)
      raiz_eqm_rho_ch[16+n,k]=round(sqrt(eqm_rho_ch[16+n,k]),4)
      raiz_eqm_rho_sp[16+n,k]=round(sqrt(eqm_rho_sp[16+n,k]),4)}
  
   }
  }
 }   
   
     
     #############  Impressao dos resultados  #################
     
     media_rho=round(rbind(media_rho_mc,media_rho_ch,media_rho_sp),4)
     media_rep=cbind(media_rho,rep(c(10,30,50,70),4))# colocando os níveis de perturbacao
     media_order=media_rep[order(media_rep[,6]),]# ordenando
     
     vies_rho=round(rbind(vies_rho_mc,vies_rho_ch,vies_rho_sp),4)
     vies_rep=cbind(vies_rho,rep(c(10,30,50,70),4))
     vies_order=vies_rep[order(vies_rep[,6]),]
     
     var_rho=round(rbind(var_rho_mc,var_rho_ch,var_rho_sp),4)
     var_rep=cbind(var_rho,rep(c(10,30,50,70),4))
     var_order=var_rep[order(var_rep[,6]),]
 
 
     raiz_eqm_rho=round(rbind(raiz_eqm_rho_mc,raiz_eqm_rho_ch,raiz_eqm_rho_sp),4)
     raiz_eqm_rep=cbind(raiz_eqm_rho,rep(c(10,30,50,70),4))
     raiz_eqm_order=raiz_eqm_rep[order(raiz_eqm_rep[,6]),]   
  
     
 if (m==1){
   ## Saindo tabela
   write.table(media_rho,file="EstimativaPerturbacaoRho=01.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
   write.table(vies_rho,file="ViésPerturbacaoRho=01.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
   write.table(var_rho,file="VarianciaPerturbacaoRho=01.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
   write.table(raiz_eqm_rho,file="ReqmPerturbacaoRho=01.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
   write.table(media_order,file="TabelaEstimativaPerturbacaoRho=01.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500","Pertur"))
   write.table(vies_order,file="TabelaViésPerturbacaoRho=01.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500","Pertur"))
   write.table(var_order,file="TabelaVarianciaPerturbacaoRho=01.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500","Pertur"))
   write.table(raiz_eqm_order,file="TabelaReqmPerturbacaoRho=01.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500","Pertur"))
 }
 
 if (m==2){
   ## Saindo tabela
   write.table(media_rho,file="EstimativaPerturbacaoRho=03.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
   write.table(vies_rho,file="ViésPerturbacaoRho=03.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
   write.table(var_rho,file="VarianciaPerturbacaoRho=03.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
   write.table(raiz_eqm_rho,file="ReqmPerturbacaoRho=03.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
   write.table(media_order,file="TabelaEstimativaPerturbacaoRho=03.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500","Pertur"))
   write.table(vies_order,file="TabelaViésPerturbacaoRho=03.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500","Pertur"))
   write.table(var_order,file="TabelaVarianciaPerturbacaoRho=03.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500","Pertur"))
   write.table(raiz_eqm_order,file="TabelaReqmPerturbacaoRho=03.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500","Pertur"))
 }
 
 if (m==3){
   ## Saindo tabela
   write.table(media_rho,file="EstimativaPerturbacaoRho=05.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
   write.table(vies_rho,file="ViésPerturbacaoRho=05.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
   write.table(var_rho,file="VarianciaPerturbacaoRho=05.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
   write.table(raiz_eqm_rho,file="ReqmPerturbacaoRho=05.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
   write.table(media_order,file="TabelaEstimativaPerturbacaoRho=05.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500","Pertur"))
   write.table(vies_order,file="TabelaViésPerturbacaoRho=05.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500","Pertur"))
   write.table(var_order,file="TabelaVarianciaPerturbacaoRho=05.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500","Pertur"))
   write.table(raiz_eqm_order,file="TabelaReqmPerturbacaoRho=05.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500","Pertur"))
 }
 
 if (m==4){
   ## Saindo tabela
   write.table(media_rho,file="EstimativaPerturbacaoRho=07.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
   write.table(vies_rho,file="ViésPerturbacaoRho=07.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
   write.table(var_rho,file="VarianciaPerturbacaoRho=07.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
   write.table(raiz_eqm_rho,file="ReqmPerturbacaoRho=07.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
   write.table(media_order,file="TabelaEstimativaPerturbacaoRho=07.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500","Pertur"))
   write.table(vies_order,file="TabelaViésPerturbacaoRho=07.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500","Pertur"))
   write.table(var_order,file="TabelaVarianciaPerturbacaoRho=07.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500","Pertur"))
   write.table(raiz_eqm_order,file="TabelaReqmPerturbacaoRho=07.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500","Pertur"))
 }
 
 if (m==5){
   ## Saindo tabela
   write.table(media_rho,file="EstimativaPerturbacaoRho=09.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
   write.table(vies_rho,file="ViésPerturbacaoRho=09.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
   write.table(var_rho,file="VarianciaPerturbacaoRho=09.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
   write.table(raiz_eqm_rho,file="ReqmPerturbacaoRho=09.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500"))
   write.table(media_order,file="TabelaEstimativaPerturbacaoRho=09.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500","Pertur"))
   write.table(vies_order,file="TabelaViésPerturbacaoRho=09.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500","Pertur"))
   write.table(var_order,file="TabelaVarianciaPerturbacaoRho=09.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500","Pertur"))
   write.table(raiz_eqm_order,file="TabelaReqmPerturbacaoRho=09.csv"
               ,row.names=FALSE,sep=";",dec=".",col.names=c("10","50","100","250","500","Pertur"))}
}
## Tempo
proc.time()-ptm


