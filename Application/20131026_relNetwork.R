#############################################################################
##
## Script to run relevance networks 

## Loading maigesPack 
library(maigesPack)
source("corChinchilli.R")
source("bootstrapCor.R")
source("relNetworkB.R")



#############################################################################
## 
## Loading normalized and summarized from sáb 12 out 2024 19:42:16 -03

##
load("../sarahSummSingle_Ras_KEGG_20130330.RData")


#########################################################################
##
## Doing Relevance Networks for Ras small gene groupo
##
## Given the bad representation for big Ras, will do all comparisons
## made for DE genes, using the small Ras (genesetRas)
##

## Comparing selecting only stomach samples
idx <- sarahSummSingle@Slabels$Orgao == "Estomago"

relNet.pearson <- relNetworkB(sarahSummSingle[,idx],
                                   sLabelID="DS_DIAG_PROJ", geneGrp=NULL,
                                   path='Renin-angiotensin system',
                                   type="pearson")


  

relNet.spearman <- relNetworkB(sarahSummSingle[,idx],
                              sLabelID="DS_DIAG_PROJ", geneGrp=NULL,
                              path='Renin-angiotensin system',
                              type="spearman")


relNet.chi.gamma1.0 <- relNetworkB(sarahSummSingle[,idx],
                                   sLabelID="DS_DIAG_PROJ", geneGrp=NULL,
                                   path='Renin-angiotensin system',
                                   type="chinchilli", gamma=1)

relNet.chi.gamma0.8 <- relNetworkB(sarahSummSingle[,idx],
                                   sLabelID="DS_DIAG_PROJ", geneGrp=NULL,
                                   path='Renin-angiotensin system',
                                   type="chinchilli", gamma=.8)

relNet.chi.gamma0.5 <- relNetworkB(sarahSummSingle[,idx],
                                   sLabelID="DS_DIAG_PROJ", geneGrp=NULL,
                                   path='Renin-angiotensin system',
                                   type="chinchilli", gamma=.5)

relNet.chi.gamma0.3 <- relNetworkB(sarahSummSingle[,idx],
                                   sLabelID="DS_DIAG_PROJ", geneGrp=NULL,
                                   path='Renin-angiotensin system',
                                   type="chinchilli", gamma=.3)

relNet.chi.gamma0.0 <- relNetworkB(sarahSummSingle[,idx],
                                   sLabelID="DS_DIAG_PROJ", geneGrp=NULL,
                                   path='Renin-angiotensin system',
                                   type="chinchilli", gamma=0)



png("figuras/teste_estomago_normal_new.png", width=1200, height=800)
par(mfrow=c(3,3))

plot(relNet.pearson , cutPval=NULL, cutCor=.5, main="pearson", ylab="pearson")
plot(relNet.spearman, cutPval=NULL, cutCor=.5, main="spearman", ylab="spearman")

plot(relNet.chi.gamma1.0, cutPval=NULL, cutCor=.5, main="gamma 1.0",ylab="g=1")
plot(relNet.chi.gamma0.8, cutPval=NULL, cutCor=.5, main="gamma 0.8",ylab="g=0.8")
plot(relNet.chi.gamma0.5, cutPval=NULL, cutCor=.5, main="gamma 0.5",ylab="g=0.5")
plot(relNet.chi.gamma0.3, cutPval=NULL, cutCor=.5, main="gamma 0.3",ylab="g=0.3")
plot(relNet.chi.gamma0.0, cutPval=NULL, cutCor=.5, main="gamma 0.0",ylab="0.0")
dev.off()
