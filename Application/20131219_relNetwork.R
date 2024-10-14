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

## Definindo os valores de gama 
gamas <- seq(1, 0, length.out=8)

## Definindo uma lista para armazenar as redes de relevância
result <- list()

for(i in 1:8)
  result[[i]] <- relNetworkB(sarahSummSingle[,idx],
                             sLabelID="DS_DIAG_PROJ", geneGrp=NULL,
                             path='Renin-angiotensin system',
                             type="chinchilli", gamma=gamas[i])

## Fazendo a rede para Spearman
result[[9]] <- relNetworkB(sarahSummSingle[,idx],
                                   sLabelID="DS_DIAG_PROJ", geneGrp=NULL,
                                   path='Renin-angiotensin system',
                                   type="spearman")



png("figuras_20131219/estomago_normal_cor_pvalor.png", width=1200, height=1200)
#pdf("figuras_20131219/estomago_normal_cor_pvalor.pdf", width=1200, height=1200)
par(mfrow=c(3,3))
for(i in 1:8) {
  plot(result[[i]], cutPval=NULL, cutCor=.5, 
       name=paste("gamma", round(gamas[i], 2), "|cor| > 0.5"))
}
plot(result[[9]], cutPval=NULL, cutCor=.5, name="|Spearman| > 0.5")
dev.off()
