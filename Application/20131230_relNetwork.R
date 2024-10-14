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



for(i in 1:8) {
  relNet2TGF(result[[i]], dir="mamafiguras_20131230/", pValue=NULL, corC=.5, 
       filename=paste("RelNet_gamma_", round(gamas[i], 2), ".tgf", sep=""))
}
relNet2TGF(result[[9]], dir="mamafiguras_20131230/", pValue=NULL, corC=.5, 
           filename="RelNet_spearman.tgf")

