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

# relNet.pearson <- relNetworkB(sarahSummSingle[,idx],
#                                    sLabelID="DS_DIAG_PROJ", geneGrp=NULL,
#                                    path='Renin-angiotensin system',
#                                    type="pearson")
  
relNet.chi.gamma1.0 <- relNetworkB(sarahSummSingle[,idx],
                                   sLabelID="DS_DIAG_PROJ", geneGrp=NULL,
                                   path='Renin-angiotensin system',
                                   type="chinchilli", gamma=1)

# relNet.chi.gamma0.8 <- relNetworkB(sarahSummSingle[,idx],
#                                    sLabelID="DS_DIAG_PROJ", geneGrp=NULL,
#                                    path='Renin-angiotensin system',
#                                    type="chinchilli", gamma=.8)

relNet.chi.gamma0.5 <- relNetworkB(sarahSummSingle[,idx],
                                   sLabelID="DS_DIAG_PROJ", geneGrp=NULL,
                                   path='Renin-angiotensin system',
                                   type="chinchilli", gamma=.5)

# relNet.chi.gamma0.3 <- relNetworkB(sarahSummSingle[,idx],
#                                    sLabelID="DS_DIAG_PROJ", geneGrp=NULL,
#                                    path='Renin-angiotensin system',
#                                    type="chinchilli", gamma=.3)

relNet.chi.gamma0.0 <- relNetworkB(sarahSummSingle[,idx],
                                   sLabelID="DS_DIAG_PROJ", geneGrp=NULL,
                                   path='Renin-angiotensin system',
                                   type="chinchilli", gamma=0)



png("figuras/estomago_normal_cor_pvalor.png", width=1200, height=800)
par(mfrow=c(2,3))
plot(relNet.chi.gamma1.0, cutPval=NULL, cutCor=.5, name="gamma 1.0, |cor| > 0.5")
plot(relNet.chi.gamma0.5, cutPval=NULL, cutCor=.5, name="gamma 0.5, |cor| > 0.5")
plot(relNet.chi.gamma0.0, cutPval=NULL, cutCor=.5, name="gamma 0.0, |cor| > 0.5")
plot(relNet.chi.gamma1.0, cutPval=0.001, name="gamma 1.0, p-value < 0.001")
plot(relNet.chi.gamma0.5, cutPval=0.001, name="gamma 0.5, p-value < 0.001")
plot(relNet.chi.gamma0.0, cutPval=0.001, name="gamma 0.0, p-value < 0.001")
dev.off()