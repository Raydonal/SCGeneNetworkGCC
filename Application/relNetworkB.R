## Function to construct relevance networks for 1 biological type1 (Butte's
## Relevance Network)
##
##
## Parameters: data       -> Object of class maiges
##             gLabelID   -> Identification of gene name label ID.
##             sLabelID   -> sample label identification to be used
##             geneGrp    -> character string (or numeric index) specifying the
##                           gene group where to calculate the correlation
##                           values between them. If NULL (together with path)
##                           all genes are used.
##             path       -> character string (or numeric index) specifying the
##                           pathway where to calculate the correlation values
##                           between them. If NULL (together with geneGrp) all
##                           genes are used.
##             samples    -> a character vector specifying the group to be
##                           compared
##             type       -> type of measure to be calculated. May be
##                           'Rpearson' (default), 'pearson', 'kendall',
##                           'spearman' or 'MI'
##             bRep       -> number of bootstraps for correlation values
##             ...        -> additional parameters for functions robustCorr or
##                           cor
##
##
##
##    Gustavo H. Esteves & Raydonal Ospina
##
##


relNetworkB <- function(data=NULL, gLabelID="GeneName",
sLabelID="Classification", geneGrp=NULL, path=NULL, samples=NULL,
type="Rpearson", bRep=1000, ...) {
    
    
    ## Doing a simple test...
    if(is.null(data))
        stop("The data MUST be specified!!!")
    
    if(!(type %in% c("Rpearson", "pearson", "kendall", "spearman", "MI",
                     "chinchilli")))
        stop("Parameter 'type' must be 'Rpearson', 'pearson', 'kendall',
        'spearman', 'MI' or 'chinchilli'.")
    
    if(!is.null(geneGrp) & !is.null(path))
        stop("You must specify only one of geneGrp and path, or leave
        both NULL.")
    
    
    ## Getting gene and sample labels
    allGenes <- getLabels(data, gLabelID, FALSE)
    allSamples <- getLabels(data, sLabelID)
    
    
    ## removing genes and samples not named or NA from the data object
    idxGtmp <- which(allGenes != "" & !is.na(allGenes))
    idxStmp <- which(allSamples != "" & !is.na(allSamples))
    data <- data[idxGtmp, idxStmp]
    
    
    ## atualizing all gene and sample labels
    allGenes <- getLabels(data, gLabelID, FALSE)
    allGenes[data@BadSpots] <- paste(allGenes[data@BadSpots], "(*)")
    allSamples <- getLabels(data, sLabelID)
    
    
    ## Constructing the table to be used
    table <- calcW(data)
    if(is.null(geneGrp) & is.null(path))
        genes <- allGenes
    else if(!is.null(geneGrp)) {
        if(!is.numeric(geneGrp)) {
            geneGrp <- which(colnames(data@GeneGrps) == geneGrp)
            genes <- allGenes[data@GeneGrps[, geneGrp]]
        }
        else
            genes <- allGenes[data@GeneGrps[, geneGrp]]
    }
    else if(!is.null(path)) {
        genePath <- which(nodes(data@Paths[[path]]) %in% getLabels(data,
        data@Paths[[1]], FALSE))
        
        genes <- allGenes[genePath]
    }
    
    
    ## Testing the length of genes
    if (length(genes) < 2)
        stop("There is less than 2 elements in this group of genes.")
    
    
    ##  Getting the first 2 sample types to use in case of NULL samples parameter
    if(is.null(samples))
        samples <- unique(allSamples)[1]
    
    
    ## Finding indexes of the samples specified and respective lengths
    idxTmp <- allSamples %in% samples
    table <- table[, idxTmp]
    sTypes <- allSamples[idxTmp]
    
    
    ## Constructing the table to be used, averaging eventual replicates...
    tableTmp <- NULL
    for(i in unique(genes)) {
        idx <- allGenes %in% i
        if(sum(idx) > 1)
            tableTmp <- rbind(tableTmp, apply(table[idx, ], 2, mean,
            na.rm=TRUE))
        
        else
            tableTmp <- rbind(tableTmp, table[idx, ])
    }
    
    
    ## Naming tableTmp
    rownames(tableTmp) <- unique(genes)
    colnames(tableTmp) <- sTypes
    
    
    ## Finding indexes of the samples specified and respective lengths
    sample1 <- which(sTypes %in% samples)
    ##n1 <- length(sample1)
    
    
    ## Calculating the correlation coefficients
    if(type == "Rpearson") {
        corObj <- robustCorr(tableTmp[, sample1], ...)[[1]]
        pValue <- bootstrapCor(tableTmp[, sample1], bRep=bRep, type=type,
        "p-value")
        bootM <- bootstrapCor(tableTmp[, sample1], bRep=bRep, type=type, "max")
    }
    else if(type == "MI") {
        corObj <- MI(tableTmp[, sample1], ...)
        pValue <- bootstrapMI(tableTmp[, sample1], bRep=bRep, "p-value")
        bootM <- bootstrapMI(tableTmp[, sample1], bRep=bRep, "max")
    }
    else if(type == "chinchilli") {
        corObj <- corChinchilli(t(tableTmp[, sample1]), ...)
        pValue <- bootstrapCor(tableTmp[, sample1], bRep=bRep, type=type,
                              "p-value")
        bootM <- bootstrapCor(tableTmp[, sample1], bRep=bRep, type=type, "max")      
    }
    else {
        corObj <- cor(t(tableTmp[, sample1]), method=type, ...)
        pValue <- bootstrapCor(tableTmp[, sample1], bRep=bRep, type=type,
        "p-value")
        bootM <- bootstrapCor(tableTmp[, sample1], bRep=bRep, type=type, "max")
    }
    
    
    ## Picking R and packages version information
    tmp <- sessionInfo()
    vInfo <- list()
    vInfo$R.version <- tmp$R.version$version.string
    vInfo$BasePacks <- tmp$basePkgs
    tmp1 <- NULL
    for (i in 1:length(tmp$otherPkgs))
        tmp1 <- c(tmp1, paste(tmp$otherPkgs[[i]]$Package, "version",
        tmp$otherPkgs[[i]]$Version))
    
    vInfo$AddPacks <- tmp1
    
    
    ## Defining the object to return
    result <- new("maigesRelNetB", W=tableTmp, Corr=corObj, Pval=pValue,
    maxB=bootM, Date=date(), type=paste(samples, collapse=", "),
    Slabel=sLabelID, V.info=vInfo)
    
    return(result)
    
}
