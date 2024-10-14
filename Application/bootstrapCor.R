## Function to calculate robust correlation coeficients (minimizing the
## influences of one outlier). This strategy uses an idea similar to
## leave-one-out of cross-validation, where we make a decision based on the
## "rule of three numbers".
##
## Parameters: x           -> matrix to calculate the bootstraped p-values of
##                            correlations (or maximum boot cor) between all
##                            pairs of rows from x. Return a square matrix with
##                            all p-values correlations (or maximum). If x is a
##                            matrix, the function calculates a matrix of
##                            p-values or maximuns, else x must be a vector and
##                            y must be specified
##             y           -> optional numerical vector. Must be specified if x
##                            is a vector. If x is a matrix, y is ignored
##             bRep        -> number of bootstraps samples (by permutation) to
##                            generate
##             type        -> Character string specifying the type of
##                            correlation statistic to be used. Possible values:
##                            Rpearson, pearson, spearman or kendall.
##             ret         -> of value to return? Must be 'p-value' ou 'max'
##             alternative -> type of test to do, must be 'two.sided' (default),
##                            'less' or 'greater'
##
## Gustavo H. Esteves
## 14/05/07
##
##


bootstrapCor <- function(x, y=NULL, bRep, type="Rpearson", ret="p-value",
alternative="two.sided", ...) {
    
    
    ## Doing basic tests
    if(!(type %in% c("Rpearson", "pearson", "kendall", "spearman",
                     "chinchilli")))
        stop("Parameter 'type' must be 'Rpearson', 'pearson', 'kendall',
        'spearman' or 'chinchilli'.")
    
    if(!(ret %in% c("p-value", "max")))
        stop("Parameter 'ret' must be 'p-value' or 'max'.")
    
    if(!(alternative %in% c("two.sided", "less", "greater")))
        stop("Parameter 'alternative' must be 'two-sided', 'less' or
        'greater'.")
    
    
    if(length(dim(x)) == 2) {
        
        ## Getting the dimensions of the matrix
        ns <- dim(x)
        
        ## Doing the cor calculations
        tmp <- array(NA, c(dim(x)[1], dim(x)[1], bRep))
        if (type == "Rpearson") {
            orig <- robustCorr(x)[[1]]
            for(j in 1:bRep) {
                
                if(ret == "max")
                    tmp[, , j] <- abs(robustCorr(t(apply(x, 1, sample)))[[1]])
                
                else if((ret == "p-value") & (alternative == "two.sided"))
                    tmp[, , j] <-
                    abs(robustCorr(t(apply(x, 1, sample)))[[1]]) >= abs(orig)
                
                else if((ret == "p-value") & (alternative == "less"))
                    tmp[, , j] <-
                    robustCorr(t(apply(x, 1, sample)))[[1]] <= orig
                
                else if((ret == "p-value") & (alternative == "greater"))
                    tmp[, , j] <-
                    robustCorr(t(apply(x, 1, sample)))[[1]] >= orig
                
            }
        }
        else if (type == "chinchilli") {
          orig <- corChinchilli(t(x), ...)
          
          for(j in 1:bRep) {
            
            if(ret == "max")
              tmp[, , j] <- abs(corChinchilli(apply(x, 1, sample), ...))
            
            else if((ret == "p-value") & (alternative == "two.sided"))
              tmp[, , j] <-
              abs(corChinchilli(apply(x, 1, sample), ...)) >= abs(orig)
            
            else if((ret == "p-value") & (alternative == "less"))
              tmp[, , j] <-
              corChinchilli(apply(x, 1, sample), ...) <= orig
            
            else if((ret == "p-value") & (alternative == "greater"))
              tmp[, , j] <-
              corChinchilli(apply(x, 1, sample), ...) >= orig
            
          }
        }
        else {
            orig <- cor(t(x), method=type)
            
            for(j in 1:bRep) {
                
                if(ret == "max")
                    tmp[, , j] <- abs(cor(apply(x, 1, sample), method=type))
            
                else if((ret == "p-value") & (alternative == "two.sided"))
                    tmp[, , j] <-
                    abs(cor(apply(x, 1, sample), method=type)) >= abs(orig)
            
                else if((ret == "p-value") & (alternative == "less"))
                    tmp[, , j] <-
                    cor(apply(x, 1, sample), method=type) <= orig
            
                else if((ret == "p-value") & (alternative == "greater"))
                    tmp[, , j] <-
                    cor(apply(x, 1, sample), method=type) >= orig

            }
        }
        
        if(ret == "p-value")
            ## Calculating p-value
            res <- apply(tmp, c(1, 2), sum)/bRep
        else if(ret == "max")
            ## Calculating maximum values
            res <- apply(tmp, c(1, 2), max)
        
        return(res)
    }
    else if(is.vector(x) & is.vector(y)) {
        tmp <- NULL
        if (type == "Rpearson") {
            orig <- robustCorr(x, y)[[1]]
            for(j in 1:bRep) {
                tmp <- c(tmp, robustCorr(sample(x), sample(y))[[1]])
            }
        }
        else if (type == "chinchilli") {
          orig <- corChinchilli(x, y, ...)
          for(j in 1:bRep) {
            tmp <- c(tmp, corChinchilli(sample(x), sample(y), ...))
          }
        }
        else {
            orig <- cor(x, y, method=type)
            for(j in 1:bRep) {
                tmp <- c(tmp, cor(sample(x), sample(y), method=type))
            }
        }
        
        ## Calculating p-value
        if(ret == "max")
            res <- max(abs(tmp))
        else if((ret == "p-value") & (alternative == "two.sided"))
            res <- sum(abs(tmp) >= abs(orig))/bRep
        else if((ret == "p-value") & (alternative == "less"))
            res <- sum(tmp <= orig)/bRep
        else if((ret == "p-value") & (alternative == "greater"))
            res <- sum(tmp >= orig)/bRep
        
        return(res)
        
    }
    else
        stop("x must be a numerical matrix or vector, case where you must
        specify another vector y of same length.")
    
}
