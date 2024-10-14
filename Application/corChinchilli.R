## Function to calculate correlation coeficients from Chinchilli et al
## (2005).
##
## Gustavo H. Esteves & Raydonal O. Martinez
## sáb 12 out 2024 19:42:16 -03

##


corChinchilli <- function(x, y=NULL, gamma=1) {
  ## Parameters: x -> matrix to calculate the correlations between all
  ##                  pairs of rows from x. Return a square matrix with
  ##                  all correlations. If x is a matrix, the function
  ##                  calculates a matrix of robust correlations, else
  ##                  x must be a vector and y must be specified as
  ##                  another vector of same length and the correlation
  ##                  between both are calculate
  ##             y -> optional numeric vector
  ##         gamma -> Parâmetro que especifica o \gamma do artigo de
  ##                  Chinchilli, et al. Padrão gamma=1 (que coincide com
  ##                  a correlação de Pearson)
  ##
  ##
  
  ## Catching the x dimension
  ns <- dim(x)
  
  if(length(ns) == 2) {
    
    ## Carregando a biblioteca implementada em C (o .so)
    dyn.load("./corChinchilli.so")
    
    ## Rodando a implementação compilada em C
    resC <- .C("corCh", as.double(as.vector(x)), as.integer(ns[1]),
               as.integer(ns[2]), as.double(gamma),
               as.double(rep(1, (ns[2])^2)))
    
    result <- matrix(resC[[5]], nrow=ns[2])
    if(!is.null(colnames(x))) {
      rownames(result) <- colnames(x)
      colnames(result) <- colnames(x)
    }
    
    
    return(result)

  }
  else if(is.vector(x) & !is.null(y)) {
    
    ## Doing simple tests
    if(!is.vector(y))
      stop("y must be a vector of same length as x!")
    if(length(x) != length(y))
      stop("x and y must have same length!")
    
    ## calculating the correlation value
    x <- cbind(x, y)
    ns <- dim(x) ## Atualiza os valores das dimensoes de x
    
    
    ## Carregando a biblioteca implementada em C (o .so)
    dyn.load("../../chinchilli/corChinchilli.so")
    
    ## Rodando a implementação compilada em C
    resC <- .C("corCh", as.double(as.vector(x)), as.integer(ns[1]),
               as.integer(ns[2]), as.double(gamma),
               as.double(rep(1, (ns[2])^2)))
    
    #return(matrix(resC[[5]], nrow=ns[2]))
    return(resC[[5]][2])
    
  }
  else
    stop("x must be a numerical matrix or a vector. In this case you must
        specify another vector y of same length.")

}

