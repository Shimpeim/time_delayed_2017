#=================================================================================
#   varSfl

#
#=================================================================================
varSfl <- function(data,sfl,startCol,endCol){  
  target <- sample(dimnames(data)[[2]][startCol:endCol],endCol-startCol+1)
  out <- data
  for(i in 0:(endCol-startCol)){
    if(sfl==TRUE){
      out[startCol+i] <- data[target[i+1]]
    }else{out[startCol+i] <- data[startCol+i]}    
  }
  return(out)  
}

#=================================================================================
#   cor.test_by.id

#
#=================================================================================
cor.test_by.id <- function(data,x,y,method.cor='spearman'){
  
  data <- data %>%
    filter(
      length(data[,x])>2 & length(data[,y])>2
    )
  
  cor.test.res <- cor.test(
    data[,x],data[,y],method=method.cor
  )
  pvalue <- cor.test.res$p.value
  
  rho.est <- cor.test.res$estimate
  output <- data.frame(pvalue,rho.est)#,data$V1,data$V4)
  return(output)
}

#=================================================================================
#   my.rank

#
#=================================================================================
my.rank <- function(data,var){
  data$rank <- rank(data[,var])
  return(data)
}

#=================================================================================
#   my.qtl

#
#=================================================================================
my.qtl  <- function(data,var){
  data$qtl  <- quantile(
    data[,var], prob=c(0.025,0.05,0.95,0.975),names=FALSE
  )
}