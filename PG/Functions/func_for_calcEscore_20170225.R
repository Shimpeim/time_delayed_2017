
#===================================================================
#    gq_method

#  Qian,2001 
#===================================================================

gq_method <- function(vec_i,vec_j,timepoint,rep=1,pos_neg='pos'){
  gq   <- list()
  for(rep in 1:rep){
    vec_i <- c(append(
      c(vec_i),0,after=((timepoint+1)*(rep-1))
    )
    )
    vec_j <- c(append(
      c(vec_j),0,after=((timepoint+1)*(rep-1))
    )
    )
  }
  veci <- vec_i
  vecj <- vec_j
  
  timePoint <- timepoint+1
  tp_i <- length(veci)
  tp_j <- length(vecj)
  
  mat   <- matrix(rep(0,tp_i*tp_j),ncol=tp_j) # initiate the matrix
  mat_M <- matrix(rep(0,tp_i*tp_j),ncol=tp_j) # initiate the matrix
  
    
  for(i in 2:length(veci)){
    for(j in 2:length(vecj)){
      
      if(pos_neg == 'pos'){mat[i,j] <- mat[i-1,j-1] + as.numeric(veci[i])*as.numeric(vecj[j]) }
      if(pos_neg == 'neg'){mat[i,j] <- mat[i-1,j-1] - as.numeric(veci[i])*as.numeric(vecj[j]) }
      if (mat[i,j]<0) {mat[i,j] <- 0}
      
      if( is.element(i-1,c(timePoint* 1:(rep-1)))  || is.element(j-1, c(timePoint* 1:(rep-1)))) {
        mat[i,j] <- 0        
      }
    }
  }
  for(i in 2:length(veci)){
    for(j in 2:length(vecj)){
      
      if(pos_neg == 'pos'){mat_M[i,j] <- as.numeric(veci[i])*as.numeric(vecj[j]) }
      if(pos_neg == 'neg'){mat_M[i,j] <- as.numeric(veci[i])*as.numeric(vecj[j]) }

      if( is.element(i-1,c(timePoint* 1:(rep-1)))  || is.element(j-1, c(timePoint* 1:(rep-1)))) {
        mat_M[i,j] <- 0        
      }
    }
  }
  
  posix <- ceiling(which.max(t(mat))/tp_j)
  posiy <- ceiling(which.max(mat)/tp_i)
  
  posixVec <- c()
  posiyVec <- c()
  maxVec   <- c()
  for(k in 1:rep){
    start <- 1 + (k-1)*timePoint
    edge  <- k*timePoint
    k_mat <- mat[start:edge,start:edge]
    k_posix <- ceiling(which.max(t(k_mat))/tp_j) ; k_posiy <- ceiling(which.max(k_mat)/tp_i)
    posixVec[k] <- k_posix ; posiyVec[k] <- k_posiy
    maxVec[k]   <- max(k_mat)
  }
  posiDisVec <-  posixVec-posiyVec
  
  gq$mat   <- data.frame(mat)
  gq$mat_M   <- data.frame(mat_M)
  gq$index <- c(posix,posiy)
  gq$score <- c(posix-posiy,max(mat)) 
  gq$subPosiDisVec <- posiDisVec
  gq$subMaxVec     <- maxVec
  
  return(gq)
}


gq_method2 <- function(data){
  vec_i <- data[data$dtname=='prt',3:length(data)]
  vec_j <- data[data$dtname=='trc',3:length(data)]
  pos <- gq_method(vec_i,vec_j,timepoint=timePoint,rep=repet,'pos')
  neg <- gq_method(vec_i,vec_j,timepoint=timePoint,rep=repet,'neg')
  res <- data.frame(t(pos$score),t(neg$score))
  #dimnames(res[[2]]) <- c('id','posit_Emax','val_Emax','posit_Dmax','val_Dmax')
  return(res)
}

#===================================================================
# gq_scoreMatOut
# gq_scoreMatOut_withName

#
#
#===================================================================
gq_scoreMatOut <- function(data){
  vec_i <- data[1,3:length(data)]
  vec_j <- data[2,3:length(data)]
  pos <- gq_method(vec_i,vec_j,timepoint=timePoint,rep=repet,'pos')
  mat <- t(c(as.matrix(pos$mat)))#,c(1:length(as.matrix(pos$mat))))  
  #  mat <- as.data.frame(t(c(as.matrix(pos$mat))))#,c(1:length(as.matrix(pos$mat))))
  return(mat)
}

gq_scoreMatOut_withName <- function(data){
  vec_i <- data[1,3:length(data)]
  vec_j <- data[2,3:length(data)]
  pos <- gq_method(vec_i,vec_j,timepoint=timePoint,rep=repet,'pos')
  w.mat <- t(c(as.matrix(pos$mat)))
  mat   <- matrix(w.mat,nrow=timePoint+1,byrow=TRUE)
  out <- list()
  out$mat <- mat
  out$id  <- data$id
  return(out)
}

#===================================================================
#    my.heatmap.3

#  input data from FUNC:gq_scoreMatOut_withName(data), 
#  data in which 'id'  
#  when used in FUNC:llply(.(id),arg:func), 
#===================================================================
my.heatmap.3 <- function(data){
  mat   <- matrix(data$mat,nrow=timePoint+1,byrow=TRUE)
  out   <- heatmap.3(
    mat,
    trace="none", 
    dendrogram="none", 
    Rowv=F,
    Colv=F, 
    color.FUN="redgreen", 
    cluster.by.row=F,
    cluster.by.col=F,
    mapratio=1,
    mapsize=4,
    main=data$id
    )
  return(out)
}
