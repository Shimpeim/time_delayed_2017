#func.directry <- '/Users/mos/Dropbox/D/Thesis001/PG/Functions/'
#outprefix     <- 'my.pvclust_20151122'


packages.my.pvclust <- c(
  'matrixcalc'
) 
new.packages.my.pvclust <-
  packages.my.pvclust[!(packages.my.pvclust %in% installed.packages())] 

require('matrixcalc')


#=================================================================================
#   my.msfit

# check singularity of the  crossprod(X, X/vv),
# where X <- cbind(sqrt(r), 1/sqrt(r)),
# and   zz <- -qnorm(bp)
#       ,vv <- ((1 - bp) * bp)/(dnorm(zz)^2 * nboot)

#=================================================================================

my.msfit <- function (bp, r, nboot){

  if (length(bp) != length(r)) 
    stop("bp and r should have the same length")
  nboot <- rep(nboot, length = length(bp))
  use <- bp > 0 & bp < 1
  p <- se <- c(0, 0)
  names(p) <- names(se) <- c("au", "bp")
  coef <- c(0, 0)
  names(coef) <- c("v", "c")
  a <- list(p = p, se = se, coef = coef, df = 0, rss = 0, pchi = 0)
  class(a) <- "msfit"
  if (sum(use) < 2) {
    if (mean(bp) < 0.5) 
      a$p[] <- c(0, bp[r == 1])
    else a$p[] <- c(1, bp[r == 1])
    return(a)
  }
  bp <- bp[use]
  r <- r[use]
  nboot <- nboot[use]
  zz <- -qnorm(bp)
  vv <- ((1 - bp) * bp)/(dnorm(zz)^2 * nboot)
  a$use <- use
  a$r <- r
  a$zz <- zz
  X <- cbind(sqrt(r), 1/sqrt(r))
  dimnames(X) <- list(NULL, c("v", "c"))
  fit <- lsfit(X, zz, 1/vv, intercept = FALSE)
  a$coef <- coef <- fit$coef
  h.au <- c(1, -1)
  h.bp <- c(1, 1)
  z.au <- drop(h.au %*% coef)
  z.bp <- drop(h.bp %*% coef)
  a$p["au"] <- pnorm(-z.au)
  a$p["bp"] <- pnorm(-z.bp)
  if(!is.singular.matrix(crossprod(X, X/vv))){
    V <- solve(crossprod(X, X/vv))
    vz.au <- drop(h.au %*% V %*% h.au)
    vz.bp <- drop(h.bp %*% V %*% h.bp)
  }else{
    vz.au <- NA
    vz.bp <- NA
  }

  a$se["au"] <- dnorm(z.au) * sqrt(vz.au)
  a$se["bp"] <- dnorm(z.bp) * sqrt(vz.bp)
  a$rss <- sum(fit$residual^2/vv)
  if ((a$df <- sum(use) - 2) > 0) {
    a$pchi <- pchisq(a$rss, lower.tail = FALSE, df = a$df)
  }
  else a$pchi <- 1
  return(a)
}

#=================================================================================
#   my.pvrect

# modified:FUNC:pvclust::pvrect
# for picking up clusters satisfying AU *AND* BP criteria
#=================================================================================

my.pvrect <- function ( 
  x, alpha.au = 0.95, alpha.bp = 0.80, pv = "au", type = "geq", max.only = TRUE){
  hc2split <- function(x){
    A <- x$merge # (n-1,n) matrix
    n <- nrow(A) + 1
    B <- list()
    
    for(i in 1:(n-1)){
      ai <- A[i,1]
      if(ai < 0) B[[i]] <- -ai
      else B[[i]] <- B[[ai]]        
      ai <- A[i,2]
      if(ai < 0) B[[i]] <- sort(c(B[[i]],-ai))
      else B[[i]] <- sort(c(B[[i]],B[[ai]]))
    }
    
    CC <- matrix(rep(0,n*(n-1)),nrow=(n-1),ncol=n)
    
    for(i in 1:(n-1)){
      bi <- B[[i]]
      m <- length(bi)
      for(j in 1:m) CC[i,bi[j]] <- 1
    }
    split <- list(pattern=apply(CC,1,paste,collapse=""), member=B)
    return(split)
  }
  
  len <- nrow(x$edges)
  member <- hc2split(x$hclust)$member
  order <- x$hclust$order
  ht <- c()
  a <- list(clusters = list(), edges = c())
  j <- 1
  if (is.na(pm <- pmatch(type, c("geq", "leq", "gt", "lt")))) 
    stop("Invalid type argument: see help(pickup)")
  for (i in (len - 1):1) {
    if (pm == 1) 
      wh <- (x$edges[i, "au"] >= alpha.au & x$edges[i, "bp"] >= alpha.bp)
    else if (pm == 2) 
      wh <- (x$edges[i, "au"] <= alpha.au & x$edges[i, "bp"] <= alpha.bp)
    else if (pm == 3) 
      wh <- (x$edges[i, "au"] > alpha.au & x$edges[i, "bp"] > alpha.bp)
    else if (pm == 3) 
      wh <- (x$edges[i, "au"] > alpha.au & x$edges[i, "bp"] > alpha.bp)
    else if (pm == 4) 
      wh <- (x$edges[i, "au"] < alpha.au & x$edges[i, "bp"] < alpha.bp)
    if (wh) {
      mi <- member[[i]]
      ma <- match(mi, order)
      if (max.only == FALSE || (max.only && sum(match(ma, 
                                                      ht, nomatch = 0)) == 0)) {
        a$clusters[[j]] <- x$hclust$labels[mi]
        a$edges <- c(a$edges, i)
        j <- j + 1
      }
      ht <- c(ht, ma)
    }
  }
  a$edges <- a$edges[length(a$edges):1]
  
  a$clusters <- a$clusters[length(a$edges):1]
  return(a)
}

#=================================================================================
#   silent.pvclust

# 
#=================================================================================


silent.pvclust <- function(
  data,methodHclust,methodDist,nboot=nboot){
  pv <- pvclust(
    t(annMstScoreMat),
    method.hclust=methodHclust,
    method.dist=methodDist,
    nboot=nboot,
    quiet = TRUE
  )
  return(pv)
}

#save(
#  msfit,
#  my.pvrect,
#  file = sprintf('%s%s_func.RData',
#                          func.directry,
#                          outprefix)
#           )
