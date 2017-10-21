#==================================
# renamed from 'OldFunc_20151125.R'
#==================================

#```{r}
makeDifData <- function(data,var){
  for (i in 2:length(var)){
    end   <- var[i]
    start <- var[i-1]
    rescol <- paste('d',end,sep='_')
    ddata  <- data[,end] - data[,start]
    data[,rescol] <- ddata
  }
  return(data)
}
#```



#```{r defFUNC det_out (determinater of outlier)}
det_out <- function(q_dat,hinge){
  IQR       <- q_dat[4,]  - q_dat[2,]
  upp_lim   <- q_dat[4,]  + hinge * (IQR)
  low_lim   <- q_dat[2,]  - hinge * (IQR)
  outBounds <- data.frame(upp_lim,low_lim)
  return(outBounds) 
}
#```

#Define 'out_detect' function : OUTPUT the value is outlier or not as 0/1 under the input criteria 'outBounds' 
#```{r detecting outlier}
outDetectUpp <- function( data, uppBounds)
{
  if (data > uppBounds)
  {return(c(data,1))}else return(c(data,0)) 
}

outDetectLow <- function( data, lowBounds)
{
  if (data < lowBounds)
  {return(c(data,1))}else return(c(data,0)) 
}
#```

#```{r defining corPrt(function)}
corPrt <- function(i,j,timeLag,indata,rep=repet,timePoint=timePoint){
  outdata <- list()
  IDx <- i
  IDy <- j
  
  w_PrtX   <- indata %>% filter(prtNo==IDx)
  w_PrtY   <- indata %>% filter(prtNo==IDy)
  
  # make a key variable 'seq', for joinning x_data and y_data.
  # For y_data, 'w_seq' is set interimly,  and adjust it to 'seq'
  # after deletion of the first timepoint observation of y_data.
  
  w_PrtX$seq   <- c(1:nrow(w_PrtX)) 
  w_PrtY$w_seq <- c(1:nrow(w_PrtY))
  
  PrtX <- w_PrtX %>% 
    filter(seq <= rep*(timePoint-1)-timeLag) %>%
    rename(xVal1=vec )
  
  PrtY <- w_PrtY %>% 
    mutate(seq=w_seq-timeLag) %>%
    rename(yVal1=vec )
  
  outdata$PrtTab     <- PrtX %>% inner_join(PrtY,by='seq')
  outdata$PrtXtab    <- as.data.frame(ftable(yVal1~xVal1,data=outdata$PrtTab))
  
  return(outdata)
}
#```