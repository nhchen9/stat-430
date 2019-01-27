library(lubridate)
options(digits.secs=3)
setwd("~/Development/stat 430")
dat <- read.csv("./datasets/AMZN_2012-06-21_34200000_57600000_message_10.csv", header = F)
names(dat) <- c("Time", "Type", "OrderID", "Size", "Price", "Direction")

dat$Size <- as.numeric(dat$Size)
dat$Price <- as.numeric(dat$Price)

demodate <- "2012-06-21"
dat$tStamp <- as_datetime(demodate, tz="US/Eastern") + dat$Time
dat_exc <- subset(dat, Type %in% c(4,5))

tick_imbalance <- function(data)
{
  p_d = diff(data$Price)
  ticks = rep(0, (length(data$Price)))
  for( i in 1:(length(data$Price)-1)){
    ticks[i+1] = ifelse(p_d[i] == 0, ticks[i], sign(p_d[i]))
  }
  return(ticks)
}

tstar_trb <- function(dat, w0=10, bkw_T = 5, bkw_Pb1 = 5)
{
  b_t = tick_imbalance(dat)
  l = length(b_t)
  th_T = sapply(1:l, function(i){
    b_t_tmp = b_t[1:i]
    if (sum(b_t_tmp %in% c(-1, 1)) == 0){out = 0}else{
      out = max(cumsum(b_t_tmp[b_t_tmp == 1]), -cumsum(b_t_tmp[b_t_tmp==-1]))
    }
    out
  })
  
  w0 = max(min(which(th_T != 0)), w0)
  w0 = max(min(which(b_t ==1)), w0)
  tvec = w0
  e0t = tvec
  last = tvec
  pb1 = sum(b_t[1:w0]==1)/w0
  pb1vec = pb1
  th_T_expected = e0t*max(pb1, 1-pb1)
  new = FALSE
  while(last < l){
    #print(last)
    new = FALSE
    last = sum(tvec)

    for(i in 1:(l-last)){
      b_t_tmp = b_t[(last+1):(last+i)]
      th_T_tmp <- max(cumsum(b_t_tmp[b_t_tmp==1]), -cumsum(b_t_tmp[b_t_tmp==-1]))
        if (th_T_tmp > th_T_expected){
          t_new = i
          tvec = c(tvec, t_new)
          last = sum(tvec)
          pb1_new = sum(b_t_tmp==1) / i
          pb1vec = c(pb1vec, pb1_new)
          new = TRUE
          break
        }
      }
    
    
    if(new){
      ltv = length(tvec)
      if(ltv <=2){
        e0t = mean(tvec)
        pb1 = mean(pb1vec)
      }else{
        nT = min(bkw_T, ltv-1)
        e0t = tail(pracma::movavg(tvec[(ltv-nT):ltv], n= nT, type = 'e'),1)
        npb = min(bkw_Pb1, ltv-1)
        pb1 = tail(pracma::movavg(pb1vec[(ltv-nT):ltv], n = nT, type = 'e'),1)
      }
      th_T_expected = e0t*max(pb1, 1-pb1)
    }else{
      break
    }
    
  }
  
  return(tvec)
  
}
y = tstar_trb(dat_exc)

bar_tick_imbalance <- function(dat, w0=10, bkw_T=5, bkw_b=5)
{
  T_tib <- tstar_tib(dat)
  T_tib <- c(T_tib, nrow(dat)-sum(T_tib)) # the remaining data is treated as a bar
  
  winEnd <- cumsum(T_tib)
  winIdx <- as.factor(unlist(sapply(1:length(winEnd), function(i){rep(winEnd[i], T_tib[i])})))
  
  H <- stats::aggregate(dat$Price, by = list(winIdx), max)$x
  L <- stats::aggregate(dat$Price, by = list(winIdx), min)$x
  O <- stats::aggregate(dat$Price, by = list(winIdx), function(x){x[1]})$x
  C <- stats::aggregate(dat$Price, by = list(winIdx), function(x){x[length(x)]})$x
  V <- stats::aggregate(dat$Size, by = list(winIdx), sum)$x
  list(H=H,L=L,O=O,C=C,V=V)
}
library(fmlr)
fmlr::bar_tick_imbalance(dat_exc)
bar_tick_imbalance(dat_exc)
