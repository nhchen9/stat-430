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

tstar_tib <- function(data, w0=20, nbt0=30, nT=20)
{
  b_t = tick_imbalance(data)
  w0 <- max(min(which(cumsum(b_t) != 0)), w0)
  tvec = w0
  e0t = w0
  l = length(b_t)
  
  repeat
  {
    last = sum(tvec)
    pbt = tail(pracma::movavg(b_t[1:last], n = min(nbt0, last-1), type = "e"), 1)
    #print(pbt)
    e0t = max(10, e0t)
    b_t_expected = e0t*abs(pbt)
    b_t_psum = abs(cumsum(b_t[-(1:last)]))
    #b_t_psum = abs(cumsum(b_t))
    
    if(max(b_t_psum) < b_t_expected){
      break
    }
    else{
      tnew = min(which(b_t_psum >= b_t_expected))
    }
    
    last = last + tnew
    
    if(last >= l){ 
      break
    }
    else{
      tvec = c(tvec, tnew)
      if(length(tvec) <= 2){
        e0t = mean(tvec)
      }else{
        e0t = tail(pracma::movavg(tvec, n=min(nT, length(tvec)-1), type = "e"),1)
      }
    }
  }
  return(tvec)
}

x = tstar_tib(dat_exc)

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
