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

volume_imbalance <- function(data)
{
  p_d = diff(data$Price)
  ticks = rep(0, (length(data$Price)))
  for( i in 1:(length(data$Price)-1)){
    ticks[i+1] = ifelse(p_d[i] == 0, ticks[i], sign(p_d[i]))
  }
  theta = rep(0, (length(data$Price)))
  for( i in 1:(length(data$Size))){
    theta[i] = ticks[i] * data$Size[i]
  }
  return(theta)
}


tstar_vib <- function(data, w0=10, vtgpbt0 = 100)
{
  b_t = tick_imbalance(data)
  theta = volume_imbalance(data)
  w0 <- max(min(which(cumsum(b_t) != 0)), w0)
  tvec = w0
  e0t = w0
  l = dim(dat)[1]
  vtgpbt = vtgpbt0
  
  
  repeat
  {
    last = sum(tvec)
    datv = dat$Size[1:last]
    datb1 = datv[which(b_t[1:last] == 1)]
    pbt = tail(pracma::movavg(b_t[1:last], n = last-1, type = "e"), 1)
    vtgpbt = tail(pracma::movavg(datb1, n = length(datb1)-1, type = 'e'), 1)
    e0vt = tail(pracma::movavg(datv, n = length(datv)-1, type= 'e'), 1)
    
    
    vplus = pbt*vtgpbt 
    
    theta_expected = abs(e0t*(2*vplus - e0vt))
    theta_sum = abs(cumsum(theta[-(1:last)]))
    #theta_sum = abs(cumsum(theta))
    if(max(theta_sum) < theta_expected){
      break
    }
    else{
      tnew = min(which(theta_sum >= theta_expected))
    }
    
    last = last + tnew
    
    if(last > l){ 
      break
    }
    else{
      tvec = c(tvec, tnew)
      if(length(tvec) <= 2){
        e0t = mean(tvec)
      }else{
        e0t = tail(pracma::movavg(tvec, n= length(tvec)-1, type = "e"),1)
      }
    }
  }
  return(tvec)
}
y = tstar_vib(dat_exc)


theta = volume_imbalance(dat_exc)
b_t = tick_imbalance(dat_exc)
tstar_tib(dat_exc)
