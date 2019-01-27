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

vrb_theta <- function(btvt, bt){
  vplus = rep(0, length(bt))
  vminus = vplus
  vrb_theta = data.frame(vplus,vminus)
  for( i in 1: length(bt)){
    if(bt[i] ==1){
      vrb_theta$vplus[i] = vrb_theta$vplus[i-1] + btvt[i]
      vrb_theta$vminus[i] = vrb_theta$vminus[i-1]
    }
    if(bt[i] == -1){
      vrb_theta$vplus[i] = vrb_theta$vplus[i-1]
      vrb_theta$vminus[i] = vrb_theta$vminus[i-1] - btvt[i]
    }
  }
  
  return (vrb_theta)
}

tstar_vrb <- function(data, w0 = 10, bkw_T = 5, bkw_Pb1 = 5){
  btvt = volume_imbalance(data)
  b_t = tick_imbalance(data)
  theta = vrb_theta(btvt, b_t)
  
  l = length(b_t)
  w0 = max(min(which(theta$vplus != 0)), w0)
  w0 = max(min(which(b_t ==1)), w0)
  tvec = w0
  e0t = tvec
  last = tvec
  pb1 = sum(b_t[1:w0]==1)/w0
  btvt_temp = btvt[1:w0]
  e0vtbp = sum(btvt_temp[which(btvt_temp < 0)])/length(which(btvt_temp[1:w0] < 0))
  e0vtbm = sum(btvt_temp[which(btvt_temp > 0)])/length(which(btvt_temp[1:w0] > 0))
  
  e0vtbpvec = e0vtbp
  e0vtbmvec = e0vtbm
  
  pb1vec = pb1
  
  theta_expected = e0t * max(pb1*e0vtbp, abs((1-pb1)*e0vtbm))
  print('first theta_exp')
  print(theta_expected)
  print(e0vtbp)
  print(e0vtbm)
  
  while(last < l){
    new = FALSE
    last= sum(tvec)
    for(i in 1:(l-last)){
      b_t_tmp = b_t[(last+1):(last+i)]
      btvt_tmp = btvt[(last+1):(last+i)]
      theta_tmp <- max(cumsum(btvt_tmp[b_t_tmp==1]), -cumsum(btvt_tmp[b_t_tmp==-1]))
      if(theta_tmp > theta_expected){
        t_new = i
        tvec = c(tvec, t_new)
        last = sum(tvec)
        pb1_new = sum(b_t_tmp==1) / i
        pb1vec = c(pb1vec, pb1_new)
        
        e0vtbp_new = sum(btvt_tmp[b_t_tmp==1]) / max(length(which(b_t_tmp > 0)), 1)
        e0vtbm_new = sum(btvt_tmp[b_t_tmp==-1]) / max(length(which(b_t_tmp < 0)), 1)
        
        e0vtbpvec = c(e0vtbpvec, e0vtbp_new)
        e0vtbmvec = c(e0vtbmvec, e0vtbm_new)
        
        new = TRUE
        break
      }
    }
    
    if(new){
      ltv = length(tvec)
      if(ltv <=2){
        e0t = mean(tvec)
        pb1 = mean(pb1vec)
        e0vtbp = mean(e0vtbpvec)
        e0vtbm = mean(e0vtbmvec)
      }else{
        nT = min(bkw_T, ltv-1)
        e0t = tail(pracma::movavg(tvec[(ltv-nT):ltv], n= nT, type = 'e'),1)
        npb = min(bkw_Pb1, ltv-1)
        pb1 = tail(pracma::movavg(pb1vec[(ltv-nT):ltv], n = nT, type = 'e'),1)
        e0vtbp = mean(e0vtbpvec)
        e0vtbm = mean(e0vtbmvec)
        #e0vtbp = tail(pracma::movavg(e0vtbpvec[(ltv-nT):ltv], n = nT, type = 'e'),1)
        #e0vtbm = tail(pracma::movavg(e0vtbmvec[(ltv-nT):ltv], n = nT, type = 'e'),1)
      }
      theta_expected = e0t * max(pb1*e0vtbp, abs((1-pb1)*e0vtbm))
      #print('new iteration')
      #print(ltv)
      #print(e0vtbp)
      #print(e0vtbm)
      #print(theta_expected)
    }else{
      break
    }
    
  }
  
  return(tvec)
  
}
x = tstar_vrb(dat_exc)

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

plot(dat_exc$Price, pch=20, xlab="ticks", ylab="Price", main="Where to sample tick runs bars?")
abline(v=cumsum(y), lwd=0.2)
