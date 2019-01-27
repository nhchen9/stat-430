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



tstar_tib <- function(data, w0=10, nbt0=20, nT=15)
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

tstar_vrb <- function(data, w0 = 10, bkw_T = 10, bkw_Pb1 = 10){
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

tib = tstar_tib(dat_exc)
vib = tstar_vib(dat_exc)
trb = tstar_trb(dat_exc)
vrb = tstar_vrb(dat_exc)

plot(dat_exc$Price, pch=20, xlab="ticks", ylab="Price", main="Tick Runs bars")
abline(v=cumsum(trb), lwd=0.2)
