#loading the necessary libraries
library(wavelets)
library(WaveletComp)
library(cmvnorm)
library(MASS)
library(mixtools)
library(jocre)
library(matlib)
library(tidyverse)


return_df_dates<- seq(as.Date("2019-03-11"), as.Date("2022-12-30"), by="days")
WTI <- WTI %>%
  mutate(Date = as.Date(observation_date)) %>%
  complete(Date = seq(as.Date("2019-03-11"), as.Date("2022-12-30"), by="days") )
wti_price<- imputeTS::na.interpolation(WTI$DCOILWTICO,option = "stine")
wti_return<- 100*c(NA,diff(log(wti_price)))
wti_df<- data.frame(as.Date(WTI$Date),wti_return)
colnames(wti_df)<- c("Date","wti")
wti_df$wti<- imputeTS::na.interpolation(wti_df$wti,option="stine")
#wti_df<- wti_df[-1,]

#Getting the GS statistics
sigma_mat<- matrix(c(1,(-(1/sqrt(6))),(-5/sqrt(60)),(-1/sqrt(6)),1,(2/sqrt(360)),(-5/sqrt(60)),(2/sqrt(360)),1),nrow = 3)
sigma_inv<- solve(sigma_mat)
mu_vec<- as.matrix(c(0,0,0))
crit_value_z2<- qchisq(p=.05, df=3, lower.tail=FALSE) #z2=t(GS)*sigma_inv*GS #crit_value_z2=5.9915


sigtest_func<- function(svec)
{
  svect<- as.matrix(t(svec),nrow=length(svec),ncol=1)
  z2<- sum(svec*colSums(svec*t(sigma_inv)))
  return(z2)
}


square_func<- function(x){
  return(x^2)
}

test_func<- function(series)
{
  ts<- as.ts(series)
  n<- length(series)
  ts_modwt<- modwt(ts,filter = "haar",n.levels=3)
  scale1_wc<- ts_modwt@W[["W1"]] #scale 1 wavelet coefficients
  scale2_wc<- ts_modwt@W[["W2"]] #scale 2 wavelet coefficients
  scale3_wc<- ts_modwt@W[["W3"]] #scale 3 wavelet coefficients
  ss<- sum(unlist(lapply(series,square_func)))
  E1<- sum(unlist(lapply(scale1_wc,square_func)))/ss
  E2<- sum(unlist(lapply(scale2_wc,square_func)))/ss
  E3<- sum(unlist(lapply(scale3_wc,square_func)))/ss
  GS1<- sqrt(4*n)*(E1-(1/2))
  GS2<- sqrt(32*(n/3))*(E2-(1/4))
  GS3<- sqrt(256*(n/15))*(E3-(1/8))
  GS<- c(GS1,GS2,GS3)
  z2<- sigtest_func(GS)
  return(z2)
}




#now we would create rolling window of 30 days and get z2 value for each of them

library(roll)
library(runner)
id_func<- function(x)
{
  return(x)
}
rolling_sample_func<- function(series,fun)
{
  return(list(runner(series,k=120,f=fun)))
}


wti_roll_return<- rolling_sample_func(wti_df$wti,id_func)

test_stat_wti<- vector()

for(i in 1:length(wti_df$wti))
{
  test_stat_wti[i]<- test_func(wti_roll_return[[1]][[i]])
}

test_stat_wti<- test_stat_wti[-(1:119)]



wti_stat_df<- data.frame(wti_df$Date[-(1:119)],test_stat_wti)
colnames(wti_stat_df)<- c("Date","wti")


library(ggplot2)
wti_me_plot_120rw<- ggplot(data=wti_stat_df)+
  geom_line(aes(x=as.Date(Date),y=wti),linetype= "solid")+
  geom_hline(yintercept=crit_value_z2, linetype= "dashed")+
  xlab("Date")+ylab(NULL)+
  scale_y_continuous(limits = c(0, 50),
                     breaks = seq(0,50,5))+theme(axis.text = element_text(size = 10),axis.title=element_text(size=12,face="plain"))+scale_x_date(breaks = scales::breaks_pretty(8))+
  theme(plot.title = element_text(size = 15))+
  scale_x_date(breaks = scales::breaks_pretty(5))+
  ggtitle("WTI Crude Oil")

brent_me_plot_120rw<- brent_me_plot_120rw+scale_y_continuous(limits = c(0, 50),
                                                             breaks = seq(0,50,5))

drw120_plot<- grid.arrange(stock_me_plot_120rw,gold_me_plot_120rw,
                           btc_me_plot_120rw,gas_me_plot_120rw,
                           brent_me_plot_120rw,wti_me_plot_120rw,
                           nrow=3)
drw30_plot<- grid.arrange(oil_me_plot_30rw,btc_me_plot_30rw,btc_me_plot_30rw,gold_me_plot_30rw,gas_me_plot_30rw,wti_me_plot_30rw,nrow=3)

ggsave(filename = 'drw90_ee.pdf',
       plot = drw90_plot,
       width = 300,    
       height = 200,
       dpi = 1000,
       units = 'mm')
