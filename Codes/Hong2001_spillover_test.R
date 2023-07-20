library(vars)
library(mFilter)
library(tseries)
library(TSstudio)
library(forecast)
library(tidyverse)
library(wavelets)
library(rugarch)
library(rmgarch)
library(mgarchBEKK)
library(waveslim)
library(aTSA)
library(ggplot2)
library(tidyverse)
library(sjmisc)
library(runner)
library(ggplot2)
library(Hmisc)
library(boot)
library(blocklength)
library(AnalystHelper)
library(OBL)
library(doParallel)
library(parallel)
library(tidyverse)
library(tidyr)

setwd("C:/Users/Archi Roy/Desktop/Effects of Conflicts on Cryptocurrency Market Project/Code/All new things/Workspaces/Bootstrap subsetted")
#getting the return data
return_df<- as.data.frame(cbind(
  100*c(NA,diff(log(final_prices$Brent))),
  100*c(NA,diff(log(final_prices$Stock))),
  100*c(NA,diff(log(final_prices$Gold))),
  100*c(NA,diff(log(final_prices$BTC))),
  100*c(NA,diff(log(final_prices$Gas))),
  100*c(NA,diff(log(final_prices$WTI)))))
return_df<- cbind(as.Date(final_prices$Date),return_df)
#return_df$V1<- as.Date(return_df$V1)
colnames(return_df)<- c("Date","brent_return","stock_return","gold_return","btc_return","gas_return","wti_return")
return_df<- return_df[-1,]

return_df_dates<- seq(as.Date("2019-03-11"), as.Date("2022-12-30"), by="days")
return_df <- return_df %>%
  mutate(Date = as.Date(Date)) %>%
  complete(Date = seq(as.Date("2019-03-11"), as.Date("2022-12-30"), by="days") )

return_df$brent_return<- imputeTS::na.interpolation(return_df$brent_return,option = "stine")
return_df$stock_return<- imputeTS::na.interpolation(return_df$stock_return,option = "stine")
return_df$gold_return<- imputeTS::na.interpolation(return_df$gold_return,option = "stine")
return_df$btc_return<- imputeTS::na.interpolation(return_df$btc_return,option = "stine")
return_df$gas_return<- imputeTS::na.interpolation(return_df$gas_return,option = "stine")
return_df$wti_return<- imputeTS::na.interpolation(return_df$wti_return,option = "stine")


##############################################################################################
#Define some necessary functions
rolling_sample_func<- function(series,fun)
{
  return(list(runner(series,k=120,f=fun)))
}

significance_test_func<- function(series)
{
  sig_return <- ifelse(series > 1.645,"significant","insignificant")
  return(sig_return)
}

###############################################################################################
#This part is subject to change
s1wld<- mra(return_df$brent_return,wf="haar",J=7,method="modwt")
s1<- s1wld[["D1"]]

s2wld<- mra(return_df$stock_return,wf="haar",J=7,method="modwt")
s2<- s2wld[["D1"]]

y<- cbind(s1,s2)

rolling_sample1<- rolling_sample_func(y[,1],function(x){x})
rolling_sample2<- rolling_sample_func(y[,2],function(x){x})
##############################################################################################

#Define the prerequisites (a kernel of choice and M (the lag selection order))
library(changepoint)
bartlett_kernel<- function(z)
{
  if(abs(z)<1)
    return(1-abs(z))
  else
    return(0)
}
M=5
cit_junk<- vector()
dit_junk<- vector()
Q1part1<- vector()


Hong2001<- function(series1,series2)
{

#STEP 1: Check for structural breaks in variance of the input series
cpvar1<- cpt.var(series1)
cpvar1<- setdiff(cpvar1@cpts,c(1,1138))

cpvar2<- cpt.var(series2)
cpvar2<- setdiff(cpvar2@cpts,c(1,1138))

#STEP 2: Including the changepoints as dummy variables
cpn1<- length(cpvar1)
cpdates1 <- return_df$Date[cpvar1]
dummy1<-  matrix(0,nrow = length(series1),ncol = cpn1)
for (j in 1:cpn1){
  dummy1[(cpvar1[j]:length(series1)),j] = 1
}

cpn2<- length(cpvar2)
cpdates2 <- return_df$Date[cpvar2]
dummy2<-  matrix(0,nrow = length(series2),ncol = cpn2)
for (j in 1:cpn2){
  dummy2[(cpvar2[j]:length(series2)),j] = 1
}

#STEP 3: Fit a simple ARMA GARCH model of appropriate order including changepoints in variance equation
spec1<- ugarchspec(variance.model = list(model = "eGARCH",garchOrder = c(1,1),
                                      external.regressors = dummy1),
                mean.model = list(armaOrder = c(0,0,0),include.mean = T,
                                  external.regressors = NULL))
garch1<- ugarchfit(spec1,series1,solver="hybrid")
sqresid_standardized1<- (as.numeric(residuals(garch1,standardize=TRUE)))^2

spec2<- ugarchspec(variance.model = list(model = "eGARCH",garchOrder = c(1,1),
                                         external.regressors = dummy2),
                   mean.model = list(armaOrder = c(0,0,0),include.mean = T,
                                     external.regressors = NULL))
garch2<- ugarchfit(spec2,series2,solver="hybrid")
sqresid_standardized2<- (as.numeric(residuals(garch2,standardize=TRUE)))^2

#STEP 4: Compute the cross correlations between the two series of squared residuals upto lag M
rho<- ccf(sqresid_standardized1,sqresid_standardized2,lag.max = M,plot=FALSE)
rho<- rho[["acf"]][1:M]

len<- length(series1)
#STEP 5: Evaluate the other constants t<- length(return_df$Date)
for(j in 1:(len-1))
{
  cit_junk[j]<- (1-(j/len))*(bartlett_kernel(j/M)^2)
  dit_junk[j]<- (1-(j/len))*(1-((j+1)/len))*(bartlett_kernel(j/M)^4)
}
C1T<- sum(cit_junk)
D1T<- sum(dit_junk)


for(j in 1:(len-1))
{
 Q1part1[j]<- len*(bartlett_kernel(j/M)^2)*(rho[j]^2) 
}
Q1<- (sum(na.omit(Q1part1)-C1T))/sqrt(2*D1T)
return(c(Q1,significance_test_func(Q1)))
}

######################################################################################################
#Applying the test
hong2001results_forward<- array(numeric(),c(1,2,1391))
hong2001results_backward<- array(numeric(),c(1,2,1391))


for(i in 1:119)
{
  hong2001results_forward[,,i]<- matrix(NA,nrow=1,ncol=2)
  hong2001results_backward[,,i]<- matrix(NA,nrow=1,ncol=2)
}

for(a in 120:1391) {
  print(a)
  try({
    hong2001results_forward[ , ,a]<- Hong2001(unlist(rolling_sample1[[1]][[a]]),unlist(rolling_sample2[[1]][[a]]))
    hong2001results_backward[ , ,a]<- Hong2001(unlist(rolling_sample2[[1]][[a]]),unlist(rolling_sample1[[1]][[a]]))
  })
}


hong2001_stock_brent<- as.data.frame(na.omit(cbind(hong2001results_forward[1,1,],hong2001results_forward[1,2,])))
hong2001_stock_brent<- cbind(return_df$Date[-c(1:119)],hong2001_stock_brent)
colnames(hong2001_stock_brent)<- c("Date","stock-Brent-statistic","significance")

hong2001_brent_stock<- as.data.frame(na.omit(cbind(hong2001results_backward[1,1,],hong2001results_backward[1,2,])))
hong2001_brent_stock<- cbind(return_df$Date[-c(1:119)],hong2001_brent_stock)
colnames(hong2001_brent_stock)<- c("Date","Brent-stock-statistic","significance")


hong2001_gold_brent_plot<- ggplot(hong2001_stock_brent, aes(x=as.Date(Date), y=as.numeric(`stock-Brent-statistic`), color=as.factor(significance))) +
  geom_line(aes(group = 1))+
  scale_color_manual(values = c("darkgreen", "red")) + 
  ylab(NULL)+xlab("Date")+ theme(axis.text=element_text(size=5))+scale_x_date(breaks = scales::breaks_pretty(4))+theme(legend.position="none", axis.text = element_text(size = 10))+
  ggtitle("Gold->Brent (D1)")#filled circle is significant

hong2001_brent_gold_plot<- ggplot(hong2001_brent_stock, aes(x=as.Date(Date), y=as.numeric(`Brent-stock-statistic`), color=as.factor(significance))) +
  geom_line(aes(group = 1))+
  scale_color_manual(values = c("darkgreen", "red")) + 
  ylab(NULL)+xlab("Date")+ theme(axis.text=element_text(size=5))+scale_x_date(breaks = scales::breaks_pretty(4))+theme(legend.position="none", axis.text = element_text(size = 10))+
  ggtitle("Brent->Gold (D1)")#filled circle is significant

save(hong2001_gold_brent_plot,hong2001_brent_gold_plot,file="D1BG.Rdata")

rm(list = setdiff(ls(), "return_df"))