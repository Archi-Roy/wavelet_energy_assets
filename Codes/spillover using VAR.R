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
library(Spillover)



setwd("C:/Users/Archi Roy/Desktop/Effects of Conflicts on Cryptocurrency Market Project/Code/All new things/Workspaces/Bootstrap subsetted")

#getting the return data
return_df<- as.data.frame(cbind(
  100*c(NA,diff(log(final_data$oil_price))),
  100*c(NA,diff(log(final_data$stock_price))),
  100*c(NA,diff(log(final_data$gold_price))),
  100*c(NA,diff(log(final_data$btc_price))),
  100*c(NA,diff(log(final_data$gas_price))),
  100*c(NA,diff(log(final_data$wti_price)))))
return_df<- cbind(as.Date(final_data$Date),return_df)
#return_df$V1<- as.Date(return_df$V1)
colnames(return_df)<- c("Date","brent_return","stock_return","gold_return","btc_return","gas_return","wti_return")
return_df<- return_df[-1,]

return_df_dates<- seq(as.Date("2019-08-02"), as.Date("2022-09-12"), by="days")
return_df <- return_df %>%
  mutate(Date = as.Date(Date)) %>%
  complete(Date = seq(as.Date("2019-08-02"), as.Date("2022-09-12"), by="days") )

return_df$brent_return<- imputeTS::na.interpolation(return_df$brent_return,option = "stine")
return_df$stock_return<- imputeTS::na.interpolation(return_df$stock_return,option = "stine")
return_df$gold_return<- imputeTS::na.interpolation(return_df$gold_return,option = "stine")
return_df$btc_return<- imputeTS::na.interpolation(return_df$btc_return,option = "stine")
return_df$gas_return<- imputeTS::na.interpolation(return_df$gas_return,option = "stine")
return_df$wti_return<- imputeTS::na.interpolation(return_df$wti_return,option = "stine")

#analyzing spillover using https://doi.org/10.1016/j.ijforecast.2011.02.006
#we can also try replicating https://rpubs.com/rsayed/573439

#we would first divide the data in 120 days (4 months) subsets
brent_decomposed<- mra(return_df$brent_return,wf="haar",J=6,method="modwt")
stock_decomposed<- mra(return_df$stock_return,wf="haar",J=6,method="modwt")
gold_decomposed<- mra(return_df$gold_return,wf="haar",J=6,method="modwt")
btc_decomposed<- mra(return_df$btc_return,wf="haar",J=6,method="modwt")
gas_decomposed<- mra(return_df$gas_return,wf="haar",J=6,method="modwt")
wti_decomposed<- mra(return_df$wti_return,wf="haar",J=6,method="modwt")

################################################ #this part is subject to change
brent_decomposed_scaled<- brent_decomposed[["D6"]]
stock_decomposed_scaled<- stock_decomposed[["D6"]]
gold_decomposed_scaled<-  gold_decomposed[["D6"]]
btc_decomposed_scaled<- btc_decomposed[["D6"]]
gas_decomposed_scaled<- gas_decomposed[["D6"]]
wti_decomposed_scaled<- wti_decomposed[["D6"]]
#################################################################

data_scaled<- cbind(brent_decomposed_scaled,stock_decomposed_scaled,
                  gold_decomposed_scaled,btc_decomposed_scaled,
                  gas_decomposed_scaled,wti_decomposed_scaled)

library(reshape2)

n <- 120
nr <- nrow(data_scaled)
test_split<- split(as.data.frame(data_scaled), rep(1:ceiling(nr/n), each=n, length.out=nr))

#the last dataframe in the list contains 58 observations for each series, the G.spillover function does not work for that
#These observations date from July 17,2022- September 12,2022
#For the sake of convenience I am eliminating this part, since our objective is to see the impact of the war, data from July to September, 2022 might not be useful
test_split[length(test_split)]<- NULL

Gspillfunction<- function(data) #input= test_spilt[[i]]
{
  #Step 1:Select the optimal VAR order using AIC criteria
  varorder<- VARselect(data,type="const")$selection[1]
  #Step 2:Fitting the VAR model using the above order
  varmodel<- VAR(data,p=varorder,type="const")
  #Step 3:Getting the generalized spillover index
  return(G.spillover(varmodel,n.ahead = 10,standardized = TRUE))
}

D1_spilltable_list<- lapply(test_split, Gspillfunction)

tidyspilltable<- function(tab)
{
  colnames(tab)<- c("Brent","S&P500","Gold","BTC","NG","WTI","From all others")
  rownames(tab)<- c("Brent","S&P500","Gold","BTC","NG","WTI","To all others","To all others including self")
  return(tab)
}
D1_spilltable_list<- lapply(D1_spilltable_list,tidyspilltable)


heatmap_func<- function(matrix)
{
  matrix<- matrix[-c(7,8),-c(7,8)]
  melted_matrix<- melt(matrix)
  library(ggplot2)
  return(ggplot(data = melted_matrix, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile())
}

D1_spill_heatmap<- lapply(D1_spilltable_list,heatmap_func)

library(grid)
library(gridExtra)

final_plot<- grid.arrange(arrangeGrob(D1_spill_heatmap[[1]],top="August,2019- November,2019"),
                          arrangeGrob(D1_spill_heatmap[[2]],top="November,2019- March,2020"),
                          arrangeGrob(D1_spill_heatmap[[3]],top = "March,2020- July,2020"),
                          arrangeGrob(D1_spill_heatmap[[4]],top="July,2020- November,2020"),
                          arrangeGrob(D1_spill_heatmap[[5]],top="November,2020- March,2021"),
                          arrangeGrob(D1_spill_heatmap[[6]],top="March,2021- July,2021"),
                          arrangeGrob(D1_spill_heatmap[[7]],top = "July,2021- November,2021"),
                          arrangeGrob(D1_spill_heatmap[[8]],top="November,2021- March,2022"),
                          arrangeGrob(D1_spill_heatmap[[9]],top = "March,2022- July,2022"),
                          nrow=3,ncol=3)

total_dynamic <- total.dynamic.spillover(as.zoo(data_scaled),
                                         width = 120, 
                                         index="generalized",
                                         p=4)
require(ggplot2)

d6<- tibble(Date= return_df$Date[-c(1:119)],  index = total_dynamic) %>% 
  mutate(Date=as.Date(as.character(Date))) %>% 
  ggplot(aes(x=Date, y=index)) +
  geom_line()+
  labs(caption = "Total volatility spillovers")+
  theme(plot.caption = element_text(hjust = 0))

rm(list = ls()[!ls() %in% c("return_df","a","d2","d3","d4","d5","d6")])

total_spill_plot<- grid.arrange(arrangeGrob(a,top="D1"),
                                arrangeGrob(d2,top="D2"),
                                arrangeGrob(d3,top="D3"),
                                arrangeGrob(d4,top="D4"),
                                arrangeGrob(d5,top="D5"),
                                arrangeGrob(d6,top="D6"),
                                nrow=2)