library(extRemes)
library(zoo)

#Read in streamflow and dates

setwd("E:/WGEN_Paper_Figures/rg727_earthsfuture_WGEN_Decomposition/Metrics/Flood/")
streamflow=read.table("streamflow_0T_0CC.txt",sep="\t",header=TRUE)

dates=readRDS("dates.rds")

#############################################################10-Year Event (30-Year Window)##########################################################

#Aggregate to 3-day flow 
TLG_flows_rolled=rollsum(streamflow, 3, na.pad = TRUE, align = c("center"))

#Start with first ensemble member (actually is 7th column)
j=7

#Calculate peak annual flow

TLG_flow_peak_annual=aggregate(TLG_flows_rolled[,j],by=list(dates$water_year),FUN=max,na.rm=TRUE, na.action=NULL)

w=1

gevfit_TLG <- fevd(TLG_flow_peak_annual$x[w:(w+30)])
flow=return.level(gevfit_TLG, return.period = c(10))[1]
dataframe_0T_0CC=data.frame(flow)

#Create 30 -year moving window
for (w in 2:578){
  gevfit_TLG <- fevd(TLG_flow_peak_annual$x[w:(w+30)])
  flow=return.level(gevfit_TLG, return.period = c(10))[1]
  dataframe_0T_0CC=rbind(dataframe_0T_0CC,flow)
}


range=c(8:56)

for (j in range){
  #Calculate peak annual flow
  TLG_flow_peak_annual=aggregate(TLG_flows_rolled[,j],by=list(dates$water_year),FUN=max,na.rm=TRUE, na.action=NULL)
  
  #Create 30 -year moving window
  for (w in 1:578){
    gevfit_TLG <- fevd(TLG_flow_peak_annual$x[w:(w+30)])
    flow=return.level(gevfit_TLG, return.period = c(10))[1]
    dataframe_0T_0CC=rbind(dataframe_0T_0CC,flow)
    
  }
}

years=rep(c(1399:1976),50)
dataframe_0T_0CC_final=data.frame(years,dataframe_0T_0CC)
colnames(dataframe_0T_0CC_final)=c('years','flow')



write.csv(dataframe_0T_0CC_final,file = "Tuolumne_GEV_10yr_30yrMW.csv",row.names = FALSE)

######################################################100 Yr Event#############################################################
w=1

gevfit_TLG <- fevd(TLG_flow_peak_annual$x[w:(w+30)])
flow=return.level(gevfit_TLG, return.period = c(100))[1]
dataframe_0T_0CC=data.frame(flow)

#Create 30 -year moving window
for (w in 2:578){
  gevfit_TLG <- fevd(TLG_flow_peak_annual$x[w:(w+30)])
  flow=return.level(gevfit_TLG, return.period = c(100))[1]
  dataframe_0T_0CC=rbind(dataframe_0T_0CC,flow)
}


range=c(8:56)

for (j in range){
  #Calculate peak annual flow
  TLG_flow_peak_annual=aggregate(TLG_flows_rolled[,j],by=list(dates$water_year),FUN=max,na.rm=TRUE, na.action=NULL)
  
  #Create 30 -year moving window
  for (w in 1:578){
    gevfit_TLG <- fevd(TLG_flow_peak_annual$x[w:(w+30)])
    flow=return.level(gevfit_TLG, return.period = c(100))[1]
    dataframe_0T_0CC=rbind(dataframe_0T_0CC,flow)
    
  }
}

years=rep(c(1399:1976),50)
dataframe_0T_0CC_final=data.frame(years,dataframe_0T_0CC)
colnames(dataframe_0T_0CC_final)=c('years','flow')


write.csv(dataframe_0T_0CC_final,file = "Tuolumne_GEV_100yr_30yrMW.csv",row.names = FALSE)
