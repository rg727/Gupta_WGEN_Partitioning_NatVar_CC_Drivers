library(mvtnorm)
library(extRemes)
library(mnormt)
library(TLMoments)

#Read in annual observed flows 
TLG_flows=read.table("D:/copula/FNF_TLG_mm.txt")
MIL_flows=read.table("D:/copula/FNF_MIL_mm.txt")
NML_flows=read.table("D:/copula/FNF_NML_mm.txt")

#Calculate water year column
TLG_flows$water_year=0

for (i in 1:dim(TLG_flows)[1]){
  if (TLG_flows$V2[i]==10|TLG_flows$V2[i]==11|TLG_flows$V2[i]==12){
    TLG_flows$water_year[i]=TLG_flows$V1[i]+1
  } else{
    TLG_flows$water_year[i]=TLG_flows$V1[i]}
}


MIL_flows$water_year=0

for (i in 1:dim(MIL_flows)[1]){
  if (MIL_flows$V2[i]==10|MIL_flows$V2[i]==11|MIL_flows$V2[i]==12){
    MIL_flows$water_year[i]=MIL_flows$V1[i]+1
  } else{
    MIL_flows$water_year[i]=MIL_flows$V1[i]}
}


NML_flows$water_year=0

for (i in 1:dim(NML_flows)[1]){
  if (NML_flows$V2[i]==10|NML_flows$V2[i]==11|NML_flows$V2[i]==12){
    NML_flows$water_year[i]=NML_flows$V1[i]+1
  } else{
    NML_flows$water_year[i]=NML_flows$V1[i]}
}



#Calculate 3-day total flow


TLG_flows$three_day=0

for (i in 1:length(TLG_flows$V4)){
  if (i == 1){
    TLG_flows$three_day[i]=sum(TLG_flows$V4[i],TLG_flows$V4[i+1])
  }
  else
    TLG_flows$three_day[i]=sum(TLG_flows$V4[i-1],TLG_flows$V4[i],TLG_flows$V4[i+1])
}

MIL_flows$three_day=0

for (i in 1:length(MIL_flows$V4)){
  if (i == 1){
    MIL_flows$three_day[i]=sum(MIL_flows$V4[i],MIL_flows$V4[i+1])
  }
  else
    MIL_flows$three_day[i]=sum(MIL_flows$V4[i-1],MIL_flows$V4[i],MIL_flows$V4[i+1])
}

NML_flows$three_day=0

for (i in 1:length(NML_flows$V4)){
  if (i == 1){
    NML_flows$three_day[i]=sum(NML_flows$V4[i],NML_flows$V4[i+1])
  }
  else
    NML_flows$three_day[i]=sum(NML_flows$V4[i-1],NML_flows$V4[i],NML_flows$V4[i+1])
}



#Calculate peak annual flow
TLG_flow_peak_annual=aggregate(TLG_flows,by=list(TLG_flows$water_year),FUN=max,na.rm=TRUE, na.action=NULL)
MIL_flow_peak_annual=aggregate(MIL_flows,by=list(MIL_flows$water_year),FUN=max,na.rm=TRUE, na.action=NULL)
NML_flow_peak_annual=aggregate(NML_flows,by=list(NML_flows$water_year),FUN=max,na.rm=TRUE, na.action=NULL)

#Remove extra year (1986) from TLG
TLG_flow_peak_annual=TLG_flow_peak_annual[2:33,]



#Determine level of a 10-year return period over the historical period
gevfit_MIL <- fevd(MIL_flow_peak_annual$three_day)
ten_yr_event_MIL=return.level(gevfit_MIL, return.period = c(10))[1]

gevfit_NML <- fevd(NML_flow_peak_annual$three_day)
ten_yr_event_NML=return.level(gevfit_NML, return.period = c(10))[1]

gevfit_TLG <- fevd(TLG_flow_peak_annual$three_day)
ten_yr_event_TLG=return.level(gevfit_TLG, return.period = c(10))[1]



#Convert the thresholds to u's for each basin 
u_TLG_event = pgev(ten_yr_event_TLG, loc=gevfit_TLG$results$par[1], scale=gevfit_TLG$results$par[2], shape=gevfit_TLG$results$par[3])
u_MIL_event = pgev(ten_yr_event_MIL, loc=gevfit_MIL$results$par[1], scale=gevfit_MIL$results$par[2], shape=gevfit_MIL$results$par[3])
u_NML_event = pgev(ten_yr_event_NML, loc=gevfit_NML$results$par[1], scale=gevfit_NML$results$par[2], shape=gevfit_NML$results$par[3])


#Now takes these u's and push them through a standard normal to get Z scores 

Z_TLG_event=qnorm(u_TLG_event)
Z_MIL_event=qnorm(u_MIL_event)
Z_NML_event=qnorm(u_NML_event)


#Now we want to fit the copula on our historical flows. First create u's

u_TLG = pgev(TLG_flow_peak_annual$three_day, loc=gevfit_TLG$results$par[1], scale=gevfit_TLG$results$par[2], shape=gevfit_TLG$results$par[3])
u_MIL = pgev(MIL_flow_peak_annual$three_day, loc=gevfit_MIL$results$par[1], scale=gevfit_MIL$results$par[2], shape=gevfit_MIL$results$par[3])
u_NML = pgev(NML_flow_peak_annual$three_day, loc=gevfit_NML$results$par[1], scale=gevfit_NML$results$par[2], shape=gevfit_NML$results$par[3])

#Then create Z-scores 
Z_TLG=qnorm(u_TLG)
Z_MIL=qnorm(u_MIL)
Z_NML=qnorm(u_NML)

#Create a pearson correlation matrix to capture the correlation structure across all the basins
cor_matrix=cor(cbind(Z_TLG,Z_MIL,Z_NML),method = "spearman")


#Now calculate the historical probability of exceeding the threshold in all basins
exceedance=1-pnorm(Z_TLG_event)-pnorm(Z_MIL_event)-pnorm(Z_NML_event)+pmnorm(c(Z_TLG_event,Z_MIL_event),mean=rep(0,2),varcov = cor(cbind(Z_TLG,Z_MIL)))+pmnorm(c(Z_TLG_event,Z_NML_event),mean=rep(0,2),varcov = cor(cbind(Z_TLG,Z_NML)))+pmnorm(c(Z_MIL_event,Z_NML_event),mean=rep(0,2),varcov = cor(cbind(Z_MIL,Z_NML)))-pmnorm(c(Z_TLG_event,Z_MIL_event,Z_NML_event),mean=rep(0,3),varcov = cor(cbind(Z_TLG,Z_MIL,Z_NML)))

#exceedance probability is 0.0564


