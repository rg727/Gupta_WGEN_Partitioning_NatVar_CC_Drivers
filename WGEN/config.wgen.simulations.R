#Load in functions
setwd("/home/fs02/pmr82_0001/rg727/WGEN")
files.sources = list.files("./helpers",full.names = TRUE)
my.functions <- sapply(files.sources, source)


library(zoo)
library(mvtnorm)
library(MASS)
library(abind)



############################################# Load Data and Functions (Rohini) #####################################################
num.iter=1
basin.cnt=11
num.states=5

prct.prob <- 0.95 #GPD threshold

processed.data.filename.meteohydro <- paste0("/home/fs02/pmr82_0001/rg727/WGEN/Data/processed.data.files/processed.meteohydro/processed.meteohydro.",basin.cnt,".RData")
load(processed.data.filename.meteohydro) # processed meteohydro. data
n.sites <-380
ncol=380

prcp_sacsma=read.table("/home/fs02/pmr82_0001/rg727/WGEN/Final/TLG/TLG_Prcp.txt",sep="\t",header=TRUE)
tmin_sacsma=read.table("/home/fs02/pmr82_0001/rg727/WGEN/Final/TLG/TLG_Temp.txt",sep="\t",header=TRUE)
tmax_sacsma=read.table("/home/fs02/pmr82_0001/rg727/WGEN/Final/TLG/TLG_Temp.txt",sep="\t",header=TRUE)
n.sites=380
num.sites=380
prcp.basin.sacsma=rowMeans(prcp_sacsma)
tmin.basin.sacsma=rowMeans(tmin_sacsma)
tmax.basin.sacsma=rowMeans(tmax_sacsma)

#Cut days accordingly to remove warm season
start_date="1952-11-01"; end_date="2017-03-31"
first.date.weather <- as.Date(start_date)
last.date.weather <- as.Date(end_date)

months <- c(11,12,1,2,3,4)
first.month <- months[1]
last.month <- months[length(months)]
first.day <- 1
last.day <- 30

keep.dates.weather <- which(dates.weather>=first.date.weather & dates.weather<=last.date.weather & as.numeric(format(dates.weather,"%m"))%in%months)
leap.indices <- which(format(dates.weather,"%m")=="02" & format(dates.weather,"%d")=="29")

for (i in 1:length(keep.dates.weather)){
  if (keep.dates.weather[i]==leap.indices[1]){
    keep.dates.weather=keep.dates.weather[-c(i)]
  }
}


dates.weather <- dates.weather[keep.dates.weather]
prcp.basin <- prcp.basin.sacsma[keep.dates.weather]
tmax.basin <- tmax.basin.sacsma[keep.dates.weather]
tmin.basin <- tmin.basin.sacsma[keep.dates.weather]
prcp.site <- prcp_sacsma[keep.dates.weather,]
tmax.site <- tmax_sacsma[keep.dates.weather,]
tmin.site <- tmin_sacsma[keep.dates.weather,]

years.weather <- as.numeric(format(dates.weather,"%Y"))
months.weather <- as.numeric(format(dates.weather,"%m"))
days.weather <- as.numeric(format(dates.weather,"%d"))


#time series data for historical and simulated WRs
sNHMM <- readRDS("/home/fs02/pmr82_0001/rg727/WGEN/fit.modHMMs.depmix.paleo.april.rds")
dates.synoptic=readRDS("/home/fs02/pmr82_0001/rg727/WGEN/dates.truncated.rds")
sNHMM =sNHMM[1:11117] 
weather.state.assignments=sNHMM
dates.synoptic=dates.synoptic[1:11117]
#first.date.synoptic <- as.Date("1952-11-01")   
#last.date.synoptic <- as.Date("2017-3-31")
#dates.synoptic <- seq(first.date.synoptic,last.date.synoptic, by="days")  #this defines the date range of the weather regimes
#weather.state.assignments <- sNHMM$viterbi.seq[dates.synoptic%in%dates.weather]
#sim.weather.state.assignments <- sNHMM$matrix.hmms.seq.states[dates.synoptic%in%dates.weather,1:num.iter]
#rm(sNHMM) # for memory
#weather.state.assignments=sNHMM[304:10134]
#dates.synoptic=dates.synoptic[304:10134]
#months.weather=months.weather[304:10134]
#prcp.basin=prcp.basin[304:10134]
#prcp.site=prcp.site[304:10134,]
#tmax.site=tmax.site[304:10134,]
#tmin.site=tmin.site[304:10134,]
#Read in AR files

#AR <- read.table("./Data/simulated.data.files/AR_Presence_Absence_CA_only_only.txt",sep="")
AR <- read.table("/home/fs02/pmr82_0001/rg727/WGEN/Data/simulated.data.files/AR_Presence_Absence_categorical.txt",sep="")
AR$V1=as.character(AR$V1)
AR$V1=as.Date(AR$V1,"%Y-%m-%d")

AR_cut <- numeric(length(dates.synoptic))
for (i in 1:(length(dates.synoptic))) {
  AR_cut[i] <- which(AR[,1]==dates.synoptic[i])
}

predictor <- AR[AR_cut,2]



#The spearman correlation between basin and site precipitation, used in the copula-based jitters
prcp.basin.site <- cbind(prcp.basin,prcp.site)
#prcp.basin.site <- prcp.basin.site[1:(which(is.na(prcp.basin.site[,1]))[1]-1),]
S <- cor(prcp.basin.site,method="spearman")




#Quantile mapping choices
qq <- .999              # quantile to anchor the CC-scaling
perc.q.per.C <- 0.07    # (0.07) percent change (as decimal) in the qqth quantile for each month of non-zero prcp per degree C warming


#window size for bootstrapping 
window.size <- rep(20,length(months))   #the size of the window (in days) from which runs can be bootstrapped around the current day of simulation, specified separately by month



##############################Define perturbations#######################################
# ## For not HISTORICAL and with JITTER:
# change.list <- data.frame("hc"=  c(TRUE),      #do we use the historical WRs? TRUE if yes, FALSE indicates to use simulated WRs
#                           "tc"=  c( 3),          #the mean temperature change to apply
#                           "jc"=  c( TRUE),       #use the jitter algorithm? TRUE is yes, FALSE indicates to not jitter the bootstrapped data
#                           "pccc"=c( 2),          #values of 1, 2, 3 etc represent the multiplier on the C-C scaling rate on theoretical 7% (perc.q.per.C) rate
#                           "pmuc"=c( 0)           #ex. 0.03 would indicate mean precipitation thermodynamically scales at 3% per degree warming
# )
###############################################################################################



##########################Simulate model with perturbations#######################################

#currently the code is designed to run the weather generator over the same length of the historical record (can increase iter to lengthen the run)
df=readRDS("/home/fs02/pmr82_0001/rg727/WGEN/df_with_year_leap.rds")
dates.sim <- df$date

# #subset the change.list for those changes that require new runs of the weather generator (i.e., dynamic changes), rather than just different purturbations to the simulations (i.e., thermodynamic changes)
# hc.new.sim.list <- data.frame('hc'=change.list$hc)
# #find the unique changes for WR dynamics from the list above
# hc.new.sim.list.unique <- unique(hc.new.sim.list)

# #loop through the unique changes for WR dynamics in the list above
#   #the current, unique set of changes that require new simulations for WR dynamics
#   hc.new.sim=1
#   cur.hc <- hc.new.sim.list.unique[hc.new.sim,]
#   #the list of thermodynamic changes to apply in post-processing for this particular set of WR dynamics 
#   change.list.CC <- change.list[which(apply(hc.new.sim.list,1,function(x) {all.equal(as.numeric(x),as.numeric(cur.hc))})==TRUE),]
  
 
  ###############################################################################################################################################
  mc_simulations_replaced=readRDS("/home/fs02/pmr82_0001/rg727/WGEN/mcsimulations_with_year_leap.rds")
  prcp.site.sim.total <-lapply(1:50, matrix, data= NA, nrow=112008, ncol=380)
  tmax.site.sim.total <- lapply(1:50, matrix, data= NA, nrow=112008, ncol=380)
  tmin.site.sim.total <- lapply(1:50, matrix, data= NA, nrow=112008, ncol=380)
  resampled.date.sim.total=readRDS("/home/fs02/pmr82_0001/rg727/WGEN/Final/MIL/resampled.date.sim.total.rds")
  nsites=380
  
  for (i in 1:50){
    sim.weather.state.assignments <- mc_simulations_replaced[,i+1]
    markov.chain.sim <- as.list(data.frame(sim.weather.state.assignments))
    
    #run the daily weather generate iter times using the iter Markov chains
    #note: if this gets caught in an infinite loop, it occurs at line 94 in wgen.simulator, and suggests your simulated WRs are inconsistent with the observed WRs. the appropriate fix is to redo your simulated WRs
    my.sim.output <- run.simulations(weather.state.assignments=weather.state.assignments,markov.chain.sim=markov.chain.sim,
                                     predictor=predictor,prcp.basin=prcp.basin,dates.weather=dates.weather,
                                     first.month=first.month,last.month=last.month,num.iter=num.iter,dates.sim=dates.sim,
                                     months=months,window.size=window.size,
                                     prcp.site=prcp.site,tmax.site=tmax.site,tmin.site=tmin.site, resampled.dates=resampled.date.sim.total[[i]])    
    
    #each of these is a list of length iter
    mc.sim <- my.sim.output[[1]]
    resampled.date.sim <- my.sim.output[[2]]
    prcp.site.sim <- my.sim.output[[3]]
    tmax.site.sim <- my.sim.output[[4]]
    tmin.site.sim <- my.sim.output[[5]]
    
    prcp.site.sim.total[[i]]=prcp.site.sim
    tmax.site.sim.total[[i]]=tmax.site.sim
    tmin.site.sim.total[[i]]=tmin.site.sim
    
  }



#############################################################WARM SEASON BOOTSTRAP###################################################
#Create a warm season set of dates: 

dates_total=data.frame(dates.vec.all,as.numeric(format(dates.vec.all,"%m")))
colnames(dates_total)=c("Dates","month")
dates_warm=subset(dates_total, dates_total$month==5 |dates_total$month==6 |dates_total$month==7| dates_total$month==8 |dates_total$month==9|dates_total$month==10)
dates_warm$year=as.numeric(format(dates_warm$Dates,"%y"))
rownames(dates_warm) <- seq(length=nrow(dates_warm))

#Create a unique list of warm WRs 

unique_warm_WRs=unique(dates_warm$year)



all_season_precip=prcp.site
all_season_tmin=tmin.site
all_season_tmax=tmax.site

start_date="1950-01-01"; end_date="2013-12-31"
first.date.weather.warm <- as.Date(start_date)
last.date.weather.warm <- as.Date(end_date)

months.warm <- c(5,6,7,8,9,10)
first.month.warm <- months.warm[1]
last.month.warm <- months.warm[length(months.warm)]


keep.dates.weather.warm <- which(dates.vec.all>=first.date.weather.warm & dates.vec.all<=last.date.weather.warm & as.numeric(format(dates.vec.all,"%m"))%in%months.warm)
dates.weather.warm <- dates.vec.all[keep.dates.weather.warm]

prcp.site.warm <- prcp_sacsma[keep.dates.weather.warm,]
tmax.site.warm <- tmax_sacsma[keep.dates.weather.warm,]
tmin.site.warm <- tmin_sacsma[keep.dates.weather.warm,]


#Create total date sequence for paleo-period
all.dates=seq(as.Date("1399/11/1"),as.Date("2017/4/30"), by = "day")
start_date_paleo="1399-11-01"; end_date_paleo="2017/4/30"
first.date.weather.warm.paleo <- as.Date(start_date_paleo)
last.date.weather.warm.paleo <- as.Date(end_date_paleo)

keep.dates.weather.warm.paleo <- which(all.dates>=first.date.weather.warm.paleo & all.dates<=last.date.weather.warm.paleo & as.numeric(format(all.dates,"%m"))%in%months.warm)

all.dates.warm=all.dates[keep.dates.weather.warm.paleo]

months.cold <- c(11,12,1,2,3,4)
first.month.cold <- months.cold[1]
last.month.cold <- months.cold[length(months.cold)]

keep.dates.weather.cold.paleo <- which(all.dates>=first.date.weather.warm.paleo & all.dates<=last.date.weather.warm.paleo & as.numeric(format(all.dates,"%m"))%in%months.cold)

all.dates.cold=all.dates[keep.dates.weather.cold.paleo]


#Test how the boostrap will work 

#length(unique(years.sim))


iterations=50
#Create precipitation, min t, and max t lists
prcp.site.sim.warm <- lapply(1:50, matrix, data= NA, nrow=113528, ncol=ncol)

for (n in 1:iterations){
  counter=1
  warm.season.index=sample(1:64, 617, replace=TRUE)
  warm.season.samples=unique_warm_WRs[warm.season.index]
  for (i in 1:length(warm.season.samples)){
    example_index=which(dates_warm$year==warm.season.samples[i])
    test=data.matrix(prcp.site.warm[example_index,])
    prcp.site.sim.warm[[n]][counter:(counter+length(example_index)-1),]=test
    counter=counter+length(example_index)
  }
  
}




tmax.site.sim.warm <- lapply(1:50, matrix, data= NA, nrow=113528, ncol=ncol)

for (n in 1:iterations){
  counter=1
  warm.season.index=sample(1:64, 617, replace=TRUE)
  warm.season.samples=unique_warm_WRs[warm.season.index]
  for (i in 1:length(warm.season.samples)){
    
    example_index=which(dates_warm$year==warm.season.samples[i])
    tmax.site.sim.warm[[n]][counter:(counter+length(example_index)-1),]=data.matrix(tmax.site.warm[example_index,])
    counter=counter+length(example_index)
  }
  
}

tmin.site.sim.warm <- lapply(1:50, matrix, data= NA, nrow=113528, ncol=ncol)

for (n in 1:iterations){
  counter=1
  warm.season.index=sample(1:64, 617, replace=TRUE)
  warm.season.samples=unique_warm_WRs[warm.season.index]
  for (i in 1:length(warm.season.samples)){
    
    example_index=which(dates_warm$year==warm.season.samples[i])
    tmin.site.sim.warm[[n]][counter:(counter+length(example_index)-1),]=data.matrix(tmin.site.warm[example_index,])
    counter=counter+length(example_index)
  }
  
}

#Cut days accordingly to remove warm season
start_date="1952-11-01"; end_date="2017-03-31"
first.date.weather.cold <- as.Date(start_date)
last.date.weather.cold <- as.Date(end_date)

months.cold <- c(11,12,1,2,3,4)
first.month.cold <- months[1]
last.month.cold <- months[length(months)]
first.day.cold <- 1
last.day.cold <- 30

keep.dates.weather.cold <- which(dates.weather>=first.date.weather.cold & dates.weather<=last.date.weather.cold & as.numeric(format(dates.weather,"%m"))%in%months.cold)
dates.weather.cold <- dates.weather[keep.dates.weather.cold]
prcp.basin.cold <- prcp.basin[keep.dates.weather.cold]
tmax.basin.cold <- tmax.basin[keep.dates.weather.cold]
tmin.basin.cold <- tmin.basin[keep.dates.weather.cold]
prcp.site.cold <- prcp.site[keep.dates.weather.cold,]
tmax.site.cold <- tmax.site[keep.dates.weather.cold,]
tmin.site.cold <- tmin.site[keep.dates.weather.cold,]


###########Merge lists####################################
#dates.vec.all[1036] starts at 1952
prcp.total.sim.compiled <- lapply(1:50, matrix, data= NA, nrow=225536, ncol=ncol)
for (j in 1:50){
  
  #For example simulation 
  cold_season=data.frame(all.dates.cold,prcp.site.sim.total[[j]])
  colnames(cold_season)[1] <- "Date"
  warm_season=data.frame(all.dates.warm,prcp.site.sim.warm[[j]])
  colnames(warm_season)[1] <- "Date"
  
  names=colnames(cold_season)
  colnames(warm_season)=names
  new <- rbind(cold_season,warm_season)
  final=new[order(as.Date(new$Date, format="%m/%d/%Y")),]
  row.names(final) <- NULL
  
  #Remove the beginning dates
  
  dates=as.Date(final$Date, format="%m/%d/%Y")
  years.weather=as.numeric(format(dates,"%Y"))
  months.weather=as.numeric(format(dates,"%m"))
  days.weather=as.numeric(format(dates,"%d"))
  final=final[,2:(ncol+1)]
  prcp.total.sim.compiled[[j]]=final
  
  }

rm(prcp.site.sim.warm)
rm(prcp.site.sim.total)
saveRDS(prcp.total.sim.compiled,"/home/fs02/pmr82_0001/rg727/WGEN/Final/TLG/prcp.site.sim.total.whole.rds")
rm(prcp.total.sim.compiled)

tmax.total.sim.compiled <- lapply(1:50, matrix, data= NA, nrow=225536, ncol=ncol)
for (j in 1:50){
  
  #For example simulation 
  cold_season=data.frame(all.dates.cold,tmax.site.sim.total[[j]])
  colnames(cold_season)[1] <- "Date"
  warm_season=data.frame(all.dates.warm,tmax.site.sim.warm[[j]])
  colnames(warm_season)[1] <- "Date"
  
  names=colnames(cold_season)
  colnames(warm_season)=names
  new <- rbind(cold_season,warm_season)
  final=new[order(as.Date(new$Date, format="%m/%d/%Y")),]
  row.names(final) <- NULL
  
  #Remove the beginning dates
  
  dates=as.Date(final$Date, format="%m/%d/%Y")
  years.weather=as.numeric(format(dates,"%Y"))
  months.weather=as.numeric(format(dates,"%m"))
  days.weather=as.numeric(format(dates,"%d"))
  final=final[,2:(ncol+1)]
  tmax.total.sim.compiled[[j]]=final
  }

rm(tmax.site.sim.warm)
rm(tmax.site.sim.total)
saveRDS(tmax.total.sim.compiled,"/home/fs02/pmr82_0001/rg727/WGEN/Final/TLG/tmax.site.sim.total.whole.rds")
rm(tmax.total.sim.compiled)

tmin.total.sim.compiled <- lapply(1:50, matrix, data= NA, nrow=225536, ncol=ncol)
for (j in 1:50){
  
  #For example simulation 
  cold_season=data.frame(all.dates.cold,tmin.site.sim.total[[j]])
  colnames(cold_season)[1] <- "Date"
  warm_season=data.frame(all.dates.warm,tmin.site.sim.warm[[j]])
  colnames(warm_season)[1] <- "Date"
  
  names=colnames(cold_season)
  colnames(warm_season)=names
  new <- rbind(cold_season,warm_season)
  final=new[order(as.Date(new$Date, format="%m/%d/%Y")),]
  row.names(final) <- NULL
  
  #Remove the beginning dates
  
  dates=as.Date(final$Date, format="%m/%d/%Y")
  years.weather=as.numeric(format(dates,"%Y"))
  months.weather=as.numeric(format(dates,"%m"))
  days.weather=as.numeric(format(dates,"%d"))
  final=final[,2:(ncol+1)]
  tmin.total.sim.compiled[[j]]=final
  }

rm(tmin.site.sim.warm)
rm(tmin.site.sim.total)
saveRDS(tmin.total.sim.compiled,"/home/fs02/pmr82_0001/rg727/WGEN/Final/TLG/tmin.site.sim.total.whole.rds")
rm(tmin.total.sim.compiled)























  
  #saveRDS(prcp.site.sim.total,"/home/fs02/pmr82_0001/rg727/WGEN/prcp.site.sim.total.rds")
  #saveRDS(tmax.site.sim.total,"/home/fs02/pmr82_0001/rg727/WGEN/tmax.site.sim.total.rds")
  #saveRDS(tmin.site.sim.total,"/home/fs02/pmr82_0001/rg727/WGEN/tmin.site.sim.total.rds")
  
  # ####################################################################################################################################################  
  # # #remove for memory
  # # rm(my.sim.output)

  # prcp.site.sim.total=readRDS("/home/fs02/pmr82_0001/rg727/WGEN/prcp.site.sim.total.rds")
  # tmax.site.sim.total=readRDS("/home/fs02/pmr82_0001/rg727/WGEN/tmax.site.sim.total.rds")
  # tmin.site.sim.total=readRDS("/home/fs02/pmr82_0001/rg727/WGEN/tmin.site.sim.total.rds")
  
  # prcp.site.sim.total.perturbed <-lapply(1:50, matrix, data= NA, nrow=112008, ncol=380)
  # tmax.site.sim.total.perturbed <- lapply(1:50, matrix, data= NA, nrow=112008, ncol=380)
  # tmin.site.sim.total.perturbed <- lapply(1:50, matrix, data= NA, nrow=112008, ncol=380)
  
  
  # #once the simulations are created, we now apply post-process climate changes (and jitters)
 
  #   change=1
  #   cur.tc <- change.list.CC$tc[change]
  #   cur.jc <- change.list.CC$jc[change]
  #   cur.pccc <- change.list.CC$pccc[change]
  #   cur.pmuc <- change.list.CC$pmuc[change]
    
  #   #precipitation scaling (temperature change dependent)
  #   perc.mu <- (1 + cur.pmuc)^cur.tc       #scaling in the mean for each month of non-zero prcp
  #   perc.q <- (1 + perc.q.per.C*cur.pccc)^cur.tc    #CC scaling in the qqth quantile for each month of non-zero prcp
    
  #   for (i in 1:50){
      
  #     #perturb the climate from the simulations above (the longest procedure in this function is saving the output files)
  #     set.seed(1)   #this ensures the copula-based jitterrs are always performed in the same way for each climate change
  #     perturbed.sim <-  perturb.climate(prcp.site.sim.total[[i]],tmax.site.sim.total[[i]],tmin.site.sim.total[[i]],
  #                                       emission.fit.site=emission.fit.site,
  #                                       months=months,sim.weather.state.assignments=sim.weather.state.assignments,dates.sim=dates.sim,n.sites=n.sites,
  #                                       qq=qq,perc.mu=perc.mu,perc.q=perc.q,S=S,jitter=cur.jc,temp.change=cur.tc,num.iter=num.iter)
  #     set.seed(NULL)
  #     prcp.site.sim.total.perturbed[[i]] <- perturbed.sim[[1]]
  #     tmax.site.sim.total.perturbed[[i]] <- perturbed.sim[[2]]
  #     tmin.site.sim.total.perturbed[[i]] <- perturbed.sim[[3]]
      
    
    
  #   #remove for memory
  #   rm(perturbed.sim)
    
   
  #   }
    
  #   saveRDS(prcp.site.sim.total.perturbed,"/home/fs02/pmr82_0001/rg727/WGEN/prcp.site.sim.total.perturbed_3T_2CC_0M.rds")
  #   saveRDS(tmax.site.sim.total.perturbed,"/home/fs02/pmr82_0001/rg727/WGEN/tmax.site.sim.total.perturbed_3T_2CC_0M.rds")
  #   saveRDS(tmin.site.sim.total.perturbed,"/home/fs02/pmr82_0001/rg727/WGEN/tmin.site.sim.total.perturbed_3T_2CC_0M.rds")