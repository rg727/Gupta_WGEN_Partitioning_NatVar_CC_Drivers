#Load in functions
setwd("/home/fs02/pmr82_0001/rg727/WGEN")
files.sources = list.files("./Programs/quick_2",full.names = TRUE)
my.functions <- sapply(files.sources, source)


library(zoo)
library(mvtnorm)
library(MASS)
library(abind)
library(arrow)
library(Matrix)
library(data.table)
library(eva)
library(fExtremes)
library(evmix)
library(rapportools)


args <- commandArgs()
C=as.numeric(args[6])




############################################# Load Data and Functions (Rohini) #####################################################
num.iter=50
basin.cnt=11
num.states=5



processed.data.filename.meteohydro <- paste0("/home/fs02/pmr82_0001/rg727/WGEN/Data/processed.data.files/processed.meteohydro/processed.meteohydro.",basin.cnt,".RData")
load(processed.data.filename.meteohydro) # processed meteohydro. data
n.sites <-380

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

months <- c(11,12,1,2,3,4,5,6,7,8,9,10)
first.month <- months[1]
last.month <- months[length(months)]
first.day <- 1
last.day <- 31

keep.dates.weather <- which(dates.weather>=first.date.weather & dates.weather<=last.date.weather & as.numeric(format(dates.weather,"%m"))%in%months)
leap.indices <- which(format(dates.weather,"%m")=="02" & format(dates.weather,"%d")=="29")

for (i in 1:length(keep.dates.weather)){
  if (keep.dates.weather[i]==leap.indices[1]){
    keep.dates.weather=keep.dates.weather[-c(i)]
  }
}


dates.weather <- dates.vec.all
prcp.basin <- prcp.basin.sacsma
tmax.basin <- tmax.basin.sacsma
tmin.basin <- tmin.basin.sacsma
prcp.site <- prcp_sacsma
tmax.site <- tmax_sacsma
tmin.site <- tmin_sacsma



#years.weather <- as.numeric(format(dates.vec.all,"%Y"))
months.weather <- as.numeric(format(dates.vec.all,"%m"))
#days.weather <- as.numeric(format(dates.vec.all,"%d"))


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
#AR <- read.table("/home/fs02/pmr82_0001/rg727/WGEN/Data/simulated.data.files/AR_Presence_Absence_categorical.txt",sep="")
#AR$V1=as.character(AR$V1)
#AR$V1=as.Date(AR$V1,"%Y-%m-%d")

#AR_cut <- numeric(length(dates.synoptic))
#for (i in 1:(length(dates.synoptic))) {
 # AR_cut[i] <- which(AR[,1]==dates.synoptic[i])
#}

#predictor <- AR[AR_cut,2]



#The spearman correlation between basin and site precipitation, used in the copula-based jitters
prcp.basin.site <- cbind(prcp.basin,prcp.site)
#prcp.basin.site <- prcp.basin.site[1:(which(is.na(prcp.basin.site[,1]))[1]-1),]
S <- cor(prcp.basin.site,method="spearman")


#fit emission distributions to each site by weather regimes and month (currently hard-coded for gamma distribution). sites along the columns, parameters for each weather regime and month down the rows

#Quantile mapping choices
perc.q.per.C <- 0.07    # (0.07) percent change (as decimal) in the qqth quantile for each month of non-zero prcp per degree C warming

qq <- .99              # percentile threshold to separate Gamma and GPD distributions
thshd.prcp <- apply(prcp.site,2,function(x) {quantile(x[x!=0],qq,na.rm=T)})


#Bootstrapping choices###
window.size <- rep(10,length(months))   #the size of the window (in days) from which runs can be bootstrapped around the current day of simulation, by month


emission.fits.site <- fit.emission(prcp.site=prcp.site,
                                   months=months,
                                   months.weather=months.weather,
                                   n.sites=n.sites,
                                   thshd.prcp=thshd.prcp)

#how often is prcp under threshold by month and site
qq.month <- sapply(1:n.sites,function(i,x=prcp.site,m,mm) {
  sapply(m,function(m) {
    length(which(as.numeric(x[x[,i]!=0 & mm==m & !is.na(x[,i]),i])<=thshd.prcp[i]))/length(as.numeric(x[x[,i]!=0 & mm==m & !is.na(x[,i]),i]))
  })},
  m=months, mm=months.weather)

##############################Define perturbations#######################################
## For not HISTORICAL and with JITTER:
change.list <- data.frame("hc"=  c(TRUE),      #do we use the historical WRs? TRUE if yes, FALSE indicates to use simulated WRs
                          "tc"=  c( 1),          #the mean temperature change to apply
                          "jc"=  c( TRUE),       #use the jitter algorithm? TRUE is yes, FALSE indicates to not jitter the bootstrapped data
                          "pccc"=c( C),          #values of 1, 2, 3 etc represent the multiplier on the C-C scaling rate on theoretical 7% (perc.q.per.C) rate
                          "pmuc"=c( 0)           #ex. 0.03 would indicate mean precipitation thermodynamically scales at 3% per degree warming
)
###############################################################################################



##########################Simulate model with perturbations#######################################

#currently the code is designed to run the weather generator over the same length of the historical record (can increase iter to lengthen the run)
dates.sim=readRDS("/home/fs02/pmr82_0001/rg727/WGEN/dates.sim.rds")


#subset the change.list for those changes that require new runs of the weather generator (i.e., dynamic changes), rather than just different purturbations to the simulations (i.e., thermodynamic changes)
hc.new.sim.list <- data.frame('hc'=change.list$hc)
#find the unique changes for WR dynamics from the list above
hc.new.sim.list.unique <- unique(hc.new.sim.list)

#loop through the unique changes for WR dynamics in the list above
  #the current, unique set of changes that require new simulations for WR dynamics
  hc.new.sim=1
  cur.hc <- hc.new.sim.list.unique[hc.new.sim,]
  #the list of thermodynamic changes to apply in post-processing for this particular set of WR dynamics 
  change.list.CC <- change.list[which(apply(hc.new.sim.list,1,function(x) {all.equal(as.numeric(x),as.numeric(cur.hc))})==TRUE),]
  
 
  ###############################################################################################################################################
  #mc_simulations_replaced=readRDS("/home/fs02/pmr82_0001/rg727/WGEN/mcsimulations_with_year_leap.rds")
  # prcp.site.sim.total <-lapply(1:50, matrix, data= NA, nrow=112008, ncol=380)
  # tmax.site.sim.total <- lapply(1:50, matrix, data= NA, nrow=112008, ncol=380)
  # tmin.site.sim.total <- lapply(1:50, matrix, data= NA, nrow=112008, ncol=380)
  # nsites=380
  
  # for (i in 1:50){
  #   sim.weather.state.assignments <- mc_simulations_replaced[,i+1]
  #   markov.chain.sim <- as.list(data.frame(sim.weather.state.assignments))
    
  #   #run the daily weather generate iter times using the iter Markov chains
  #   #note: if this gets caught in an infinite loop, it occurs at line 94 in wgen.simulator, and suggests your simulated WRs are inconsistent with the observed WRs. the appropriate fix is to redo your simulated WRs
  #   my.sim.output <- run.simulations(weather.state.assignments=weather.state.assignments,markov.chain.sim=markov.chain.sim,
  #                                    predictor=predictor,prcp.basin=prcp.basin,dates.weather=dates.weather,
  #                                    first.month=first.month,last.month=last.month,num.iter=num.iter,dates.sim=dates.sim,
  #                                    months=months,window.size=window.size,
  #                                    prcp.site=prcp.site,tmax.site=tmax.site,tmin.site=tmin.site)    
    
  #   #each of these is a list of length iter
  #   mc.sim <- my.sim.output[[1]]
  #   resampled.date.sim <- my.sim.output[[2]]
  #   prcp.site.sim <- my.sim.output[[3]]
  #   tmax.site.sim <- my.sim.output[[4]]
  #   tmin.site.sim <- my.sim.output[[5]]
    
  #   prcp.site.sim.total[[i]]=prcp.site.sim
  #   tmax.site.sim.total[[i]]=tmax.site.sim
  #   tmin.site.sim.total[[i]]=tmin.site.sim
    
  # }
  
  # saveRDS(prcp.site.sim.total,"/home/fs02/pmr82_0001/rg727/WGEN/prcp.site.sim.total.rds")
  # saveRDS(tmax.site.sim.total,"/home/fs02/pmr82_0001/rg727/WGEN/tmax.site.sim.total.rds")
  # saveRDS(tmin.site.sim.total,"/home/fs02/pmr82_0001/rg727/WGEN/tmin.site.sim.total.rds")
  
  # ####################################################################################################################################################  
  # #remove for memory
  # rm(my.sim.output)
  
  lat=read.table("/home/fs02/pmr82_0001/rg727/WGEN/Final/TLG/lat.txt",header=TRUE)
  lon=read.table("/home/fs02/pmr82_0001/rg727/WGEN/Final/TLG/lon.txt",header=TRUE)

  prcp.site.sim.total=readRDS("/home/fs02/pmr82_0001/rg727/WGEN/Final/TLG/prcp.site.sim.total.whole.rds")
  #tmax.site.sim.total=readRDS("/home/fs02/pmr82_0001/rg727/WGEN/Final/TLG/tmax.site.sim.total.whole.rds")
  tmin.site.sim.total=readRDS("/home/fs02/pmr82_0001/rg727/WGEN/Final/TLG/tmin.site.sim.total.whole.rds")
  
  #prcp.site.sim.total.perturbed <-lapply(1:50, matrix, data= NA, nrow=225536, ncol=380)
  #tmax.site.sim.total.perturbed <- lapply(1:50, matrix, data= NA, nrow=225536, ncol=380)
  #tmin.site.sim.total.perturbed <- lapply(1:50, matrix, data= NA, nrow=225536, ncol=380)
  
  n.sites=380
  num.sites=380
  #once the simulations are created, we now apply post-process climate changes (and jitters)
 
    change=1
    cur.tc <- change.list.CC$tc[change]
    cur.jitter <- change.list.CC$jc[change]
    cur.pccc <- change.list.CC$pccc[change]
    cur.pmuc <- change.list.CC$pmuc[change]
    
    #precipitation scaling (temperature change dependent)
    perc.mu <- (1 + cur.pmuc)^cur.tc       #scaling in the mean for each month of non-zero prcp
    perc.q <- (1 + perc.q.per.C*cur.pccc)^cur.tc    #CC scaling in the qqth quantile for each month of non-zero prcp
    

    
    # for (i in 1:50){
      
    #   #perturb the climate from the simulations above (the longest procedure in this function is saving the output files)
    #   set.seed(1)   #this ensures the copula-based jitterrs are always performed in the same way for each climate change
    #   perturbed.sim <-  perturb.climate(prcp.site.sim.total,tmax.site.sim.total,tmin.site.sim.total,
    #                                     emission.fit.site=emission.fit.site,
    #                                     months=months,sim.weather.state.assignments=sim.weather.state.assignments,dates.sim=dates.sim,n.sites=n.sites,
    #                                     qq=qq,perc.mu=perc.mu,perc.q=perc.q,S=S,jitter=cur.jc,temp.change=cur.tc,num.iter=num.iter)
    #   set.seed(NULL)
    #   prcp.site.sim.total.perturbed[[i]] <- perturbed.sim[[1]]
    #   tmax.site.sim.total.perturbed[[i]] <- perturbed.sim[[2]]
    #   tmin.site.sim.total.perturbed[[i]] <- perturbed.sim[[3]]
      
    
    
    # #remove for memory
    # rm(perturbed.sim)
    
   
    # }

      num.iter=50
      #perturb the climate from the simulations above (the longest procedure in this function is saving the output files)
      set.seed(1)   #this ensures the copula-based jitterrs are always performed in the same way for each climate change
      
      perturbed.sim <-  perturb.climate(prcp.site.sim=prcp.site.sim.total,
                                    tmin.site.sim=tmin.site.sim.total,
                                    emission.fits.site=emission.fits.site,
                                    months=months,dates.sim=dates.sim,n.sites=n.sites,
                                    qq=qq,perc.mu=perc.mu,perc.q=perc.q,S=S,cur.jitter=cur.jitter,cur.tc=cur.tc,
                                    num.iter=num.iter,thshd.prcp=thshd.prcp,qq.month=qq.month)
                                    
                                    
   
      set.seed(NULL)
      
      
      
      #rm(prcp.site.sim.total)
      #rm(tmax.site.sim.total)
      #rm(tmin.site.sim.total)
      prcp.site.sim.total.perturbed <- perturbed.sim[[1]]
      #tmax.site.sim.total.perturbed <- perturbed.sim[[2]]
      tmin.site.sim.total.perturbed <- perturbed.sim[[2]]
      
    
    
    #remove for memory
    rm(perturbed.sim)
    
   
  

for (j in 29:50){
  
    final=prcp.site.sim.total.perturbed[[j]]
    for (i in 1:380){
      test=data.frame(final[ ,i])
      write_parquet(test,paste0("/home/fs02/pmr82_0001/rg727/WGEN/Final/TLG/perturbations_2.0/1T_",C,"CC/precip/",j,"/meteo_",sprintf("%1.6f",lat[i,]),"_",sprintf("%1.6f",lon[i,]),".parquet",compression = "gzip"))
    }
}


for (j in 29:50){
  
  
    final=tmin.site.sim.total.perturbed[[j]]
    for (i in 1:380){
      test=data.frame(final[ ,i])
      write_parquet(test,paste0("/home/fs02/pmr82_0001/rg727/WGEN/Final/TLG/perturbations_2.0/1T_",C,"CC/temp/",j,"/meteo_",sprintf("%1.6f",lat[i,]),"_",sprintf("%1.6f",lon[i,]),".parquet",compression = "gzip"))
    }
}