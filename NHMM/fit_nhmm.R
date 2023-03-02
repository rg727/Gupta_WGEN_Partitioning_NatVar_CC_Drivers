library(depmixS4)
library(ncdf4)
library(DirichletReg)
library(maps)
library(parallel)
library(markovchain)


############################Fit NHMM with paleo-reconstructed instrumental covariates#################################
#First grab the 500 hpa geopotential height and cut down to the right dates

setwd("E:/Box Sync/Weather_Generator/Data/processed.data.files")

processed.data.filename.atmosclim <- "./processed.data.atmosclim.RData"
load(processed.data.filename.atmosclim) # processed atmosclim. data

hgt500.synoptic.region.pca <- prcomp(hgt500.synoptic.region,center=T,scale=T)
num_eofs <- 9 # number of PCs
synoptic.pcs <- hgt500.synoptic.region.pca$x[,1:num_eofs]

synoptic.pcs.truncated=synoptic.pcs[304:dim(synoptic.pcs)[1],]

###############################Gupta 2022 Method for Getting PC_WR#################################################
#Next we calculate the PC_WR as in Gupta et al., 2022

setwd("E:/DriveF_Backup/NHMM/WRRecon_2021")

source("DirichReg_optim.R")
source("DirichReg_predict_cv.R")
source("DirichReg_reconstruct.R")

param=c(11,0.013415209)

use.smooth=TRUE 
num.pcs <- round(param[1])        #number of SPI PCs to use in reconstruction
cur.span <- param[2]      #key parameter for smoothing of SPI. set to enable a 10-year smoothing



my.nc <- nc_open("SPI_recon_NDJFM_0.5deg_Master_Extended_v2.nc")
print(my.nc)
spi_recon_org <- ncvar_get(my.nc,"spi_recon")
yr_recon <- ncvar_get(my.nc,"yr_recon")
vrsq_org <- ncvar_get(my.nc,"vrsq")
lat <- ncvar_get(my.nc,"lat")
lon <- ncvar_get(my.nc,"lon")
nc_close(my.nc)
lon_lat_org <- expand.grid(lon,lat)


#process SPI data into a matrix
spi_recon <- aperm(spi_recon_org,c(3,1,2))
dim(spi_recon) <- c(dim(spi_recon)[1],prod(dim(spi_recon)[2:3]))
spi_recon <- apply(spi_recon,2,function(x){scale(x)[,1]})
vrsq <- aperm(vrsq_org,c(3,1,2))
dim(vrsq) <- c(dim(vrsq)[1],prod(dim(vrsq)[2:3]))

#drop any grid cells that don't go back to 1400
keep <- apply(spi_recon,c(2),function(x) {length(which(is.na(x)))==0})
spi_recon <- spi_recon[,keep]
lon_lat <- lon_lat_org[keep,]
vrsq <- vrsq[,keep]

#PCs of SPI data
spi_recon.pca <- prcomp(spi_recon,center=TRUE,scale=TRUE)
spi_recon.pcs <- spi_recon.pca$x
colnames(spi_recon.pcs) <- paste("PC",1:ncol(spi_recon.pcs),sep="")

#######################################################################################################
num.states=5
synoptic.state.assignments <- read.table("synoptic_state_assignments.txt")
yr <- unique(synoptic.state.assignments$V2)
n.year <- length(yr)
WR <- array(NA,c(n.year,num.states))
for (i in 1:n.year) {
  WR[i,] <- tabulate(synoptic.state.assignments$V3[synoptic.state.assignments$V2==yr[i]],nbins=num.states)
}

#define WR fractions
WR_frac <- WR / apply(WR,FUN=sum,1)
WR_frac.DR <- DR_data(WR_frac)

#PCA on WRs
mypc <- prcomp(WR,center=TRUE,scale=TRUE)
(mypc$sdev^2)/sum(mypc$sdev^2)
WR.PCs <- mypc$x

WR_frac.df <- data.frame(WR.PCs,WR_frac)
WR_frac.df$Y <- DR_data(WR_frac.df[,grep("X",names(WR_frac.df))])
WR_frac.df <- WR_frac.df[,-grep("X",names(WR_frac.df))]
my.DirichReg <- DirichReg(Y ~ PC1 + PC2 + PC3 + PC4, data=WR_frac.df)   #currently hardcoded for 4 WR PCs

#smooth PSI
if(use.smooth) {
  spi_recon.pcs.smoothed <- apply(spi_recon.pcs,2,function(x){
    loess(x~yr_recon,span=cur.span)$fitted
  })
} else {
  spi_recon.pcs.smoothed <- spi_recon.pcs
}


#instrumental, smoothed SPI PCs for fitting reconstruction
spi.pcs.instr <- spi_recon.pcs.smoothed[yr_recon%in%yr,1:num.pcs]
#full smoothed SPI PCs for fitting reconstruction
spi.pcs.recon <- spi_recon.pcs.smoothed[,1:num.pcs]

#predict WR PCs based on SPI PCs
n.year.recon=621
pred.pcs <- array(NA,c(n.year.recon,num.states))
my.lm <- list(num.states)
for (i in 1:num.states) {
  cur.var <- WR.PCs[,i]
  fit.data <- data.frame('cur.var'=cur.var,'x'=spi.pcs.instr)
  my.lm[[i]] <- lm(cur.var~.,data=fit.data)
  pred.pcs[,i] <- predict(my.lm[[i]],newdata=data.frame('x'=spi.pcs.recon))    
}

#########################################Fit NHMM#########################################################################

#Keep paleo-reconstructed covariates during the instrumental period because those are the covariates that we will condition the NHMM on  

setwd("E:/DriveF_Backup/NHMM")
years.synoptic=readRDS("years.synoptic.rds")
nino34.annual=readRDS("nino34.annual.rds")

my.covariate.df=data.frame(yr,pred.pcs[554:618,1:4])

#Remove 1951 and 1952 from the record and create new synoptic states and years time series 
synoptic.state.assignments <- readRDS("HMMs_5states.rds")
nino34.annual_truncated=nino34.annual[3:67,]
years.synoptic.truncated=years.synoptic[304:10134]


#Create a daily covariate (just repeats the covariates for all days of the specific year)
my.match <- numeric(length(years.synoptic.truncated))
for (i in 1:length(years.synoptic.truncated)) {
  my.match[i] <- which(nino34.annual_truncated[ ,1]==years.synoptic.truncated[i])
}
predictor <- my.covariate.df[my.match,2:5]



modHMMs <- depmix(list(PC1~1,PC2~1,PC3~1, PC4~1, PC5~1, PC6~1, PC7~1, PC8~1,
                       PC9~1),
                  nstates = num.states,
                  family=list(gaussian(),gaussian(),gaussian(),gaussian(),gaussian(),gaussian(),gaussian(),gaussian(),
                              gaussian()),
                  ntimes =  nrow(synoptic.pcs.truncated),
                  data = data.frame(synoptic.pcs.truncated),transition = ~predictor$X1+predictor$X2+predictor$X3+predictor$X4)
fit.modHMMs.depmix.paleo <- fit(modHMMs)


synoptic.state.assignments.paleo.april <- posterior(fit.modHMMs.depmix.paleo)$state # state sequence (using the Viterbi algorithm)


#check that the states match the current way we've been defining WRs (that they aren't mixed up).
#Then plot the WR frequencies over the instrumental period. 

reordered_states_paleo_april=numeric(length(synoptic.state.assignments.paleo.april))


#This is an example reordering

for (i in 1:length(synoptic.state.assignments.paleo.april)){
  if (synoptic.state.assignments.paleo.april[i]==1){
    reordered_states_paleo_april[i]=4}
  if (synoptic.state.assignments.paleo.april[i]==2){
    reordered_states_paleo_april[i]=1}
  if (synoptic.state.assignments.paleo.april[i]==3){
    reordered_states_paleo_april[i]=2}
  if (synoptic.state.assignments.paleo.april[i]==4){
    reordered_states_paleo_april[i]=5}
  if (synoptic.state.assignments.paleo.april[i]==5){
    reordered_states_paleo_april[i]=3}
  
}


num.states <- 5

#test <- read.table("F:/NHMM/synoptic_state_assignments.txt")
#yr.test <- unique(test$V2)
n.year <- length(yr)
WR <- array(NA,c(n.year,num.states))
paleo.informed.states=data.frame(years.synoptic.truncated,reordered_states_paleo_april)
#synoptic.state.assignments$V3=reordered_states
for (i in 1:n.year) {
  WR[i,] <- tabulate(paleo.informed.states$reordered_states_paleo_april[paleo.informed.states$years.synoptic.truncated==yr[i]],nbins=num.states)
}

#define WR fractions
WR_frac_NHMM <- WR / apply(WR,FUN=sum,1)
WR_frac_NHMM.DR <- DR_data(WR_frac_NHMM)


#################################################Create Daily Matrices################################################################################################
#Set up transition matrix formula 
transition_matrix=matrix(, nrow = 5, ncol = 5)

#Row 1 = State 1 
for (i in 6:10){
  print(getpars(fit.modHMMs.depmix.paleo)[i])
  print(getpars(fit.modHMMs.depmix.paleo)[i+5])
  print(getpars(fit.modHMMs.depmix.paleo)[i+10])
  print(getpars(fit.modHMMs.depmix.paleo)[i+15])
  print(getpars(fit.modHMMs.depmix.paleo)[i+20])
}

#Row 2 = State 2 
for (i in 31:35){
  print(getpars(fit.modHMMs.depmix.paleo)[i])
  print(getpars(fit.modHMMs.depmix.paleo)[i+5])
  print(getpars(fit.modHMMs.depmix.paleo)[i+10])
  print(getpars(fit.modHMMs.depmix.paleo)[i+15])
  print(getpars(fit.modHMMs.depmix.paleo)[i+20])
}

#Row 3 = State 3 
for (i in 56:60){
  print(getpars(fit.modHMMs.depmix.paleo)[i])
  print(getpars(fit.modHMMs.depmix.paleo)[i+5])
  print(getpars(fit.modHMMs.depmix.paleo)[i+10])
  print(getpars(fit.modHMMs.depmix.paleo)[i+15])
  print(getpars(fit.modHMMs.depmix.paleo)[i+20])
}


#Row 4 = State 4 
for (i in 81:85){
  print(getpars(fit.modHMMs.depmix.paleo)[i])
  print(getpars(fit.modHMMs.depmix.paleo)[i+5])
  print(getpars(fit.modHMMs.depmix.paleo)[i+10])
  print(getpars(fit.modHMMs.depmix.paleo)[i+15])
  print(getpars(fit.modHMMs.depmix.paleo)[i+20])
}

#Row 5 = State 5 
for (i in 106:110){
  print(getpars(fit.modHMMs.depmix.paleo)[i])
  print(getpars(fit.modHMMs.depmix.paleo)[i+5])
  print(getpars(fit.modHMMs.depmix.paleo)[i+10])
  print(getpars(fit.modHMMs.depmix.paleo)[i+15])
  print(getpars(fit.modHMMs.depmix.paleo)[i+20])
}
#########################################################################

#Example- create transition matrices for the 621 years 

#Create the date simulation 
first.year=1400
first.month=11
first.day=1
last.year=2017
last.month=4
last.day=30
iter=1

dates.sim <- create.date.seq.for.season(first.year=first.year,first.month=first.month,first.day=first.day,
                                        last.year=last.year,last.month=last.month,last.day=last.day,no.leap=FALSE)[[1]]

datetxt <- as.Date(dates.sim)
df <- data.frame(date = datetxt,
                 year = as.numeric(format(datetxt, format = "%Y")),
                 month = as.numeric(format(datetxt, format = "%m")),
                 day = as.numeric(format(datetxt, format = "%d")))

df$PC1=0
df$PC2=0
df$PC3=0
df$PC4=0

covariates=data.frame(yr_recon,pred.pcs[,1],pred.pcs[,2],pred.pcs[,3],pred.pcs[,4])
colnames(covariates)=c("Year","PC1","PC2","PC3","PC4")



for (i in 1:(dim(covariates)[1])){
  for (j in 1:(dim(df)[1])){
    if (df$year[j]==covariates$Year[i] & (df$month[j]==1|df$month[j]==2|df$month[j]==3|df$month[j]==4)){
      df$PC1[j]=covariates$PC1[i]
      df$PC2[j]=covariates$PC2[i]
      df$PC3[j]=covariates$PC3[i]
      df$PC4[j]=covariates$PC4[i]
    }
    if (df$year[j]==(covariates$Year[i]-1) & (df$month[j]==11|df$month[j]==12)){
      df$PC1[j]=covariates$PC1[i]
      df$PC2[j]=covariates$PC2[i]
      df$PC3[j]=covariates$PC3[i]
      df$PC4[j]=covariates$PC4[i]
    }
  }
}


#Store the daily probability matrixes in a list#
n=dim(df)[[1]]
mm<-matrix(list(), 1,n)


for (j in 1:n){
  transition_matrix=matrix(, nrow = 5, ncol = 5)
  for (i in 6:10){
    transition_matrix[1,i-5]=getpars(fit.modHMMs.depmix.paleo)[i]+(getpars(fit.modHMMs.depmix.paleo)[i+5]*df$PC1[j])+ (getpars(fit.modHMMs.depmix.paleo)[i+10]*df$PC2[j])+(getpars(fit.modHMMs.depmix.paleo)[i+15]*df$PC3[j])+(getpars(fit.modHMMs.depmix.paleo)[i+20]*df$PC4[j])
  }
  denominator=sum(exp(transition_matrix[1,]))
  for (i in 6:10){
    transition_matrix[1,i-5]=exp(transition_matrix[1,i-5])/denominator
  }
  
  for (i in 31:35){
    transition_matrix[2,i-30]=getpars(fit.modHMMs.depmix.paleo)[i]+(getpars(fit.modHMMs.depmix.paleo)[i+5]*df$PC1[j])+ (getpars(fit.modHMMs.depmix.paleo)[i+10]*df$PC2[j])+(getpars(fit.modHMMs.depmix.paleo)[i+15]*df$PC3[j])+(getpars(fit.modHMMs.depmix.paleo)[i+20]*df$PC4[j])
  }
  denominator=sum(exp(transition_matrix[2,]))
  for (i in 31:35){
    transition_matrix[2,i-30]=exp(transition_matrix[2,i-30])/denominator
  }
  for (i in 56:60){
    transition_matrix[3,i-55]=getpars(fit.modHMMs.depmix.paleo)[i]+(getpars(fit.modHMMs.depmix.paleo)[i+5]*df$PC1[j])+ (getpars(fit.modHMMs.depmix.paleo)[i+10]*df$PC2[j])+(getpars(fit.modHMMs.depmix.paleo)[i+15]*df$PC3[j])+(getpars(fit.modHMMs.depmix.paleo)[i+20]*df$PC4[j])
    
  }
  denominator=sum(exp(transition_matrix[3,]))
  for (i in 56:60){
    transition_matrix[3,i-55]=exp(transition_matrix[3,i-55])/denominator
    
  }
  for (i in 81:85){
    transition_matrix[4,i-80]=getpars(fit.modHMMs.depmix.paleo)[i]+(getpars(fit.modHMMs.depmix.paleo)[i+5]*df$PC1[j])+ (getpars(fit.modHMMs.depmix.paleo)[i+10]*df$PC2[j])+(getpars(fit.modHMMs.depmix.paleo)[i+15]*df$PC3[j])+(getpars(fit.modHMMs.depmix.paleo)[i+20]*df$PC4[j])
    
  }
  denominator=sum(exp(transition_matrix[4,]))
  for (i in 81:85){
    transition_matrix[4,i-80]=exp(transition_matrix[4,i-80])/denominator
    
  }
  
  for (i in 106:110){
    transition_matrix[5,i-105]=getpars(fit.modHMMs.depmix.paleo)[i]+(getpars(fit.modHMMs.depmix.paleo)[i+5]*df$PC1[j])+ (getpars(fit.modHMMs.depmix.paleo)[i+10]*df$PC2[j])+(getpars(fit.modHMMs.depmix.paleo)[i+15]*df$PC3[j])+(getpars(fit.modHMMs.depmix.paleo)[i+20]*df$PC4[j])
    
  }
  denominator=sum(exp(transition_matrix[5,]))
  for (i in 106:110){
    transition_matrix[5,i-105]=exp(transition_matrix[5,i-105])/denominator
    
  }
  mm[[j]]=transition_matrix
  print(j)
}
################################################Define a function that simulates Monte Carlo Chains#####################################################################################
simulate.mc <- function(mcObject,num.states,dates.sim,last.month,last.day,n.sim,iter) {
  
  #this function will simulate the Markov chain iter times
  
  #Arguments:
  #mcObject = a Markov chain object from the markovchain package
  #num.states = the number of states 
  #date.sim = a time series of month values for the simulation period
  #last.month = last month of the season
  #last.day = last day of the last month of the season
  #iter = the number of iterations
  
  #day and month sequences for the simulation period
  days.sim <- as.numeric(format(dates.sim,"%d"))
  months.sim <- as.numeric(format(dates.sim,"%m"))
  n.sim <- length(dates.sim)  #length of simulation
  
  final.mc.sim <- mclapply(X=1:iter,mc.preschedule=TRUE,mc.cores=1,FUN=function(i){  
    
    mc.sim <- as.numeric(rmarkovchain(n=1,object=mcObject[[i]][[1]]))
    end.state <- paste(mc.sim[length(mc.sim)])
    for (t in 1:n.sim) {
      mc.sim <- c(mc.sim,as.numeric(rmarkovchain(n=1,object=mcObject[[i]][[t]],t0=end.state)))
      end.state <- paste(mc.sim[length(mc.sim)])
      if(months.sim[t]==last.month & days.sim[t]==last.day) {end.state <- paste(sample(1:num.states,size=1))}
    }    
    
    #here is the final mc simulation
    final.mc.sim.iter <- mc.sim[2:(n.sim+1)]
    return(final.mc.sim.iter)
  }
  )
  return(final.mc.sim)
  
}


####################Create markov chain list############################
n.sim=n
iter=1

mcWeather.Regime <- mclapply(X=1:iter,mc.preschedule=TRUE,mc.cores=1,FUN=function(j){ 
  
  mcWeather.Regime.time.varying <- mclapply(X=1:n.sim,mc.preschedule=TRUE,mc.cores=1,FUN=function(t){
    tr.prob=as.matrix(mm[[t]])
    mcWeather.Regime.time.varying.out <- new("markovchain", states = c("1","2","3","4","5"),
                                             transitionMatrix = tr.prob, name = paste("Weather.Regime",t,sep=""))    
    return(mcWeather.Regime.time.varying.out)
  }
  )
  mcWeather.Regime.final <- new("markovchainList",markovchains = mcWeather.Regime.time.varying, name = "Weather.Regime.nh")
  return(mcWeather.Regime.final)
  
}
)

##############################################################Simulate Markov States######################################
#Create the date simulation 
first.year=1400
first.month=11
first.day=1
last.year=2017
last.month=3
last.day=31
iter=1

simulate.mc <- function(mcObject,num.states,dates.sim,last.month,last.day,n.sim,iter) {
  
  #this function will simulate the Markov chain iter times
  
  #Arguments:
  #mcObject = a Markov chain object from the markovchain package
  #num.states = the number of states 
  #date.sim = a time series of month values for the simulation period
  #last.month = last month of the season
  #last.day = last day of the last month of the season
  #iter = the number of iterations
  
  #day and month sequences for the simulation period
  days.sim <- as.numeric(format(dates.sim,"%d"))
  months.sim <- as.numeric(format(dates.sim,"%m"))
  n.sim <- length(dates.sim)  #length of simulation
  
  final.mc.sim <- mclapply(X=1:iter,mc.preschedule=TRUE,mc.cores=1,FUN=function(i){  
    
    mc.sim <- as.numeric(rmarkovchain(n=1,object=mcObject[[i]][[1]]))
    end.state <- paste(mc.sim[length(mc.sim)])
    for (t in 1:n.sim) {
      mc.sim <- c(mc.sim,as.numeric(rmarkovchain(n=1,object=mcObject[[i]][[t]],t0=end.state)))
      end.state <- paste(mc.sim[length(mc.sim)])
      if(months.sim[t]==last.month & days.sim[t]==last.day) {end.state <- paste(sample(1:num.states,size=1))}
    }    
    
    #here is the final mc simulation
    final.mc.sim.iter <- mc.sim[2:(n.sim+1)]
    return(final.mc.sim.iter)
  }
  )
  return(final.mc.sim)
  
}


mcsimulations = data.frame(matrix(vector(), 111858, 100))
for (i in 1:100){
  #start_time = Sys.time()
  mcsimulations[,i] <- simulate.mc(mcObject=mcWeather.Regime,num.states=num.states,
                                         dates.sim=dates.sim,last.month=last.month,last.day=last.day,iter=iter)
  #end_time = Sys.time()
  #end_time - start_time
}


