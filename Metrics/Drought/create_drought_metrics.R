
####Tuolumne drought dynamics### 

#Drought Occcurence 

streamflow_0T_0CC_storage=data.frame(matrix(ncol = 50, nrow = 588))

#read in paleo streamflow  
streamflow_0T_0CC=read.table("streamflow_0T_0CC.txt",header=TRUE)

#create monthly data

monthly_data=aggregate(streamflow_0T_0CC, by=list(streamflow_0T_0CC$sim_datemat_1,streamflow_0T_0CC$sim_datemat_2), FUN=mean)

#Create a date column 
monthly_data$Date=as.Date(with(monthly_data, paste(monthly_data$sim_datemat_1,monthly_data$sim_datemat_2,monthly_data$sim_datemat_3,sep="-")), "%Y-%m-%d")
#Now sort by Date column
monthly_data_sorted=monthly_data[order(monthly_data$Date),]

streamflow_0T_0CC=monthly_data_sorted


#We need to fit the distribution across all the ensemble members (since SSI is a relative index)

#This is an R version of the Matlab script from the SDAT toolbox (Farahmand & AghaKouchak, 2015) 
overall_MW=numeric((7410*50))

count=1
for(w in 9:58){
  moving_window=c(monthly_data_sorted[,w])
  overall_MW[count:(count+7410-1)]=moving_window
  count=count+7410
}


#Calculate SSI for moving window
sc=6
n=length(overall_MW);

SI=numeric(n);



if (length(overall_MW)/length(overall_MW) != 1){
  SI[n]=NA
}else{
  SI[1:(sc-1)]=NA
}


A1=matrix(0,370495, 6)
for (i in 1:sc){  
  A1[ ,i]=overall_MW[i:(length(overall_MW)-sc+i)]
}

Y=rowSums(A1)

#Compute the SPI or SSI

nn=length(Y)
SI1=numeric(nn)

for (k in 1:12){
  
  sequence=seq(k,nn,12)
  d=Y[sequence]
  #compute the empirical probability
  nnn=length(d)
  bp=numeric(nnn)
  
  
  for (i in 1:nnn){
    bp[i]=sum(d<=d[i])
  }
  
  y=(bp-0.44)/(nnn+0.12)
  sequence=seq(k,nn,12)
  SI1[sequence]=y
  print(k)
}

SI1=qnorm(SI1)
#output                   
SI[sc:length(SI)]=SI1


##################################Calculate Drought Occurrence (30-year Moving Window)###############################
temp=matrix(SI,nrow=7410)

for (w in 1:50){
  
  a=seq(1,7050,12)
  row=1
  
  
  for (p in a){
    
    moving_window=temp[p:(p+359),w]
    temp_min<-vector()
    temp_min <- length(moving_window[moving_window < -1])
    temp_min=(temp_min)/360*100
    
    
    
    streamflow_0T_0CC_storage[row,w]=as.numeric(temp_min[1])
    
    row=row+1
  }
}


#Concatenate the rows of the storage matrices to create a long vector
streamflow_0T_0CC_vector=as.vector(as.matrix(streamflow_0T_0CC_storage))
years=rep(c(1399:1986),50)
dataframe_0T_0CC_final=data.frame(years,streamflow_0T_0CC_vector)
colnames(dataframe_0T_0CC_final)=c('years','percentage')


write.csv(dataframe_0T_0CC_final,"drought_occurrence.csv",row.names = FALSE)




##################################Calculate Drought Severity (30-Year Moving Window)###############################
for (w in 1:50){
  
  a=seq(1,7050,12)
  row=1
  
  
  for (p in a){
    
    moving_window=temp[p:(p+359),w]
    temp_min=min(moving_window,na.rm=TRUE)
    
    
    
    streamflow_0T_0CC_storage[row,w]=as.numeric(temp_min[1])
    
    row=row+1
  }
}


#Concatenate the rows of the storage matrices to create a long vector
streamflow_0T_0CC_vector=as.vector(as.matrix(streamflow_0T_0CC_storage))
years=rep(c(1399:1986),50)
dataframe_0T_0CC_final=data.frame(years,streamflow_0T_0CC_vector)
colnames(dataframe_0T_0CC_final)=c('years','min_ssi')


write.csv(dataframe_0T_0CC_final,"drought_severity.csv",row.names = FALSE)

##################################Calculate Drought Duration (30-Year Moving Window)###############################

for (w in 1:50){
  
  a=seq(1,7050,12)
  row=1
  
  
  for (p in a){
    condition <- function(x) x <= -1.5
    moving_window=temp[p:(p+359),w]
    # we need to filter those rle results for TRUE:
    r <- rle(condition(moving_window))
    r$length[r$values == TRUE]
    # The answer is the max of the latter
    temp_min=max(r$length[r$values],na.rm=TRUE)
    streamflow_0T_0CC_storage[row,w]=as.numeric(temp_min)
    
    row=row+1
  }
}


#Concatenate the rows of the storage matrices to create a long vector
streamflow_0T_0CC_vector=as.vector(as.matrix(streamflow_0T_0CC_storage))
years=rep(c(1399:1986),50)
dataframe_0T_0CC_final=data.frame(years,streamflow_0T_0CC_vector)
colnames(dataframe_0T_0CC_final)=c('years','months')


write.csv(dataframe_0T_0CC_final,"drought_duration.csv",row.names = FALSE)


