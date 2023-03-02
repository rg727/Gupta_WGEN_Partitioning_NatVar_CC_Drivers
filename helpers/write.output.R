
write.output <- function(prcp.site.sim,tmax.site.sim,tmin.site.sim,mc.sim,resampled.date.sim,dates.sim,file.suffix) {
  
  #this function writes out the simulated data as .RData files
  
  #Arguments:
  #prcp.site.sim = a list of length iter, with each list element containing a matrix of simulated daily prcp (rows: days; cols: sites)
  #tmax.site.sim = a list of length iter, with each list element containing a matrix of simulated daily tmax (rows: days; cols: sites)
  #tmin.site.sim = a list of length iter, with each list element containing a matrix of simulated daily tmin (rows: days; cols: sites)
  #mc.sim = a list of length iter, with each list element containing a Markov chain vector
  #resampled.date.sim = a list of length iter, with each list element containing a vector of resampled dates
  #dates.sim = a vector of dates over the simulation period
  #file.suffix = a string to append to the end of each file name to track perturbations in each set of simulations
  
  save(prcp.site.sim,file= paste("./Data/simulated.data.files/prcp.site.sim",file.suffix,".RData",sep=""))          
  save(tmax.site.sim,file= paste("./Data/simulated.data.files/tmax.site.sim",file.suffix,".RData",sep=""))          
  save(tmin.site.sim,file= paste("./Data/simulated.data.files/tmin.site.sim",file.suffix,".RData",sep=""))  
  save(mc.sim,file = paste("./Data/simulated.data.files/mc.sim",file.suffix,".RData",sep=""))
  save(resampled.date.sim,file=  paste("./Data/simulated.data.files/resampled.date.sim.sim",file.suffix,".RData",sep=""))  
  save(dates.sim,file = paste("./Data/simulated.data.files/dates.sim",file.suffix,".RData",sep=""))
  
}  