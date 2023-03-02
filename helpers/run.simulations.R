
run.simulations <- function(weather.state.assignments,markov.chain.sim,
                            predictor,prcp.basin,dates.weather,
                            first.month,last.month,num.iter,dates.sim,
                            months,window.size,
                            prcp.site,tmax.site,tmin.site,resampled.dates) {
  
  #this function runs the daily weather generator iter times given iter simulated Markov chains
  
  #Arguments:
  #weather.state.assignments = the historical assignments of weather states, cut to the length of the weather record
  #markov.chain.sim = a time series of weather states for the simulation
  #basin.prcp =  time series of basin-averaged precipitation
  #dates.weather = time series of historic dates for the basin weather
  #first.month = the first month of the record over which to simulate
  #last.month = the last month of the record over which to simulate
  #num.iter = number of iterations (model simulations)
  #dates.sim = a vector of dates over the simulation period
  #months = a vector of the calendar months that define the season (max length == 12) 
  #window.size = the window size, by calendar month, to use for bootstrapping based on days into the season  
  #prcp.site = observed matrix of prcp
  #tmax.site = observed matrix of tmax
  #tmin.site = observed matrix of tmin
  ####################################################################
  
  
  
  #####################Run Simulations################################
  
  #create simulations
  mc.sim <- list()
  resampled.date.sim <- list()
  prcp.site.sim <- list()
  tmax.site.sim <- list()
  tmin.site.sim <- list()
  for (k in 1:num.iter) {

    my.systime <- Sys.time()
    my.sim <- wgen.simulator(weather.state.assignments=weather.state.assignments,mc.sim=markov.chain.sim[[k]],
                             predictor=predictor,prcp.basin=prcp.basin,dates.weather=dates.weather,
                             first.month=first.month,last.month=last.month,dates.sim=dates.sim,
                             months=months,window.size=window.size,
                             prcp.site=prcp.site,tmax.site=tmax.site,tmin.site=tmin.site,resampled.dates)
    print(paste(k,":", Sys.time()-my.systime))
    mc.sim[[k]] <- my.sim[[1]]
    resampled.date.sim[[k]] <- my.sim[[2]]
    prcp.site.sim[[k]] <- round(my.sim[[3]],2)   #round to limit memory storage
    tmax.site.sim[[k]] <- my.sim[[4]]
    tmin.site.sim[[k]] <- my.sim[[5]]
  }
  
  return(list(mc.sim,resampled.date.sim,prcp.site.sim,tmax.site.sim,tmin.site.sim))
}
####################################################################
