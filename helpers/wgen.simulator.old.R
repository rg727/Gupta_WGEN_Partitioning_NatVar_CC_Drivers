wgen.simulator <- function(weather.state.assignments,mc.sim,
                           predictor,prcp.basin,dates.weather,
                           first.month,last.month,dates.sim,
                           months,window.size,
                           prcp.site,tmax.site,tmin.site,resampled.dates) {
  
  #This function uses a simulation of weather states to conduct 
  #the primary daily weather generation, which is based on a block-bootstrap
  #procedure. It returns the simulated weather regimes it used, the resampled
  #dates from the block bootstrap, and matricies of resampled data for prcp, tmax, and tmin
  #at all sites in the basin.
  
  #Arguments:
  #weather.state.assignments = observed weather regimes for the historical period
  #mc.sim = a simulation of weather regimes to use in the daily simulation
  #prcp.basin = a time series of basin averaged precipitation (same length as weather.state.assignments)
  #dates.weather = a vector of dates associated with the historical record (same length as weather.state.assignments)
  #first.month = first month of the season
  #last.month = last month of the season
  #dates.sim = a vector of dates associated with the simulation (same length as mc.sim)
  #months = a vector of the calendar months that define the season (max length == 12) 
  #window.size = the window size, by calendar month, to use for bootstrapping based on days into the season
  #prcp.site = a matrix of observed daily precipitation at all sites across the basin
  #tmax.site = a matrix of observed daily tmax at all sites across the basin
  #tmin.site = a matrix of observed daily tmin at all sites across the basin
  

  prcp.site.sim <- prcp.site[as.numeric(unlist(resampled.dates)),]
  tmax.site.sim <- tmax.site[as.numeric(unlist(resampled.dates)),]
  tmin.site.sim <- tmin.site[as.numeric(unlist(resampled.dates)),]
  sampled.date.sim=as.numeric(unlist(resampled.dates))
  return(list(mc.sim,sampled.date.sim,prcp.site.sim,tmax.site.sim,tmin.site.sim))
  
}