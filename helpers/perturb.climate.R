
perturb.climate <- function(prcp.site.sim,tmax.site.sim,tmin.site.sim,
                            emission.fits.site,prct.prob,months,dates.sim,n.sites,qq,perc.mu,perc.q,S,jitter,temp.change,iter) {
  
  
  #this function takes the simulated prcp, tmax, and tmin from the daily weather generator and perturbs the data to 
  #1) impose thermodynamic climate changes
  #2) impose copula-based jitters
  
  #Arguemnts:
  #prcp.site.sim = a list of length iter, with each list element containing a matrix of simulated daily prcp (rows: days; cols: sites)
  #tmax.site.sim = a list of length iter, with each list element containing a matrix of simulated daily tmax (rows: days; cols: sites)
  #tmin.site.sim = a list of length iter, with each list element containing a matrix of simulated daily tmin (rows: days; cols: sites)
  #emission.fits.site = a list of length 2 (gamma and gpd fits), with each list containing an array of dimension (n.parameters x n.months x n.sites) that contains emission distribution parameters for each site and and month
  #emission.type = a string specifying the emission distribution
  #months = a vector of the calendar months included in each year
  #dates.sim = a vector of dates over the simulation period
  #n.sites = the number of sites across the basin
  #qq = quantile to anchor the CC-scaling
  #perc.mu = percent change in the mean for each month of non-zero prcp for CC-scaling
  #perc.q = percent change in the qqth quantile for each month of non-zero prcp for CC-scaling
  #S = the spearman correlation matrix (n.site+1 x nsite+1) between basin averaged precipitation and site precipitation
  #jitter = TRUE/FALSE indicating whether we randomly perturb the non-exceedance probabilities at each site conditional on the basin average value?
  #temp.change = a step change to apply to all temperatures (tmax and tmin)
  #iter = the number of iterations for the simulation

  #create time series of gamma parameters for both original and CC-scaled gamma distribution
  emission.fit.old.new1 <- CC.scale(emission.fits.site=emission.fits.site,months=months,dates.sim=dates.sim,
                                      n.sites=n.sites,qq=qq,perc.mu=perc.mu,perc.q=perc.q)
  emission.old1 <- emission.fit.old.new1[[1]]
  emission.new1 <- emission.fit.old.new1[[2]]
  emission.old2 <- emission.fits.site[[2]]
  #loop through each of the iterations and impose the appropriate scaling change
  for (j in 1:iter) {  
    my.systime <- Sys.time()
    #perturb the simulated time series of precipitation at each site
    prcp.site.sim[[j]] <- quantile.mapping(prcp.site=prcp.site.sim[[j]],
                                           S=S,prct.prob=prct.prob,perc.q=perc.q,
                                           emission.old1=emission.old1,
                                           emission.old2=emission.old2,
                                           emission.new1=emission.new1,
                                           n.sites=n.sites,months=months,dates.sim=dates.sim,jitter=jitter)
    
    #A placeholder for a spot to add in temperature trends  
    tmin.site.sim[[j]] <- tmin.site.sim[[j]] + temp.change
    tmax.site.sim[[j]] <- tmax.site.sim[[j]] + temp.change
    print(paste("Quant map",j,":", Sys.time()-my.systime))
  }
  return(list(prcp.site.sim,tmax.site.sim,tmin.site.sim))
}