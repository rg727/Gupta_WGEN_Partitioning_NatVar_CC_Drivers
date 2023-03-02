CC.scale <- function(emission.fits.site,months,dates.sim,n.sites,qq,perc.mu,perc.q) {
  
  #this function takes fitted emission distributions to each site by month, and returns those parameters as 
  #an array of parameters with dimensions (n.sim, n.sites,n.parameters) 
  #it also enables a Clausius-Clapyeron scaling for the mean and the qqth quantile, and returns a similar
  #array of parameters for emission distributions with that scaling
  
  #Arguments:
  #emission.fits.site = a list of length 2 (gamma and gpd fits), with each list containing an array of dimension (n.parameters x n.months x n.sites) that contains emission distribution parameters for each site and and month
  #emission.type = a string specifying the emission distribution
  #months = a vector of calendar months that define the season
  #dates.sim = a vector of dates associated with the simulation period
  #n.sites = the number of individual sites
  #qq = the quantile for CC scaling
  #perc.mu = the percentage by which to scale the mean
  #perc.q = the percentage by which to scale the qqth quantile for CC scaling
  
  emission.type <- "gamma"
  
  #month sequence for the simulation period
  months.sim <- as.numeric(format(dates.sim,"%m"))
  n.sim <- length(dates.sim)  #length of the simulation period
  
  #if the quantile isn't specified, we use the maximum empirical quantile
  if(is.null (qq)) {
    EP <- (1:n.sim)/(n.sim+1)
    qq <- max(EP)
  }
  #original fits  
  emission.old1 <- emission.fits.site[[1]]
  
  # monthly-only #
  #find the mean and qqth quantile of the different distributions, and concatenate them with the parameters into a single array
  q.old <- apply(emission.old1,c(2,3),function(x) {qgamma(qq,shape=x[1],rate=x[2])})
  mu.old <- apply(emission.old1,c(2,3),function(x) {x[1]/x[2]})
  q.mu.old <- abind(q.old,mu.old,emission.old1,along=1)
  
  #adjust distribution for new mean, qqth quantile
  emission.new1 <- apply(q.mu.old,c(2,3), function(x) {
    start.par <- 1
    lowerb <- 0.0001
    upperb <- 10
    opt <- optim(par=start.par,CC.scale.obj.fun,q.old=x[1],mu.old=x[2],
                 param.old=x[3:4],qq=qq,
                 perc.q=perc.q,perc.mu=perc.mu,emission.type=emission.type,
                 method="L-BFGS-B",lower=lowerb,upper=upperb)
    shape.new <- x[3]*perc.mu*opt$par
    rate.new <- x[4]*opt$par
    return(c(shape.new,rate.new))
  })
  return(list(emission.old1,emission.new1))
  
}