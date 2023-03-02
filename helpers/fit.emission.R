
fit.emission <- function(prct.prob,prcp.site,months,months.weather,n.sites) {
  
  #identify indexes for below/above threshold per site for gamma and gpd fits
  quantile.nz <- function(x,p=prct.prob){quantile(x[x!=0],p,na.rm=T)}
  thshd.prcp <- apply(prcp.site,2,FUN=quantile.nz)
  
  idx1 <- sapply(1:n.sites,function(x,v){
    (prcp.site[,x]!=0 & !is.na(prcp.site[,x]) & prcp.site[,x]<v[x])
  },v=thshd.prcp)
  idx2 <- sapply(1:n.sites,function(x,v){
    (prcp.site[,x]!=0 & !is.na(prcp.site[,x]) & prcp.site[,x]>=v[x])
  },v=thshd.prcp)
  
  #1: gamma
  prcp.site00 <- prcp.site*idx1 # tagging those ones for fitting gamma (0 or values below threshold)
  emission.fit.site1 <- apply(prcp.site00,2,function(x,m,mm) {
    sapply(m,function(m) {fitdistr(as.numeric(x[x!=0 & mm==m & !is.na(x)]),'gamma')$estimate})
  },m=months, mm=months.weather)
  
  dim(emission.fit.site1) <- c(2,length(months),n.sites)   #reformat dimensions so that shape/rate down each row, months in each column, 3rd dimension indexes the sites
  
  #2: gpd
  prcp.site0 <- rbind(prcp.site*idx2,thshd.prcp) # tagging those ones for fitting gpd (0 or values above threshold)
  emission.fit.site2 <- apply(prcp.site0,2,function(x) {eva::gpdFit(as.numeric(x[x!=0 & !is.na(x)]),threshold = tail(x,n=1))$par.ests})

  emission.fits.site <- list("emission.fit.site1"=emission.fit.site1,
                             "emission.fit.site2"=emission.fit.site2)
  return(emission.fits.site)
}