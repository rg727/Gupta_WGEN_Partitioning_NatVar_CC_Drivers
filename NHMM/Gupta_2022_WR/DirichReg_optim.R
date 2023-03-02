nvars <- 2
nobjs <-5
objs=numeric(5)
use.smooth=TRUE


DirichReg_optim <- function(X) {
  
  num.pcs <- round(X[[1]])        #number of SPI PCs to use in reconstruction
  cur.span <- X[[2]]              #key parameter for smoothing of SPI. set to enable a 10-year smoothing
            #key parameter for smoothing of SPI. set to enable a 10-year smoothing
  
  n.year <- length(yr)
  num.breaks <- round(n.year/leave.k.out)
  #set.seed(1)
  cv.blocks <- sample(cut(1:n.year,breaks=num.breaks,labels=1:num.breaks))
    
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
  
  #begin cross validation
  pred.pcs.cv <- array(NA,c(n.year,num.states))
  WR.pred.cv <- array(NA,c(n.year,num.states))
  WR.pred.alpha.cv <- array(NA,c(n.year,num.states))
  for (j in 1:num.breaks) {
    
    my.in <- -which(cv.blocks==j)
    my.out <- which(cv.blocks==j)
    
    #covariate data
    x.in <- spi.pcs.instr[my.in,]
    x.out <- spi.pcs.instr[my.out,]
    
    #predict WR PCs based on SPI PCs
    for (i in 1:num.states) {
      cur.var <- WR.PCs[my.in,i]
      fit.data <- data.frame('cur.var'=cur.var,'x'=x.in)
      my.lm <- lm(cur.var~.,data=fit.data)
      test.data <- data.frame('x'=matrix(x.out,ncol=dim(x.in)[2]))
      names(test.data) <- names(fit.data)[2:dim(fit.data)[2]]
      pred.pcs.cv[my.out,i] <- predict(my.lm,newdata=test.data)        
    }
    
    #build Dirichlet regression model on WR PCs
    WR_frac.df.cv <- data.frame(WR.PCs[my.in,],WR_frac[my.in,])
    suppressWarnings(WR_frac.df.cv$Y <- DR_data(WR_frac.df.cv[,grep("X",names(WR_frac.df.cv))]))
    WR_frac.df.cv <- WR_frac.df.cv[,-grep("X",names(WR_frac.df.cv))]
    my.DirichReg <- DirichReg(Y ~ PC1 + PC2 + PC3 + PC4, data=WR_frac.df.cv)   #currently hardcoded for 4 WR PCs
    
    #predict dirichlet alphas based on predicted WR PCs
    WR.pred.out <- predict(my.DirichReg,newdata=data.frame("PC1"=pred.pcs.cv[my.out,1],"PC2"=pred.pcs.cv[my.out,2],"PC3"=pred.pcs.cv[my.out,3],"PC4"=pred.pcs.cv[my.out,4]),alpha=TRUE)
    WR.pred.cv[my.out,] <- WR.pred.out$mu
    WR.pred.alpha.cv[my.out,] <- WR.pred.out$alpha
  }
  
  log.lik <- ddirichlet(WR_frac.DR, alpha=WR.pred.alpha.cv, log = TRUE, sum.up = TRUE)
  
  #MSE of differnt rolling averages
  MSE <- array(NA,length(smooth.list))
  for (k in 1:length(smooth.list)) {
    kk <- smooth.list[k]
    obs.rollmean <- apply(WR_frac.DR,2,function(x){rollmean(x,kk)})
    pred.rollmean <- apply(WR.pred.cv,2,function(x){rollmean(x,kk)})  
    MSE[k] <- mean((obs.rollmean - pred.rollmean)^2)
  }
  
  objs[1]=MSE[1]
  objs[2]=MSE[2]
  objs[3]=MSE[3]
  objs[4]=MSE[4]
  objs[5]=num.pcs
  return(list(objs))
  
}