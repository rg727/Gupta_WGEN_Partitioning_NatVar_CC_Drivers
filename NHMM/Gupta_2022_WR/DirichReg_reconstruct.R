
DirichReg_reconstruct <- function(param,WR.PCs,num.states,WR_frac,WR_frac.DR,yr,yr_recon,spi_recon.pcs,use.smooth) {
  
  num.pcs <- round(param[1])        #number of SPI PCs to use in reconstruction
  cur.span <- param[2]              #key parameter for smoothing of SPI. set to enable a 10-year smoothing
  
  n.year <- length(yr)
  n.year.recon <- length(yr_recon)
  
  #build Dirichlet regression model on WR PCs
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
  pred.pcs <- array(NA,c(n.year.recon,num.states))
  my.lm <- list(num.states)
  for (i in 1:num.states) {
    cur.var <- WR.PCs[,i]
    fit.data <- data.frame('cur.var'=cur.var,'x'=spi.pcs.instr)
    my.lm[[i]] <- lm(cur.var~.,data=fit.data)
    pred.pcs[,i] <- predict(my.lm[[i]],newdata=data.frame('x'=spi.pcs.recon))    
  }
  
  WR.pred <- predict(my.DirichReg,newdata=data.frame("PC1"=pred.pcs[,1],"PC2"=pred.pcs[,2],"PC3"=pred.pcs[,3],"PC4"=pred.pcs[,4]))
  
  my.out <- list(WR.pred,my.DirichReg,my.lm)
  return(my.out)
  
}

