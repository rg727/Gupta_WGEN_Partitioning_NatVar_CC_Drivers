library(areaplot)

#Import in the metrics for the 10-year flood

dataframe_0T_0CC=readRDS("dataframe_0T_0CC_final_10yr.rds") #This one is included as a truncated example- has all the ensemble members stacked in a long dataframe
dataframe_1T_0CC=readRDS("dataframe_1T_0CC_final_10yr.rds")
dataframe_2T_0CC=readRDS("dataframe_2T_0CC_final_10yr.rds")
dataframe_3T_0CC=readRDS("dataframe_3T_0CC_final_10yr.rds")
dataframe_4T_0CC=readRDS("dataframe_4T_0CC_final_10yr.rds")

dataframe_0T_0.5CC=readRDS("dataframe_0T_0.5CC_final_10yr.rds")
dataframe_1T_0.5CC=readRDS("dataframe_1T_0.5CC_final_10yr.rds")
dataframe_2T_0.5CC=readRDS("dataframe_2T_0.5CC_final_10yr.rds")
dataframe_3T_0.5CC=readRDS("dataframe_3T_0.5CC_final_10yr.rds")
dataframe_4T_0.5CC=readRDS("dataframe_4T_0.5CC_final_10yr.rds")


dataframe_0T_1CC=readRDS("dataframe_0T_1CC_final_10yr.rds")
dataframe_1T_1CC=readRDS("dataframe_1T_1CC_final_10yr.rds")
dataframe_2T_1CC=readRDS("dataframe_2T_1CC_final_10yr.rds")
dataframe_3T_1CC=readRDS("dataframe_3T_1CC_final_10yr.rds")
dataframe_4T_1CC=readRDS("dataframe_4T_1CC_final_10yr.rds")

dataframe_0T_1.5CC=readRDS("dataframe_0T_1CC_final_10yr.rds")
dataframe_1T_1.5CC=readRDS("dataframe_1T_1CC_final_10yr.rds")
dataframe_2T_1.5CC=readRDS("dataframe_2T_1CC_final_10yr.rds")
dataframe_3T_1.5CC=readRDS("dataframe_3T_1CC_final_10yr.rds")
dataframe_4T_1.5CC=readRDS("dataframe_4T_1CC_final_10yr.rds")


dataframe_0T_2CC=readRDS("dataframe_0T_2CC_final_10yr.rds")
dataframe_1T_2CC=readRDS("dataframe_1T_2CC_final_10yr.rds")
dataframe_2T_2CC=readRDS("dataframe_2T_2CC_final_10yr.rds")
dataframe_3T_2CC=readRDS("dataframe_3T_2CC_final_10yr.rds")
dataframe_4T_2CC=readRDS("dataframe_4T_2CC_final_10yr.rds")


#Set up a dataframe to store the fractions

ANOVA_values=data.frame(matrix(vector(),588,4))
years=c(1399:1986)

#conduct ANOVA decomposition by year 

for (i in 1:588){
  
  newdata_0T_0CC <-dataframe_0T_0CC[ which(dataframe_0T_0CC$years==years[i]),]
  newdata_1T_0CC <-dataframe_1T_0CC[ which(dataframe_1T_0CC$years==years[i]),]
  newdata_2T_0CC <-dataframe_2T_0CC[ which(dataframe_2T_0CC$years==years[i]),]
  newdata_3T_0CC <-dataframe_3T_0CC[ which(dataframe_3T_0CC$years==years[i]),]
  newdata_4T_0CC <-dataframe_4T_0CC[ which(dataframe_4T_0CC$years==years[i]),]
  
  newdata_0T_0.5CC <-dataframe_0T_0.5CC[ which(dataframe_0T_0.5CC$years==years[i]),]
  newdata_1T_0.5CC <-dataframe_1T_0.5CC[ which(dataframe_1T_0.5CC$years==years[i]),]
  newdata_2T_0.5CC <-dataframe_2T_0.5CC[ which(dataframe_2T_0.5CC$years==years[i]),]
  newdata_3T_0.5CC <-dataframe_3T_0.5CC[ which(dataframe_3T_0.5CC$years==years[i]),]
  newdata_4T_0.5CC <-dataframe_4T_0.5CC[ which(dataframe_4T_0.5CC$years==years[i]),]
  
  newdata_0T_1CC <-dataframe_0T_1CC[ which(dataframe_0T_1CC$years==years[i]),]
  newdata_1T_1CC <-dataframe_1T_1CC[ which(dataframe_1T_1CC$years==years[i]),]
  newdata_2T_1CC <-dataframe_2T_1CC[ which(dataframe_2T_1CC$years==years[i]),]
  newdata_3T_1CC <-dataframe_3T_1CC[ which(dataframe_3T_1CC$years==years[i]),]
  newdata_4T_1CC <-dataframe_4T_1CC[ which(dataframe_4T_1CC$years==years[i]),]
  
  newdata_0T_1.5CC <-dataframe_0T_1.5CC[ which(dataframe_0T_1.5CC$years==years[i]),]
  newdata_1T_1.5CC <-dataframe_1T_1.5CC[ which(dataframe_1T_1.5CC$years==years[i]),]
  newdata_2T_1.5CC <-dataframe_2T_1.5CC[ which(dataframe_2T_1.5CC$years==years[i]),]
  newdata_3T_1.5CC <-dataframe_3T_1.5CC[ which(dataframe_3T_1.5CC$years==years[i]),]
  newdata_4T_1.5CC <-dataframe_4T_1.5CC[ which(dataframe_4T_1.5CC$years==years[i]),]
  
  
  newdata_0T_2CC <-dataframe_0T_2CC[ which(dataframe_0T_2CC$years==years[i]),]
  newdata_1T_2CC <-dataframe_1T_2CC[ which(dataframe_1T_2CC$years==years[i]),]
  newdata_2T_2CC <-dataframe_2T_2CC[ which(dataframe_2T_2CC$years==years[i]),]
  newdata_3T_2CC <-dataframe_3T_2CC[ which(dataframe_3T_2CC$years==years[i]),]
  newdata_4T_2CC <-dataframe_4T_2CC[ which(dataframe_4T_2CC$years==years[i]),]
  
  
  newdata=rbind(newdata_0T_0CC,newdata_1T_0CC,newdata_2T_0CC,newdata_3T_0CC,newdata_4T_0CC,newdata_0T_0.5CC,newdata_1T_0.5CC,newdata_2T_0.5CC,newdata_3T_0.5CC,newdata_4T_0.5CC,newdata_0T_1CC,newdata_1T_1CC,newdata_2T_1CC,newdata_3T_1CC,newdata_4T_1CC,newdata_0T_1.5CC,newdata_1T_1.5CC,newdata_2T_1.5CC,newdata_3T_1.5CC,newdata_4T_0.5CCnewdata_0T_2CC,newdata_1T_2CC,newdata_2T_2CC,newdata_3T_2CC,newdata_4T_2CC)
  
  newdata <- newdata[!is.infinite(rowSums(newdata)),]
  
  #conduct the ANOVA
  
  newdata$cc=as.factor(newdata$cc)
  newdata$tt=as.factor(newdata$tt)
  res.aov <- aov(return ~ cc * tt, data = newdata)
  summary(res.aov)
  total_SS=sum(summary(res.aov)[[1]][, 'Sum Sq'])
  
  CC_trend_SS=summary(res.aov)[[1]][1,2]/total_SS
  TT_trend_SS=summary(res.aov)[[1]][2,2]/total_SS
  res_trend_SS=summary(res.aov)[[1]][3,2]/total_SS
  natural_trend_SS=summary(res.aov)[[1]][4,2]/total_SS
  
  ANOVA_values[i,1]=CC_trend_SS
  ANOVA_values[i,2]=TT_trend_SS
  ANOVA_values[i,3]=res_trend_SS
  ANOVA_values[i,4]=natural_trend_SS
  
}

ANOVA_values$Year=c(1399:1976)

#Plot results
cols=c('#BC6C25','#DDA15E','#FEFAE0','#283618','#606C38')

colnames(ANOVA_values)=c("CC Scaling","Temperature Trend","Interactions","Natural Variability","Year")
pdf(file="anova_10yrflood_30yrwindow.pdf")
# Stacked area chart with custom colors
areaplot(x=ANOVA_values$Year,ANOVA_values[,c(1,2,3,4)], prop = TRUE,
         col = cols,legend = FALSE,main="10-Yr Return Period Event (30-Yr Moving Average)",xlab="Year",ylab="Fraction of Variance Explained",
         args.legend = list(x = "topright", cex = 0.65,
                            bg = "white", bty = "o"))
dev.off()

