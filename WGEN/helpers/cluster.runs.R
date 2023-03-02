cluster.runs <- function (my.clust,my.date,first.month,last.month) {
  
  #This function will return a data.frame with information on cluster runs, including
  #the position of the first and last day of each run, the run length, the month
  #of the first and last day of the run, and the average day of the season for the run
  
  #Arguments:
  #my.clust = a daily vector of states 
  #my.date = a daily vector of dates of the same length as my.clust
  #first.month = the first month of the season
  #last.month = the last month of the season
  
  my.month <- as.numeric(format(my.date,"%m"))
  
  n.clust <- length(my.clust)
  
  first.of.run <- c(1,which(diff(my.clust)!=0) + 1)
  last.of.run <- c(which(diff(my.clust)!=0),n.clust)
  
  #find those events that cross over a water year
  cross <- which(my.month[first.of.run]==last.month & my.month[last.of.run]==first.month)
  while(length(cross)!=0) {
    cur.cross <- cross[1]
    months.of.cross <- my.month[first.of.run[cur.cross]:last.of.run[cur.cross]]
    first.step <- which(diff(months.of.cross)!=0)
    last.step <- length(months.of.cross)-first.step
    first.of.run <- c(first.of.run[1:(cur.cross-1)],first.of.run[cur.cross],first.of.run[cur.cross]+first.step,first.of.run[(cur.cross+1):length(first.of.run)])
    last.of.run <- c(last.of.run[1:(cur.cross-1)],last.of.run[cur.cross]-last.step,last.of.run[cur.cross],last.of.run[(cur.cross+1):length(last.of.run)])
    cross <- which(my.month[first.of.run]==last.month & my.month[last.of.run]==first.month)
  }  
  
  #run lengths and type
  run.length <- last.of.run-first.of.run + 1
  clust.type <- my.clust[first.of.run]
  
  #define the date for the start of the season with a dummy year (1901)
  start.of.season <- as.Date(paste(1901,first.month,1,sep="-"))
  #define dates of the season for that same dummy year, accounting for if the season spans calendar years
  day.of.season <- as.Date(paste(1901,as.numeric(format(my.date,"%m")),as.numeric(format(my.date,"%d")),sep="-"))
  if (first.month>last.month) {
    keep <- which(my.month<first.month)
    day.of.season[keep] <- as.Date(paste(1902,as.numeric(format(my.date[keep],"%m")),as.numeric(format(my.date[keep],"%d")),sep="-"))
  }
  #here we replace any leap years in day.of.season (which will have an NA value) with Feb 28
  day.of.season[is.na(day.of.season)] <- as.Date("1902-2-28")
  #calculate days from the first day of the season, and determine average day into season for each run
  days.into.season <- day.of.season - start.of.season
  days.into.season.first <- days.into.season[first.of.run]
  days.into.season.last <- days.into.season[last.of.run]
  days.into.season.avg <- (days.into.season.first+days.into.season.last)/2
  

  cluster.run <- data.frame('state'=clust.type,'first.run'= first.of.run, 'last.run'=last.of.run,'run.length'=run.length,
                            'month.first'=my.month[first.of.run],'month.last'=my.month[last.of.run],'avg.day.into.season'=days.into.season.avg)
  
  #checks to make sure 1) the first day of a run always immediately follows the last day of the previous run
  #and 2) we dont have a run that spans the end of one water year and the beginning of the next (assuming WYs are separated by unmodeled months)
  #w <- dim(cluster.run)[1]
  #which((cluster.run$first.run[2:w] - cluster.run$last.run[1:(w-1)])!=1)
  #which(cluster.run$month.first==last.month & cluster.run$month.last==first.month)
  
  return(cluster.run)
  
}  