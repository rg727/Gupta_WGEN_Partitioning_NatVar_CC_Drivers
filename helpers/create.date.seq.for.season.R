
create.date.seq.for.season <- function(first.year,first.month,first.day,last.year,last.month,last.day,no.leap=FALSE) {

  #create a daily series of dates, months and years over a range of years for the season of interest. years are 
  #automatically returned in water years
  
  #Arguments
  #first.year = the first year of the sequence  (in 'water year')
  #first.month = the beginning month of the season
  #first.day = the first day of the beginning month of the season
  #last.year = the last year of the sequence
  #last.month = the last month of the season
  #last.day = the last day of the last month of the season
  
  #use 1 year earlier than first.year if season spans Dec/Jan
  first.year <- first.year; if(first.month>last.month) {first.year <- first.year-1}   
  
  if(first.month>last.month) {
    n.year <- length(first.year:last.year) - 1
    date.seq <- seq(as.Date(paste(first.year,first.month,first.day,sep="-")),
                     as.Date(paste((first.year+1),last.month,last.day,sep="-")),by="day")
    year.iter <- first.year
    for (j in 2:n.year) {
      year.iter <- year.iter + 1
      date.seq <- c(date.seq,seq(as.Date(paste(year.iter,first.month,first.day,sep="-")),
                                   as.Date(paste((year.iter+1),last.month,last.day,sep="-")),by="day"))
    }
  } else {
      n.year <- length(first.year:last.year)
      date.seq <- seq(as.Date(paste(first.year,first.month,first.day,sep="-")),
                       as.Date(paste(first.year,last.month,last.day,sep="-")),by="day")
      year.iter <- first.year
      for (j in 2:n.year) {
        year.iter <- year.iter + 1
        date.seq <- c(date.seq,seq(as.Date(paste(year.iter,first.month,first.day,sep="-")),
                                     as.Date(paste(year.iter,last.month,last.day,sep="-"))
                                     ,by="day"))
      }  
  }

  #if requested, drop leap year days in the full sequence of dates
  if (no.leap) {
    to.drop <- which(format(date.seq,"%m")=="02" & format(date.seq,"%d")=="29")
    date.seq <- date.seq[-to.drop]
  }
  
  #month and water year sequences
  days.seq <- as.numeric(format(date.seq,"%d"))
  months.seq <- as.numeric(format(date.seq,"%m"))
  years.seq <- as.numeric(format(date.seq,"%Y"))
  if(first.month>last.month) {
    years.seq[months.seq>=first.month] <- years.seq[months.seq>=first.month] + 1
  }
  
  return(list(date.seq,years.seq,months.seq,days.seq))
  
}
