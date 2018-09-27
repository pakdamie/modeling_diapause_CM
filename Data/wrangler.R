###WRANGLER
###RUN THIS BEFORE RUNNING THE POMP MODEL
###THE ACTUAL DATA

library(dplyr)
library(weathermetrics)

###This is for the actual data needed for the fitting
##################################################################
###This reads in field trap data of Cydia pomonella from 1984-2016
CM_Dat <- read.csv("CM.2.csv")
###Change the format of the date into actual dates
CM_Dat$Date <- as.Date(CM_Dat$Date, format = '%A, %b %d, %Y')
###We only want to get the average across sites
CM_Dat$Mean <- rowMeans(CM_Dat[,2: 21],na.rm=TRUE)

###Missing some dates, so we create a new one
AS.DATE<-data.frame(Date= seq(as.Date("1984-01-01", format ='%Y-%m-%d'), 
                              as.Date("2016-12-31",format ='%Y-%m-%d'),'days'))

###Join the full dates with the field data.

CM_JULIAN <- left_join(AS.DATE,CM_Dat[,c(1,22)])

###The "Julian Day" is going to start from 1 to 12054
CM_JULIAN$Time <- seq(1, 12054)
###We're going to name the actual field mean data 'ARC' 
colnames(CM_JULIAN)[2]<-'ARC'
###If one of the entries is NAN- then turn into NA
CM_JULIAN$ARC[is.nan(CM_JULIAN$ARC)==TRUE]<-NA

########################################################################

###This loads in the temperature 


TEMP_DAT <- read.csv("big.temp_week.csv",stringsAsFactors = FALSE)
TEMP_DAT$Date <- as.Date(TEMP_DAT$Date, format = '%m/%d/%Y')
TEMP_DAT$Year <- as.numeric(format(TEMP_DAT$Date, '%Y'))

T_1984_2016 <- subset(TEMP_DAT, TEMP_DAT$Year > 1983 &
                        TEMP_DAT$Year < 2017)
T_1984_2016$tmin <- fahrenheit.to.celsius(T_1984_2016$tmin, round=2)
T_1984_2016$tmax <- fahrenheit.to.celsius(T_1984_2016$tmax, round=2)

T_AVG <- (T_1984_2016$tmin + T_1984_2016$tmax)/2

######################################################################

###This loads in the geopshere

PHOTO_PERIOD_GEOSPHERE <- geosphere::daylength(40, 1:366)
DIFF <- diff(PHOTO_PERIOD_GEOSPHERE)




covartable <- data.frame(
  time=CM_JULIAN$Time,
  TT=(T_1984_2016$tmax +T_1984_2016$tmin)/2,
  PP_DIFF= rep(DIFF, length =12054),
  ARC =CM_JULIAN$ARC )

save(covartable, file='covar.RData')
