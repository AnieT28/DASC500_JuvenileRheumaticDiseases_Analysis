#Prevelance calculations using david's code and unedited dataset.

#Set working directory
setwd("R:/Anie")
getwd()

#Load data
load("DataforAnie - Copy.RData")

#Load relevant libraries
library(epitools)
library(writexl)

#Empty data frame which is to be filled with results
results <- data.frame(Year = integer(),Cases = integer(),Population = integer(),
                      Prevalence = numeric(),Lower95CI = numeric(),Upper95CI = numeric())

#Prevalence 2022  - jSLE
#Inclusion criteria - regstart earlier or equal to prevalence date, and between 0-16 for JIA or 0-18 for jSLE and JDM
tmp1<-subset(JuvenileRheumatic,JuvenileRheumatic$regstartdate<="31/12/2022"&JuvenileRheumatic$yob>=2004&JuvenileRheumatic$yob<=2022)
#Exclude those from the start of prevalence year
tmp2 <- subset (tmp1, !(tmp1$regenddate < "2022-01-01"))
#Calculate prevalence
prev2022<-sum(tmp2$jSLE==1)/nrow(tmp2)*100000
#Print result
prev2022
#Compute 95% CI with pois.exact()
ci <- pois.exact(x=(sum(tmp2$jSLE==1)),pt=(nrow(tmp2)),conf.level = 0.95)
#Lower and Upper 95% CI
lower95 <- ci$lower * 100000
upper95 <- ci$upper* 100000
lower95
upper95

#Add results to empty table
results <- rbind(results, data.frame(Year = 2022, Cases = sum(tmp2$jSLE==1),Population = nrow(tmp2),Prevalence=prev2022, Lower95CI=lower95, Upper95CI=upper95))

#Prevalence 2021
tmp1<-subset(JuvenileRheumatic,JuvenileRheumatic$regstartdate<="31/12/2021"&JuvenileRheumatic$yob>=2003&JuvenileRheumatic$yob<=2021)
tmp2 <- subset (tmp1, !(tmp1$regenddate < "2021-01-01"))
prev2021<-sum(tmp2$jSLE==1)/nrow(tmp2)*100000
prev2021 
ci <- pois.exact(x=(sum(tmp2$jSLE==1)),pt=(nrow(tmp2)),conf.level = 0.95)
lower95 <- ci$lower * 100000
upper95 <- ci$upper* 100000
lower95
upper95

results <- rbind(results, data.frame(Year = 2021, Cases = sum(tmp2$jSLE==1),Population = nrow(tmp2),Prevalence=prev2021, Lower95CI=lower95, Upper95CI=upper95))


#Prevalence 2020
tmp1<-subset(JuvenileRheumatic,JuvenileRheumatic$regstartdate<="31/12/2020"&JuvenileRheumatic$yob>=2002&JuvenileRheumatic$yob<=2020)
tmp2 <- subset (tmp1, !(tmp1$regenddate < "2020-01-01"))
prev2020<-sum(tmp2$jSLE==1)/nrow(tmp2)*100000
prev2020 
ci <- pois.exact(x=(sum(tmp2$jSLE==1)),pt=(nrow(tmp2)),conf.level = 0.95)
lower95 <- ci$lower * 100000
upper95 <- ci$upper* 100000
lower95
upper95

results <- rbind(results, data.frame(Year = 2020, Cases = sum(tmp2$jSLE==1),Population = nrow(tmp2),Prevalence=prev2020, Lower95CI=lower95, Upper95CI=upper95))

#Prevalence 2019
tmp1<-subset(JuvenileRheumatic,JuvenileRheumatic$regstartdate<="31/12/2019"&JuvenileRheumatic$yob>=2001&JuvenileRheumatic$yob<=2019)
tmp2 <- subset (tmp1, !(tmp1$regenddate < "2019-01-01"))
prev2019<-sum(tmp2$jSLE==1)/nrow(tmp2)*100000
prev2019 

results <- rbind(results, data.frame(Year = 2019, Cases = sum(tmp2$jSLE==1),Population = nrow(tmp2),Prevalence=prev2019, Lower95CI=lower95, Upper95CI=upper95))


#Prevalence 2018
tmp1<-subset(JuvenileRheumatic,JuvenileRheumatic$regstartdate<="31/12/2018"&JuvenileRheumatic$yob>=2000&JuvenileRheumatic$yob<=2018)
tmp2 <- subset (tmp1, !(tmp1$regenddate < "2018-01-01"))
prev2018<-sum(tmp2$jSLE==1)/nrow(tmp2)*100000
prev2018 
ci <- pois.exact(x=(sum(tmp2$jSLE==1)),pt=(nrow(tmp2)),conf.level = 0.95)
lower95 <- ci$lower * 100000
upper95 <- ci$upper* 100000
lower95
upper95

results <- rbind(results, data.frame(Year = 2018, Cases = sum(tmp2$jSLE==1),Population = nrow(tmp2),Prevalence=prev2018, Lower95CI=lower95, Upper95CI=upper95))


#Prevalence 2017
tmp1<-subset(JuvenileRheumatic,JuvenileRheumatic$regstartdate<="31/12/2017"&JuvenileRheumatic$yob>=1999&JuvenileRheumatic$yob<=2017)
tmp2 <- subset (tmp1, !(tmp1$regenddate < "2017-01-01"))
prev2017<-sum(tmp2$jSLE==1)/nrow(tmp2)*100000
prev2017 
ci <- pois.exact(x=(sum(tmp2$jSLE==1)),pt=(nrow(tmp2)),conf.level = 0.95)
lower95 <- ci$lower * 100000
upper95 <- ci$upper* 100000
lower95
upper95

results <- rbind(results, data.frame(Year = 2017, Cases = sum(tmp2$jSLE==1),Population = nrow(tmp2),Prevalence=prev2017, Lower95CI=lower95, Upper95CI=upper95))

#Prevalence 2016
tmp1<-subset(JuvenileRheumatic,JuvenileRheumatic$regstartdate<="31/12/2016"&JuvenileRheumatic$yob>=1998&JuvenileRheumatic$yob<=2016)
tmp2 <- subset (tmp1, !(tmp1$regenddate < "2016-01-01"))
prev2016<-sum(tmp2$jSLE==1)/nrow(tmp2)*100000
prev2016 
ci <- pois.exact(x=(sum(tmp2$jSLE==1)),pt=(nrow(tmp2)),conf.level = 0.95)
lower95 <- ci$lower * 100000
upper95 <- ci$upper* 100000
lower95
upper95

results <- rbind(results, data.frame(Year = 2016, Cases = sum(tmp2$jSLE==1),Population = nrow(tmp2),Prevalence=prev2016, Lower95CI=lower95, Upper95CI=upper95))

#Prevalence 2015
tmp1<-subset(JuvenileRheumatic,JuvenileRheumatic$regstartdate<="31/12/2015"&JuvenileRheumatic$yob>=1997&JuvenileRheumatic$yob<=2015)
tmp2 <- subset (tmp1, !(tmp1$regenddate < "2015-01-01"))
prev2015<-sum(tmp2$jSLE==1)/nrow(tmp2)*100000
prev2015 
ci <- pois.exact(x=(sum(tmp2$jSLE==1)),pt=(nrow(tmp2)),conf.level = 0.95)
lower95 <- ci$lower * 100000
upper95 <- ci$upper* 100000
lower95
upper95

results <- rbind(results, data.frame(Year = 2015, Cases = sum(tmp2$jSLE==1),Population = nrow(tmp2),Prevalence=prev2015, Lower95CI=lower95, Upper95CI=upper95))

#Prevalence 2014
tmp1<-subset(JuvenileRheumatic,JuvenileRheumatic$regstartdate<="31/12/2014"&JuvenileRheumatic$yob>=1996&JuvenileRheumatic$yob<=2014)
tmp2 <- subset (tmp1, !(tmp1$regenddate < "2014-01-01"))
prev2014<-sum(tmp2$jSLE==1)/nrow(tmp2)*100000
prev2014 
ci <- pois.exact(x=(sum(tmp2$jSLE==1)),pt=(nrow(tmp2)),conf.level = 0.95)
lower95 <- ci$lower * 100000
upper95 <- ci$upper* 100000
lower95
upper95

results <- rbind(results, data.frame(Year = 2014, Cases = sum(tmp2$jSLE==1),Population = nrow(tmp2),Prevalence=prev2014, Lower95CI=lower95, Upper95CI=upper95))

#Prevalence 2013
tmp1<-subset(JuvenileRheumatic,JuvenileRheumatic$regstartdate<="31/12/2013"&JuvenileRheumatic$yob>=1995&JuvenileRheumatic$yob<=2013)
tmp2 <- subset (tmp1, !(tmp1$regenddate < "2013-01-01"))
prev2013<-sum(tmp2$jSLE==1)/nrow(tmp2)*100000
prev2013 
ci <- pois.exact(x=(sum(tmp2$jSLE==1)),pt=(nrow(tmp2)),conf.level = 0.95)
lower95 <- ci$lower * 100000
upper95 <- ci$upper* 100000
lower95
upper95

results <- rbind(results, data.frame(Year = 2013, Cases = sum(tmp2$jSLE==1),Population = nrow(tmp2),Prevalence=prev2013, Lower95CI=lower95, Upper95CI=upper95))

#Prevalence 2012
tmp1<-subset(JuvenileRheumatic,JuvenileRheumatic$regstartdate<="31/12/2012"&JuvenileRheumatic$yob>=1994&JuvenileRheumatic$yob<=2012)
tmp2 <- subset (tmp1, !(tmp1$regenddate < "2012-01-01"))
prev2012<-sum(tmp2$jSLE==1)/nrow(tmp2)*100000
prev2012 
ci <- pois.exact(x=(sum(tmp2$jSLE==1)),pt=(nrow(tmp2)),conf.level = 0.95)
lower95 <- ci$lower * 100000
upper95 <- ci$upper* 100000
lower95
upper95

results <- rbind(results, data.frame(Year = 2012, Cases = sum(tmp2$jSLE==1),Population = nrow(tmp2),Prevalence=prev2012, Lower95CI=lower95, Upper95CI=upper95))

#Prevalence 2011
tmp1<-subset(JuvenileRheumatic,JuvenileRheumatic$regstartdate<="31/12/2011"&JuvenileRheumatic$yob>=1993&JuvenileRheumatic$yob<=2011)
tmp2 <- subset (tmp1, !(tmp1$regenddate < "2011-01-01"))
prev2011<-sum(tmp2$jSLE==1)/nrow(tmp2)*100000
prev2011 
ci <- pois.exact(x=(sum(tmp2$jSLE==1)),pt=(nrow(tmp2)),conf.level = 0.95)
lower95 <- ci$lower * 100000
upper95 <- ci$upper* 100000
lower95
upper95

results <- rbind(results, data.frame(Year = 2011, Cases = sum(tmp2$jSLE==1),Population = nrow(tmp2),Prevalence=prev2011, Lower95CI=lower95, Upper95CI=upper95))

#Prevalence 2010
tmp1<-subset(JuvenileRheumatic,JuvenileRheumatic$regstartdate<="31/12/2010"&JuvenileRheumatic$yob>=1992&JuvenileRheumatic$yob<=2010)
tmp2 <- subset (tmp1, !(tmp1$regenddate < "2010-01-01"))
prev2010<-sum(tmp2$jSLE==1)/nrow(tmp2)*100000
prev2010 
ci <- pois.exact(x=(sum(tmp2$jSLE==1)),pt=(nrow(tmp2)),conf.level = 0.95)
lower95 <- ci$lower * 100000
upper95 <- ci$upper* 100000
lower95
upper95

results <- rbind(results, data.frame(Year = 2010, Cases = sum(tmp2$jSLE==1),Population = nrow(tmp2),Prevalence=prev2010, Lower95CI=lower95, Upper95CI=upper95))

#Prevalence 2009
tmp1<-subset(JuvenileRheumatic,JuvenileRheumatic$regstartdate<="31/12/2009"&JuvenileRheumatic$yob>=1991&JuvenileRheumatic$yob<=2009)
tmp2 <- subset (tmp1, !(tmp1$regenddate < "2009-01-01"))
prev2009<-sum(tmp2$jSLE==1)/nrow(tmp2)*100000
prev2009 
ci <- pois.exact(x=(sum(tmp2$jSLE==1)),pt=(nrow(tmp2)),conf.level = 0.95)
lower95 <- ci$lower * 100000
upper95 <- ci$upper* 100000
lower95
upper95

results <- rbind(results, data.frame(Year = 2009, Cases = sum(tmp2$jSLE==1),Population = nrow(tmp2),Prevalence=prev2009, Lower95CI=lower95, Upper95CI=upper95))

#Prevalence 2008
tmp1<-subset(JuvenileRheumatic,JuvenileRheumatic$regstartdate<="31/12/2008"&JuvenileRheumatic$yob>=1990&JuvenileRheumatic$yob<=2008)
tmp2 <- subset (tmp1, !(tmp1$regenddate < "2008-01-01"))
prev2008<-sum(tmp2$jSLE==1)/nrow(tmp2)*100000
prev2008 
ci <- pois.exact(x=(sum(tmp2$jSLE==1)),pt=(nrow(tmp2)),conf.level = 0.95)
lower95 <- ci$lower * 100000
upper95 <- ci$upper* 100000
lower95
upper95

results <- rbind(results, data.frame(Year = 2008, Cases = sum(tmp2$jSLE==1),Population = nrow(tmp2),Prevalence=prev2008, Lower95CI=lower95, Upper95CI=upper95))

#Prevalence 2007
tmp1<-subset(JuvenileRheumatic,JuvenileRheumatic$regstartdate<="31/12/2007"&JuvenileRheumatic$yob>=1989&JuvenileRheumatic$yob<=2007)
tmp2 <- subset (tmp1, !(tmp1$regenddate < "2007-01-01"))
prev2007<-sum(tmp2$jSLE==1)/nrow(tmp2)*100000
prev2007
ci <- pois.exact(x=(sum(tmp2$jSLE==1)),pt=(nrow(tmp2)),conf.level = 0.95)
lower95 <- ci$lower * 100000
upper95 <- ci$upper* 100000
lower95
upper95

results <- rbind(results, data.frame(Year = 2007, Cases = sum(tmp2$jSLE==1),Population = nrow(tmp2),Prevalence=prev2007, Lower95CI=lower95, Upper95CI=upper95))

#Prevalence 2006
tmp1<-subset(JuvenileRheumatic,JuvenileRheumatic$regstartdate<="31/12/2006"&JuvenileRheumatic$yob>=1988&JuvenileRheumatic$yob<=2006)
tmp2 <- subset (tmp1, !(tmp1$regenddate < "2006-01-01"))
prev2006<-sum(tmp2$jSLE==1)/nrow(tmp2)*100000
prev2006 
ci <- pois.exact(x=(sum(tmp2$jSLE==1)),pt=(nrow(tmp2)),conf.level = 0.95)
lower95 <- ci$lower * 100000
upper95 <- ci$upper* 100000
lower95
upper95

results <- rbind(results, data.frame(Year = 2006, Cases = sum(tmp2$jSLE==1),Population = nrow(tmp2),Prevalence=prev2006, Lower95CI=lower95, Upper95CI=upper95))

#Prevalence 2005
tmp1<-subset(JuvenileRheumatic,JuvenileRheumatic$regstartdate<="31/12/2005"&JuvenileRheumatic$yob>=1987&JuvenileRheumatic$yob<=2005)
tmp2 <- subset (tmp1, !(tmp1$regenddate < "2005-01-01"))
prev2005<-sum(tmp2$jSLE==1)/nrow(tmp2)*100000
prev2005 
ci <- pois.exact(x=(sum(tmp2$jSLE==1)),pt=(nrow(tmp2)),conf.level = 0.95)
lower95 <- ci$lower * 100000
upper95 <- ci$upper* 100000
lower95
upper95

results <- rbind(results, data.frame(Year = 2005, Cases = sum(tmp2$jSLE==1),Population = nrow(tmp2),Prevalence=prev2005, Lower95CI=lower95, Upper95CI=upper95))

#Prevalence 2004
tmp1<-subset(JuvenileRheumatic,JuvenileRheumatic$regstartdate<="31/12/2004"&JuvenileRheumatic$yob>=1986&JuvenileRheumatic$yob<=2004)
tmp2 <- subset (tmp1, !(tmp1$regenddate < "2004-01-01"))
prev2004<-sum(tmp2$jSLE==1)/nrow(tmp2)*100000
prev2004 
ci <- pois.exact(x=(sum(tmp2$jSLE==1)),pt=(nrow(tmp2)),conf.level = 0.95)
lower95 <- ci$lower * 100000
upper95 <- ci$upper* 100000
lower95
upper95

results <- rbind(results, data.frame(Year = 2004, Cases = sum(tmp2$jSLE==1),Population = nrow(tmp2),Prevalence=prev2004, Lower95CI=lower95, Upper95CI=upper95))

#Prevalence 2003
tmp1<-subset(JuvenileRheumatic,JuvenileRheumatic$regstartdate<="31/12/2003"&JuvenileRheumatic$yob>=1985&JuvenileRheumatic$yob<=2003)
tmp2 <- subset (tmp1, !(tmp1$regenddate < "2003-01-01"))
prev2003<-sum(tmp2$jSLE==1)/nrow(tmp2)*100000
prev2003 
ci <- pois.exact(x=(sum(tmp2$jSLE==1)),pt=(nrow(tmp2)),conf.level = 0.95)
lower95 <- ci$lower * 100000
upper95 <- ci$upper* 100000
lower95
upper95

results <- rbind(results, data.frame(Year = 2003, Cases = sum(tmp2$jSLE==1),Population = nrow(tmp2),Prevalence=prev2003, Lower95CI=lower95, Upper95CI=upper95))

#Prevalence 2002
tmp1<-subset(JuvenileRheumatic,JuvenileRheumatic$regstartdate<="31/12/2002"&JuvenileRheumatic$yob>=1984&JuvenileRheumatic$yob<=2002)
tmp2 <- subset (tmp1, !(tmp1$regenddate < "2002-01-01"))
prev2002<-sum(tmp2$jSLE==1)/nrow(tmp2)*100000
prev2002 
ci <- pois.exact(x=(sum(tmp2$jSLE==1)),pt=(nrow(tmp2)),conf.level = 0.95)
lower95 <- ci$lower * 100000
upper95 <- ci$upper* 100000
lower95
upper95

results <- rbind(results, data.frame(Year = 2002, Cases = sum(tmp2$jSLE==1),Population = nrow(tmp2),Prevalence=prev2002, Lower95CI=lower95, Upper95CI=upper95))

#Prevalence 2001
tmp1<-subset(JuvenileRheumatic,JuvenileRheumatic$regstartdate<="31/12/2001"&JuvenileRheumatic$yob>=1983&JuvenileRheumatic$yob<=2001)
tmp2 <- subset (tmp1, !(tmp1$regenddate < "2001-01-01"))
prev2001<-sum(tmp2$jSLE==1)/nrow(tmp2)*100000
prev2001 
ci <- pois.exact(x=(sum(tmp2$jSLE==1)),pt=(nrow(tmp2)),conf.level = 0.95)
lower95 <- ci$lower * 100000
upper95 <- ci$upper* 100000
lower95
upper95

results <- rbind(results, data.frame(Year = 2001, Cases = sum(tmp2$jSLE==1),Population = nrow(tmp2),Prevalence=prev2001, Lower95CI=lower95, Upper95CI=upper95))


#Prevalence 2000
tmp1<-subset(JuvenileRheumatic,JuvenileRheumatic$regstartdate<="31/12/2000"&JuvenileRheumatic$yob>=1982&JuvenileRheumatic$yob<=2000)
tmp2 <- subset (tmp1, !(tmp1$regenddate < "2000-01-01"))
prev2000<-sum(tmp2$jSLE==1)/nrow(tmp2)*100000
prev2000 
ci <- pois.exact(x=(sum(tmp2$jSLE==1)),pt=(nrow(tmp2)),conf.level = 0.95)
lower95 <- ci$lower * 100000
upper95 <- ci$upper* 100000
lower95
upper95

results <- rbind(results, data.frame(Year = 2000, Cases = sum(tmp2$jSLE==1),Population = nrow(tmp2),Prevalence=prev2000, Lower95CI=lower95, Upper95CI=upper95))

#Prevalence 1999
tmp1<-subset(JuvenileRheumatic,JuvenileRheumatic$regstartdate<="31/12/1999"&JuvenileRheumatic$yob>=1981&JuvenileRheumatic$yob<=1999)
tmp2 <- subset (tmp1, !(tmp1$regenddate < "1999-01-01"))
prev1999<-sum(tmp2$jSLE==1)/nrow(tmp2)*100000
prev1999 
ci <- pois.exact(x=(sum(tmp2$jSLE==1)),pt=(nrow(tmp2)),conf.level = 0.95)
lower95 <- ci$lower * 100000
upper95 <- ci$upper* 100000
lower95
upper95

results <- rbind(results, data.frame(Year = 1999, Cases = sum(tmp2$jSLE==1),Population = nrow(tmp2),Prevalence=prev1999, Lower95CI=lower95, Upper95CI=upper95))

#Prevalence 1998
tmp1<-subset(JuvenileRheumatic,JuvenileRheumatic$regstartdate<="31/12/1998"&JuvenileRheumatic$yob>=1980&JuvenileRheumatic$yob<=1998)
tmp2 <- subset (tmp1, !(tmp1$regenddate < "1998-01-01"))
prev1998<-sum(tmp2$jSLE==1)/nrow(tmp2)*100000
prev1998
ci <- pois.exact(x=(sum(tmp2$jSLE==1)),pt=(nrow(tmp2)),conf.level = 0.95)
lower95 <- ci$lower * 100000
upper95 <- ci$upper* 100000
lower95
upper95

results <- rbind(results, data.frame(Year = 1998, Cases = sum(tmp2$jSLE==1),Population = nrow(tmp2),Prevalence=prev1998, Lower95CI=lower95, Upper95CI=upper95))

#Prevalence 1997
tmp1<-subset(JuvenileRheumatic,JuvenileRheumatic$regstartdate<="31/12/1997"&JuvenileRheumatic$yob>=1979&JuvenileRheumatic$yob<=1997)
tmp2 <- subset (tmp1, !(tmp1$regenddate < "1997-01-01"))
prev1997<-sum(tmp2$jSLE==1)/nrow(tmp2)*100000
prev1997 
ci <- pois.exact(x=(sum(tmp2$jSLE==1)),pt=(nrow(tmp2)),conf.level = 0.95)
lower95 <- ci$lower * 100000
upper95 <- ci$upper* 100000
lower95
upper95

results <- rbind(results, data.frame(Year = 1997, Cases = sum(tmp2$jSLE==1),Population = nrow(tmp2),Prevalence=prev1997, Lower95CI=lower95, Upper95CI=upper95))

#Prevalence 1996
tmp1<-subset(JuvenileRheumatic,JuvenileRheumatic$regstartdate<="31/12/1996"&JuvenileRheumatic$yob>=1978&JuvenileRheumatic$yob<=1996)
tmp2 <- subset (tmp1, !(tmp1$regenddate < "1996-01-01"))
prev1996<-sum(tmp2$jSLE==1)/nrow(tmp2)*100000
prev1996 
ci <- pois.exact(x=(sum(tmp2$jSLE==1)),pt=(nrow(tmp2)),conf.level = 0.95)
lower95 <- ci$lower * 100000
upper95 <- ci$upper* 100000
lower95
upper95

results <- rbind(results, data.frame(Year =1996, Cases = sum(tmp2$jSLE==1),Population = nrow(tmp2),Prevalence=prev1996, Lower95CI=lower95, Upper95CI=upper95))

#Prevalence 1995
tmp1<-subset(JuvenileRheumatic,JuvenileRheumatic$regstartdate<="31/12/1995"&JuvenileRheumatic$yob>=1977&JuvenileRheumatic$yob<=1995)
tmp2 <- subset (tmp1, !(tmp1$regenddate < "1995-01-01"))
prev1995<-sum(tmp2$jSLE==1)/nrow(tmp2)*100000
prev1995 
ci <- pois.exact(x=(sum(tmp2$jSLE==1)),pt=(nrow(tmp2)),conf.level = 0.95)
lower95 <- ci$lower * 100000
upper95 <- ci$upper* 100000
lower95
upper95

results <- rbind(results, data.frame(Year = 1995, Cases = sum(tmp2$jSLE==1),Population = nrow(tmp2),Prevalence=prev1995, Lower95CI=lower95, Upper95CI=upper95))

#Plot the trend between year and prevalence
plot(results$Year, results$Prevalence)

#Save results into xlsx files - depending on juvenile rheumatic disease
write_xlsx(results, path = "Point prevalence for jSLE.xlsx")