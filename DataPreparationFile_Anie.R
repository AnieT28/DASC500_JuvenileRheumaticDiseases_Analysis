#===========================DISSERTATION========================================

#Set working directory
setwd("R:/Anie")
getwd()

#Load in data (rdata file)
#Made a copy to not distort the original data
load("DataforAnie - Copy.rdata")

#Summary 
summary(JuvenileRheumatic)

#First 5 rows of data
head(JuvenileRheumatic)

#Structure of data
str(JuvenileRheumatic)

#Check class of data
class(JuvenileRheumatic)

#Convert data to appropriate data types 
JuvenileRheumatic$gender <- factor(JuvenileRheumatic$gender) #sex 
JuvenileRheumatic$yob <- as.numeric(JuvenileRheumatic$yob) #year of birth
JuvenileRheumatic$emis_ddate <- as.Date(JuvenileRheumatic$emis_ddate, format="%d/%m/%Y") #date of death - if present, last day of follow up 
JuvenileRheumatic$regstartdate <- as.Date(JuvenileRheumatic$regstartdate, format = "%d/%m/%Y") #registration to GP
JuvenileRheumatic$regenddate <- as.Date(JuvenileRheumatic$regenddate, format = "%d/%m/%Y") #end of registration to gp
JuvenileRheumatic$acceptable <- factor(JuvenileRheumatic$acceptable)
JuvenileRheumatic$cprd_ddate <- as.Date(JuvenileRheumatic$cprd_ddate, format = "%d/%m/%Y")# Another death date
JuvenileRheumatic$region <- factor(JuvenileRheumatic$region) #region of england
JuvenileRheumatic$lcd <- as.Date(JuvenileRheumatic$lcd, format = "%d/%m/%Y") #last date the GP data sent data to CPRD (use if regend is blank)
JuvenileRheumatic$patienttypeid <- factor(JuvenileRheumatic$patienttypeid)
JuvenileRheumatic$mob <- as.numeric(JuvenileRheumatic$mob)

#Data check for correct conversion
summary(JuvenileRheumatic$emis_ddate) #first record was 14th person at 29/01/1998 so expect 1998-01-29
summary(JuvenileRheumatic$regstartdate) #first date was 25/04/2016 so expect 2016-04-25
summary(JuvenileRheumatic$lcd)
#Updated summary
summary(JuvenileRheumatic$regenddate)

#Recode region, so that it is just region name and not factor number.
# Use dplyr library
library(dplyr)

#Recode region
JuvenileRheumatic <- JuvenileRheumatic %>% 
  mutate(region = recode(
    region, `0` = "None",`1`= "North East", `2`= "North West", `3`= "Yorkshire and the Humber", 
    `4` = "East Midlands",  `5` = "West Midlands",`6` = "East of England", `7` = "London",  `8` = "South East", 
    `9` = "South West", `10` = "Wales", `11` = "Scotland", `12` = "Northern Ireland"))

# Summary of region
summary(JuvenileRheumatic$region)

# Find the last date of observation for each patient
# Found by looking for the earliest date out of the following in the code using pmin()
JuvenileRheumatic <- JuvenileRheumatic %>%
  mutate(lastobsdate = pmin(emis_ddate, cprd_ddate, regenddate,FirstJIADate, FirstJDMDate, FirstjSLEDate, lcd, na.rm = TRUE))

#Print summary of lastobsdate
summary(JuvenileRheumatic$lastobsdate)

#Within the dataset, there are 17,933 individuals with JIA, 379 with JDM and 848 with jSLE

#Make smaller dataset with relevant columns - double check later
myDT<-JuvenileRheumatic[,c("patid","pracid","gender","yob","mob","regstartdate",
                           "region","JIA","JDM","jSLE","lastobsdate",
                           "FirstJIADate", "FirstJDMDate", "FirstjSLEDate")]

# Includes relevant columns to work on 
save(myDT, file="cleaned_data.RData")

#====================DATA CLEANING = GENDER=====================================
#Summary of gender - ignoring 1701 I gender
summary(data$gender) #24197245 F and 22596943 M

#Removing I from gender - not in analysis
data <- data[gender %in% c("M", "F")]

#Drop I as an level - droplevels() drops levels without use 
data$gender <- droplevels(data$gender)

#New summary of gender
summary(data$gender)
#==================== IMPUTE DATE AND MONTH OF BIRTH ===========================
#Number of people with missing mob
sum(is.na(data$mob)) #42,319,405 individuals

#Set random seed
set.seed(42)

#Make new column with month of birth
data[, birthmonth := mob]

#Use uniform distribution to randomly assign mob
data[is.na(birthmonth),birthmonth := sample(x=1:12,size=.N,replace = TRUE)]

#Date of birth assumed to be 15th day of each month - no misclassification was greater than 3 weeks
data[,birthdate := as.Date(sprintf("%04d-%02d-15",yob,birthmonth))]

#Summary of date of birth
summary(data$birthdate)
#===============================================================================
#====================== AGE AT DIAGNOSIS =======================================
#Age at diagnosis for JIA
data[,JIA_age := floor(time_length(interval(birthdate,FirstJIADate),unit="years"))]

#Age at diagnosis for JDM 
data[,JDM_age := floor(time_length(interval(birthdate,FirstJDMDate),unit="years"))]

#Age at diagnosis for jSLE
data[,jSLE_age := floor(time_length(interval(birthdate,FirstjSLEDate),unit="years"))]

#Approximate 19th birthday - for jSLE and JDM incidence
data[, Age19Date := birthdate %m+% years(19)]

#Approximate 17th birthday - for JIA incidence
data[, Age17Date := birthdate %m+% years(17)]

#======================== CHECK AGE AT DIAGNOSIS ===============================
#Sum of the age of diagnosis between 0-18 = 16,870 individuals with JIA = CHANGED NOW
sum(!is.na(data$JIA_age) & data$JIA_age >= 0 & data$JIA_age <= 18)

#Sum of age of diagnosis between 0-18 = 369 individuals with JDM
sum(!is.na(data$JDM_age) & data$JDM_age >= 0 & data$JDM_age <= 18)

#Sum of age of diagnosis between 0-18 = 758 individuals with JDM
sum(!is.na(data$jSLE_age) & data$jSLE_age >= 0 & data$jSLE_age <= 18)

#Sum of counts with diagnosis age less than 0
sum(data$JIA_age < 0, na.rm = TRUE) #315 negative ages - data quality issues
sum(data$JDM_age < 0, na.rm = TRUE) #6 negative ages
sum(data$jSLE_age < 0, na.rm = TRUE) #19 negative ages

#==========================INVESTIGATIONS=======================================
#Check regstartdate, lastobsdate 
#These have 1899-12-30 as a first JIA date = ERROR!
data[JIA_age < 0, c("yob", "mob", "FirstJIADate", "regstartdate", "lastobsdate")]

#Summary of FirstJIADate
summary(data$FirstJIADate)

#Min diagnosis of JIA is 1899-12-30 = 269 with this
summary(data$FirstJIADate == "1899-12-30")

#Histogram of first JIA diagnosis ages
#Some with negative ages 
hist(data$JIA_age)