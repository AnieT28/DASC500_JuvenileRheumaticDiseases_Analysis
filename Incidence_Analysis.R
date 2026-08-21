#Set working directory
setwd("R:/Anie")
getwd()

#Load in data
load("cleaned_data_v2.RData") #Load in for JDM and jSLE analysis
load("R:/Anie/cleaned_data_v2_JIA.RData") #Load in for JIA incidence

#Important libraries
library(data.table)
library(lubridate)
#install.packages("epitools")
library(epitools)
library(ggplot2)
library(writexl)
library(readxl)
library(dplyr)
library(gtsummary)
library(gt)

#NOTE:
#For all the functions, need to make sure the age limits match for the specific disease analysis. e.g. JIA = 0-16
#Also, when running the function, make sure the Date of diagnosis and Diagnosis age match for the disease e.g. FirstJIADate & JIA_age

#===============================================================================
# Sum of cases in total excluding registration start date restriction
sum(data$JIA == 1 & data$JIA_age >= 0 & data$JIA_age <= 16, na.rm = TRUE) #13309
sum(data$JDM == 1 & data$JDM_age >= 0 & data$JDM_age <= 18, na.rm = TRUE) #369
sum(data$jSLE == 1 & data$jSLE_age >= 0 & data$jSLE_age <= 18, na.rm = TRUE) #758
#===============================================================================
#Data checks
summary(data)
str(data)
head(data)

#Seeing earliest to latest registration date 
summary(data$regstartdate)

#=============DETERMINING STUDY PERIOD DUE TO DATA QUALITY ISSUES===============
#Plot histogram by year
#Extract year from regstartdate
data[, start_year := as.integer(format(regstartdate, "%Y"))]

#Sum all counts all by start_year
year_counts <- data[, .N, by = start_year]

#Plot histogram using geom_hist()
ggplot(year_counts, aes(start_year, N)) +
  geom_col(fill = "steelblue", colour = "black",linewidth = 0.5) +
  scale_x_continuous(breaks = seq(min(year_counts$start_year),
                 max(year_counts$start_year), by = 10)) +
  labs(x = "Year", y = "Number of patients") +
  theme_minimal()

#I have determined to do from the year 1995 to 2023, using evidence from histogram
#and the CPRD Aurum data profile.

#================= ANNUAL INCIDENCE CALCULATION ================================
annual_incidence <- function(data, diagnosis_col, age_col, group_var = NULL,
                             start_year = 1995, end_year = 2023, multiplier = 100000) {
  
  diagnosis <- data[[diagnosis_col]]
  diagnosis_age <- data[[age_col]]
  
  CYP <- diagnosis_age >= 0 & diagnosis_age <= 16 #changed to 16 for JIA only
  
  #All years 1995-2023
  years <- start_year:end_year
  
  # If no grouping requested
  if (is.null(group_var)) {
    groups <- "Overall"
    group_values <- rep("Overall", nrow(data))
  } else {
    #If group requested, groups is equal to unique variable in groups
    groups <- unique(data[[group_var]])
    groups <- groups[!is.na(groups)]
  }
  
  #Empty list to store results for each group
  #Each element of list becomes 1 data table
  all_results <- list()
  
  #Loop through each group
  for(g in groups){
    #Overall 
    if(is.null(group_var)){
      keep <- rep(TRUE, nrow(data))
    } else {
      #g is the variable in group
      keep <- data[[group_var]] == g
    }
    
    #Empty results table for the group
    results <- data.table(Group = g, Year = years, Cases = integer(length(years)),
      PersonYears = numeric(length(years)), IncidenceRate = numeric(length(years)),
      Lower95CI = numeric(length(years)),Upper95CI = numeric(length(years))
    )
    
    #Loop through calendar year
    for(i in seq_along(years)){
      yr <- years[i]
      #Specify beginning and end of calendar year
      year_start <- as.Date(sprintf("%d-01-01", yr))
      year_end   <- as.Date(sprintf("%d-12-31", yr))
      
      #Numerator
      cases <- sum(keep & !is.na(diagnosis) & diagnosis >= year_start & diagnosis <= year_end & CYP)
      
      #Follow-up 
      entry <- pmax(data$regstartdate, year_start)
      exit <- pmin(diagnosis, data$lastobsdate, data$Age17Date, year_end, na.rm = TRUE)
      
      #Number of days contributed during this calendar year
      days <- pmax(0, as.numeric(exit - entry) + 1)
      
      #Denominator
      py <- sum(days[keep],na.rm = TRUE) / 365.25
      
      #Poisson CI
      if(is.na(py) || py <=0){
        IR <- NA
        Lower <- NA
        Upper <- NA
      } else{
        ci <- pois.exact(cases, pt = py)
        IR <- ci$rate * multiplier
        Lower <- ci$lower * multiplier
        Upper <- ci$upper * multiplier
      }
      
      #Store results into empty data table
      results[i, `:=`(Cases = cases, PersonYears = py, IncidenceRate = (cases / py) * multiplier,
              Lower95CI = Lower, Upper95CI = Upper)]
    }
    #Save group results into the list
    all_results[[as.character(g)]] <- results
  }
  #Combine every group's results into one data table
  rbindlist(all_results)
}

#============================ RESULTS ==========================================
#Calculate incidence rate of JIA annually - overall, sex, region
JIA_annual <- annual_incidence(data, diagnosis_col = "FirstJIADate", age_col = "JIA_age")
JIA_annual_sex <-annual_incidence(data, diagnosis_col = "FirstJIADate", age_col = "JIA_age", group_var = "gender")
JIA_annual_region <-annual_incidence(data, diagnosis_col = "FirstJIADate", age_col = "JIA_age", group_var = "region")

write_xlsx(list("JIA Overall Incidence" = JIA_annual,
                "JIA Sex-specific Incidence" = JIA_annual_sex,
                "JIA region-specific Incidence" = JIA_annual_region),
           path = "Annual JIA Incidence rate - overall,sex,region.xlsx")

#Calculate annual incidence rate for JDM - overall, sex, region
JDM_annual <- annual_incidence(data, diagnosis_col = "FirstJDMDate", age_col = "JDM_age")
JDM_annual_sex <- annual_incidence(data, diagnosis_col = "FirstJDMDate", age_col = "JDM_age", group_var = "gender")
JDM_annual_region <- annual_incidence(data, diagnosis_col = "FirstJDMDate", age_col = "JDM_age", group_var = "region")

#Calculate annual incidence rates for jSLE - overall, sex, region
jSLE_annual <- annual_incidence(data, diagnosis_col = "FirstjSLEDate", age_col = "jSLE_age")
jSLE_annual_sex <- annual_incidence(data, diagnosis_col = "FirstjSLEDate", age_col = "jSLE_age", group_var = "gender")
JSLE_annual_region <- annual_incidence(data, diagnosis_col = "FirstjSLEDate", age_col = "jSLE_age", group_var = "region")

#Save results to xlsx
write_xlsx(JDM_annual_sex, path = "JDM annual sex incidence rates.xlsx")
write_xlsx(jSLE_annual_sex, path = "jSLE annual sex incidence rates.xlsx")
write_xlsx(JDM_annual_region, path = "JDM annual region incidence rates.xlsx")
write_xlsx(JIA_annual_region, path = "JIA annual region incidence rates.xlsx")
write_xlsx(JSLE_annual_region, path = "jSLE annual region incidence rates.xlsx")

#===================ANNUAL AGE-SPECIFIC INCIDENCE RATES=========================
annual_age_incidence <- function(data,diagnosis_col,age_col,age_start = 0,age_end = 18, #edited for JIA 
                                 start_year = 1995,end_year = 2023,multiplier = 100000){
 
  #Diagnosis specification and age of diagnosis
  diagnosis <- data[[diagnosis_col]]
  diagnosis_age <- data[[age_col]]
  
  #Study period 
  study_start <- as.Date(sprintf("%d-01-01", start_year))
  study_end   <- as.Date(sprintf("%d-12-31", end_year))
  
  #Pre-compute study entry and exit dates ONCE
  diagnosis_exit <- diagnosis
  diagnosis_exit[is.na(diagnosis_exit)] <- study_end
  
  #Study entry is the later of:
  study_entry <- pmax(data$regstartdate, study_start)
  
  #Study exit is the earliest of:
  study_exit <- pmin(diagnosis_exit, data$lastobsdate, data$Age19Date, study_end)

  #Number of years and ages
  n_years <- end_year - start_year + 1
  n_ages  <- age_end - age_start + 1
  
  # Total number of rows
  n_rows <- n_years * n_ages
  
  #Empty results table
  results <- data.table(Year = integer(n_rows), Age = integer(n_rows),
    Cases = integer(n_rows), PersonYears = numeric(n_rows), IncidenceRate = numeric(n_rows),
    Lower95CI = numeric(n_rows), Upper95CI = numeric(n_rows))
  
  #Row counter
  row <- 1
  
  #Loop through calendar years
  for(year in start_year:end_year){
    
    #Define the beginning and end of the current year
    year_start <- as.Date(sprintf("%d-01-01", year))
    year_end   <- as.Date(sprintf("%d-12-31", year))
    
    #Identify individuals contributing person-time
    keep <- study_entry <= year_end & study_exit >= year_start
    
    #Skip this year if nobody contributes follow-up
    
    if(!any(keep))
      next
    
    #Subset only the individuals contributing person-time during this year.
    diagnosis_year <- diagnosis[keep]
    diagnosis_age_year <- diagnosis_age[keep]
    
    study_entry_year <- study_entry[keep]
    study_exit_year <- study_exit[keep]
    
    birth_year <- data$birthdate[keep]
    
    #Loop through ages
    for(age in age_start:age_end){
      
      #Calculate the birthday interval corresponding to the current age.
      age_start_date <- birth_year %m+% years(age)
      age_end_date   <- birth_year %m+% years(age + 1)
      
      #Numerator
      cases <- sum(!is.na(diagnosis_year) &
          diagnosis_year >= year_start & diagnosis_year <= year_end & diagnosis_age_year == age)
      
      #Denominator
      entry <- pmax(study_entry_year, age_start_date, year_start)
      exit <- pmin(study_exit_year, age_end_date, year_end)
      
      #Follow-up in days
      #Negative values indicate no overlap and are set to zero.
      days <- pmax(0, as.numeric(exit - entry) + 1)
      
      #Convert follow-up from days to person-years
      py <- sum(days) / 365.25

      #Incidence rate and 95% confidence interval
      if(cases == 0 || py <= 0){
        IR <- 0
        lower <- NA
        upper <- NA
        
      }else{
        #Exact Poisson confidence intervals
        ci <- pois.exact(cases, pt = py)
        IR <- ci$rate * multiplier
        lower <- ci$lower * multiplier
        upper <- ci$upper * multiplier
      }
      
      #Store results
      results[row, `:=`(Year = year,Age = age,Cases = cases, PersonYears = py,
                        IncidenceRate = IR, Lower95CI = lower, Upper95CI = upper
      )]
      row <- row + 1
    }
}
  return(results)
}

#JIA age-specific incidence rates
JIA_annual_age <- annual_age_incidence(data = data, diagnosis_col = "FirstJIADate", age_col = "JIA_age")

#MAKE SURE TO EDIT AGE 16 → 18 FOR JSLE INCIDENCE ANALYSIS
JDM_annual_age <- annual_age_incidence(data = data, diagnosis_col = "FirstJDMDate", age_col = "JDM_age")
jSLE_annual_age <- annual_age_incidence(data = data, diagnosis_col = "FirstjSLEDate", age_col = "jSLE_age")


#Write up as xlsx file
write_xlsx(JIA_annual_age, path = "JIA annual age-specific incidence rates.xlsx")
write_xlsx(jSLE_annual_age, "R:/Anie/Incidence xlsx sheets/jSLE/jSLE_annual_age-specific.xlsx") #to run

#========== OVERALL AND BY SEX AND REGION INCIDENCE RATE CALCULATIONS ==========
incidence_rates_sr <- function(data,diagnosis_col,age_col,group_var = NULL,start_year = 1995,
                              end_year = 2023,multiplier = 100000){
  #Diagnosis column
  diagnosis <- data[[diagnosis_col]]
  #Age at diagnosis
  diagnosis_age <- data[[age_col]]
  #Childhood cases only
  CYP <- diagnosis_age >= 0 & diagnosis_age <= 16
  #Study period
  study_start <- as.Date(sprintf("%d-01-01", start_year))
  study_end   <- as.Date(sprintf("%d-12-31", end_year))
  #Entry and exit dates
  entry <- pmax(data$regstartdate, study_start)
  exit <- pmin(diagnosis,data$lastobsdate,data$Age19Date,study_end,na.rm = TRUE)
  # Follow-up time (days)
  days <- pmax(0, as.numeric(exit - entry) + 1)

  #OVERALL INCIDENCE 
  if(is.null(group_var)){
    #Number of cases
    cases <- sum(!is.na(diagnosis) &diagnosis >= study_start &diagnosis <= study_end &CYP)
    #Person-years
    py <- sum(days) / 365.25
    #95% CI
    ci <- pois.exact(cases, pt = py)
    
    #Return data table with results
    return(
      data.table(Group = "Overall",Level = "Overall",Cases = cases,PersonYears = py,
        IncidenceRate = ci$rate * multiplier,Lower95CI = ci$lower * multiplier,
        Upper95CI = ci$upper * multiplier))
  }
  #GROUPED RESULTS (Sex or Region - can specify)
  #Found from factored variables
  groups <- sort(unique(na.omit(data[[group_var]])))
  
  results <- lapply(groups, function(g){
    
    #Subset by group
    subset <- data[[group_var]] == g
    #Cases
    cases <- sum(!is.na(diagnosis) & diagnosis >= study_start & diagnosis <= study_end & CYP & subset)
    #Person-years
    py <- sum(days[subset]) / 365.25
    #Check for 0 cases - will return error
    if(cases == 0 || py <= 0){
      IR <- 0
      Lower <- NA
      Upper <- NA
    } else{
      ci <- pois.exact(cases, pt = py)
      IR <- ci$rate * multiplier
      Lower <- ci$lower * multiplier
      Upper <- ci$upper * multiplier
    }
    data.table(Group = group_var, Level = g,Cases = cases, PersonYears = py,
      IncidenceRate = IR, Lower95CI = Lower, Upper95CI = Upper)
  })
   return(rbindlist(results))
}

#Incidence rates of JIA overall, by sex and region
JIA_overall <- incidence_rates_sr(data, diagnosis_col = "FirstJIADate", age_col = "JIA_age")
JIA_sex <- incidence_rates_sr(data, diagnosis_col = "FirstJIADate", age_col = "JIA_age", group_var = "gender")
JIA_region <- incidence_rates_sr(data, diagnosis_col = "FirstJIADate", age_col = "JIA_age", group_var = "region")

#Incidence rates of JDM overall, by sex and region
JDM_overall <- incidence_rates_sr(data, diagnosis_col = "FirstJDMDate", age_col = "JDM_age")
JDM_sex <- incidence_rates_sr(data, diagnosis_col = "FirstJDMDate", age_col = "JDM_age", group_var = "gender")
JDM_region <- incidence_rates_sr(data, diagnosis_col = "FirstJDMDate", age_col = "JDM_age", group_var = "region")

#Incidence rates of jSLE overall, by sex and region
jSLE_overall <- incidence_rates_sr(data, diagnosis_col = "FirstjSLEDate", age_col = "jSLE_age")
jSLE_sex <- incidence_rates_sr(data, diagnosis_col = "FirstjSLEDate", age_col = "jSLE_age", group_var = "gender")
jSLE_region <- incidence_rates_sr(data, diagnosis_col = "FirstjSLEDate", age_col = "jSLE_age", group_var = "region")

#Create data tables with all results
JIA_results  <- rbindlist(list(JIA_overall, JIA_sex, JIA_region))
JDM_results  <- rbindlist(list(JDM_overall, JDM_sex, JDM_region))
jSLE_results <- rbindlist(list(jSLE_overall, jSLE_sex, jSLE_region))

#Save data tables as csv files = All have been saved
write_xlsx(JIA_results, path = "JIA overall, sex and region.xlsx")
write_xlsx(JDM_results, path = "JDM overall, sex and region.xlsx")
write_xlsx(jSLE_results, path = "jSLE overall, sex and region.xlsx")

#=========FUNCTION TO CALCULATE INCIDENCE RATE BY INDIVIDUAL AGES 0-18==========

#Want to check whether anyone has a diagnosis at 0
#There are some cases at 0
sum(data$JIA_age == 0, na.rm = TRUE) # 183 cases
sum(data$JDM_age == 0, na.rm = TRUE) # 4 cases
sum(data$jSLE_age == 0, na.rm = TRUE) # 15 cases

#Function to calculate age incidence
age_incidence <- function(data, diagnosis_col, age_col, start_year = 1995, end_year = 2023,multiplier = 100000){
  
  #Diagnosis dates and age at diagnosis
  diagnosis <- data[[diagnosis_col]]
  diagnosis_age <- data[[age_col]]
  
  #Study period
  study_start <- as.Date(sprintf("%d-01-01", start_year))
  study_end   <- as.Date(sprintf("%d-12-31", end_year))
  
  #Empty results table
  results <- data.table(Age = 0:18,Cases = integer(19), PersonYears = numeric(19),
    IncidenceRate = numeric(19), Lower95CI = numeric(19), Upper95CI = numeric(19))

  #Loop through ages
  for(age in 0:18){
    #Numerator
    cases <- sum(!is.na(diagnosis) & diagnosis >= study_start &
                   diagnosis <= study_end & diagnosis_age == age)
    
    #Start of this age interval (birthday when individual turns age)
    age_start <- data$birthdate %m+% years(age)
    
    #End of this age interval (birthday when individual turns age + 1)
    age_end <- data$birthdate %m+% years(age + 1)
    
    #Entry date
    entry <- pmax(data$regstartdate, study_start,age_start)
    
    #Exit date
    exit <- pmin(diagnosis, data$lastobsdate, study_end, age_end, data$Age19Date, na.rm = TRUE)
    
    #Person-time while exactly this age
    days <- pmax(0,as.numeric(exit - entry) + 1)
    
    #Person-years
    py <- sum(days) / 365.25
    
    #Incidence rate and Poisson CI
    if(cases == 0 || py <= 0){
      IR <- 0
      lower <- NA
      upper <- NA
      
    } else{
      ci <- pois.exact(cases, pt = py)
      IR <- ci$rate * multiplier
      lower <- ci$lower * multiplier
      upper <- ci$upper * multiplier
    }

    #Store results
    results[Age == age,`:=`(Cases = cases, PersonYears = py, IncidenceRate = IR,
              Lower95CI = lower, Upper95CI = upper)]
  }
  return(results)
}

#Print results
JIA_age_incidence <- age_incidence(data, diagnosis_col = "FirstJIADate", age_col = "JIA_age")
JDM_age_incidence <- age_incidence(data, diagnosis_col = "FirstJDMDate", age_col = "JDM_age")
JSLE_age_incidence <- age_incidence(data, diagnosis_col = "FirstjSLEDate", age_col = "jSLE_age")

#============================SENSITIVITY ANALYSIS===============================

#copy of data w/o NI
data_excl_NI <- data[region != "Northern Ireland"]

#Number of individuals
nrow(data) #46,794,188
nrow(data_excl_NI) #46,676,062

#Calculate difference in population after exclusion of NI
46794188-46676062 #118216 individuals from NI excluded in sensitivity analysis


# Run code before for both with NI and without NI
#Just edit age limits for the appropriate age boundary depending on the disease

#Study period
study_start <- as.Date("1995-01-01")
study_end <- as.Date("2023-12-31")
#XXXX diagnosis and age at diagnosis
diagnosis <- data$FirstJIADate
diagnosis_age <- data$JIA_age
#Childhood cases only (0-16)
CYP <- diagnosis_age >=0 & diagnosis_age <=18
#Entry date
entry <- pmax(data$regstartdate, study_start)
#Exit date
exit<- pmin(diagnosis,data$lastobsdate,data$Age19Date,study_end,na.rm = TRUE)
#Follow-up time in days
days <- pmax(0, as.numeric(exit-entry)+1)
#Number of XXXX cases
cases <-sum(!is.na(diagnosis)&diagnosis>=study_start&diagnosis<=study_end&CYP)
#Person-years
py <- sum(days)/365.25
#Poisson 95% CI
ci <- pois.exact(cases, pt = py)
#Incidence rate per 100,000
incidence_rate <- ci$rate * 100000
lower95 <- ci$lower *100000
upper95 <- ci$upper *100000
#Results
results <- data.table(Cases = cases, PersonYears=py,IncidenceRate=incidence_rate,Lower95CI=lower95,Upper95CI=upper95)
#Print results
results

#Overall incidence of JIA without function = 12653 cases, 18.83 (95% CI: 18.51-19.16)
#Overall incidence of JIA without NI individuals = 12,608 cases, 18.88 (95% CI: 18.55-19.21)

#Overall incidence of JDM with NI individuals = 339, 0.45 (95% CI: 0.41-0.51) 
#Overall incidence of JDM without NI individuals = 339, 0.45 (95% CI: 0.41-0.51)

#Overall incidence of jSLE with NI individuals = 716 cases, 0.96 (95% CI: 0.89-1.03)
#Overall incidence of jSLE without NI individuals = 713 cases, 0.96 (95% CI: 0.89-1.03)

#==================================END OF CODE==================================
