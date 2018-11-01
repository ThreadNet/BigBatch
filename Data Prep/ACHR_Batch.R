##########################################################################################################
# THREADNET:  Batch processing for larger data sets

# (c) 2017 Michigan State University. This software may be used according to the terms provided in the
# GNU General Public License (GPL-3.0) https://opensource.org/licenses/GPL-3.0?
# Absolutely no warranty!
##########################################################################################################
library(tidyr)
library(data.table)
library(dplyr)
library(ThreadNet)
library(ngram)
library(lubridate)
library(expss)

cd = apply_labels(clinic_days,
                  Clinic_date =  'ClinicDate',
                  YMD_date = 'Date',
                  NEvents = 'Num Steps',
                  NumVisits = 'Num Visits',
                  NetComplexity= 'Complexity',
                  CompressRatio = 'Compressibility',
                  Clinic = 'Clinic',
                  NumUniqueProcedures = 'Num Unique Procedures',
                  NumUniqueDiagnosisGroups  = 'Num Unique Diagnosis Groups', 
                  NumPhysicians  = 'Num Physicians', 
                  Weekday  = 'Weekday', 
                  Month  = 'Month')

cd %>% 
  tab_cells(NetComplexity,NumVisits,NEvents/NumVisits,NumUniqueProcedures,NumUniqueDiagnosisGroups,
            NumUniqueProcedures/NumVisits,NumUniqueDiagnosisGroups/NumVisits) %>% 
  tab_cols(total(), Clinic) %>% 
  tab_stat_mean() %>% 
  tab_pivot()

cd %>% 
  tab_cells(Weekday) %>% 
  tab_cols(total(), Clinic) %>% 
  tab_stat_cases() %>% 
  tab_pivot()


v %>% 
  tab_cells( NEvents, VisitDuration, NetComplexity, CompressRatio, Action_count, Workstation_count, Role_count ) %>% 
  tab_cols(total(), Clinic) %>% 
  tab_stat_mean() %>% 
  tab_pivot()

b %>% 
  tab_cells( NEvents, VisitDuration, NetComplexity,  Action_count, Workstation_count, Role_count ) %>% 
  tab_cols(total(), Diagnosis_group) %>% 
  tab_stat_mean() %>% 
  tab_pivot() %>% 
  tab_transpose()

b %>% 
  tab_cells(  NEvents, VisitDuration, NetComplexity,  Action_count ) %>% 
  tab_cols(total(), Diagnosis_group) %>% 
  tab_stat_mean_sd_n() %>% 
  tab_pivot() %>% 
  tab_transpose()

b %>% 
  tab_cells( NEvents, VisitDuration, NetComplexity,   Action_count  ) %>% 
  tab_cols(total(), Physician) %>% 
  tab_stat_mean_sd_n() %>% 
  tab_pivot() %>% 
  tab_transpose()


read_ACHR_data <- function(){



  
  
  # This code is tailored for reading in new data from URMC, October 2018
  # read the file into data frame
  d <<- fread('auditfinal_10022018.csv')
  
  # Sort and convert to data.table
  vt <<-data.table(arrange(d,desc(Visit_ID,asc(Timestamps))))
  
  # move and rename the columns, as needed
  cn=colnames(vt)
  cnx = cn[c(13,1:12,14:15)]
  vt=setcolorder(vt,cnx)
  setnames(vt,'Timestamps','tStamp')
  setnames(vt,'X','seqn')
  
  # make a dataframe copy just for fun
  vdf <<- as.data.frame(vt)
  
  # call the functions that clean up the occurrences
  # Column names are hard coded from the UMRC file
  # EVENT_CF = c('Action','Workstation','Role','Clinic','Diagnosis_Group')
  EVENT_CF = c('Action','Workstation','Role')
  TN = 'Visit_ID'
  
  # This converts numbers to char and replaces spaces with underscore
  occ = cleanOccBatch(vdf)
  
  # ThreadNet code assumes these columns will be there, so we need to add them
  occ$threadNum = occ$Visit_ID
  occ$seqNum = occ$seqn
  
  # because the occurrences are already sorted, and they have sequence numbers, we should be good to go. 
  # new_occ = combineContextFactors(occ,EVENT_CF,newColName(EVENT_CF))[1:100000,]
  new_occ = combineContextFactors(occ,EVENT_CF,newColName(EVENT_CF)) 
  
  
  visits = ACHR_batch_visits(new_occ[1:1000,], TN, EVENT_CF)
  save(visits, file='ACHR_visits_with_timestamps.rData')
  
  clinic_days = ACHR_batch_clinic_days(new_occ[1:1000,], TN, EVENT_CF)
  save(clinic_days, file='ACHR_clinic_days.rData')
  
  
  
}
  
  
###################################################################
###################################################################
ACHR_batch_clinic_days <- function(occ,TN, CFs) {
  
  # Name for column that has events
  DV= newColName(CFs)
  
  # name for column that defines the buckets
  BU = 'Clinic_ymd'
  
    # first  get the date only 
    occ$ymd <- format(as.POSIXct(occ$tStamp),"%Y-%m-%d")
    
    # make a new column for clinic + day
    occ = unite(occ, 'Clinic_ymd', c('Clinic','ymd'),sep='_',remove = 'false')
    

  # pick subsets -- one visit at a time in this version, but could be more
  bucket_list <- unique(occ[['Clinic_ymd']])
  
  # get the size (number of buckets)
  N = length(bucket_list)
  print(N)
  
  
  # pre-allocate the data.table.  Tables are supposed to be faster.
  ACHR = data.table(bucket=integer(N),
                    Clinic_date =  character(N),
                    YMD_date = character(N),
                    NEvents = integer(N),
                    NumVisits = numeric(N),
                    NetComplexity=double(N),
                    CompressRatio = double(N),
                    Clinic = character(N),
                    NumUniqueProcedures = numeric(N),
                    NumUniqueDiagnosisGroups  = numeric(N), 
                    NumPhysicians  = numeric(N), 
                    Weekday  = character(N), 
                    Month  = character(N)
  )
  
  # Now add columns for the IVs.  There will be three for each IV
  
  # Add the IV columns
  for (cf in CFs){
    
    ACHR[, paste0(cf,"_count"):= double(N)]
    ACHR[, paste0(cf,"_compression"):= double(N)]
    ACHR[, paste0(cf,"_entropy"):= double(N)]
    
  }
  
  # loop through the buckets. Result will be data frame with one row per bucket
  for (i in 1:N){
    
    b = i #  as.integer(bucket_list[i])
    
    # print once every 10 buckets
    if (b%%10==0) {print(b)}
    
    # select the threads that go in this bucket
    df = occ[occ[[BU]] ==bucket_list[i],]
    
    # bucket number
    ACHR[b,bucket := b]
    
    # length of the thread (number of rows)
    ACHR[b,NEvents := nrow(df)]
    
    # only do the computations if there are more than two occurrences
    if (nrow(df) > 2) {
      
      # compute the duration of the visit in hours
 #     ACHR[b,VisitDuration := difftime(max(lubridate::ymd_hms(df$tStamp)),  min(lubridate::ymd_hms(df$tStamp)), units='hours') ]
      
      # compressibility of DV
      ACHR[b,CompressRatio := compression_index(df,DV)]
      
      # NetComplexity of DV
      # First get the network
      n = threads_to_network_original(df,TN, DV)
      ACHR[b,NetComplexity := estimate_network_complexity( n )]
      
      # compute stuff on each context factor
      for (cf in CFs){
        
        # Count the unique elements in each cf
        ACHR[b, paste0(cf,"_count") :=  length(unique(df[[cf]])) ]
        
        # get the compression
        ACHR[b, paste0(cf,"_compression") := compression_index(df,cf) ]
        
        # get the entropy
        ACHR[b, paste0(cf,"_entropy") := compute_entropy(table(df[[cf]])[table(df[[cf]])>0]) ]
        
      }
    } # kf nrows > 2
    
    # Now copy in the rest of data  
    # this works because one visit is one bucket
    # count the number of different diagnoses --> typical indicators of complexity
    
    
    ACHR[b,'Clinic_date'  := df[1,'Clinic_ymd']]
    ACHR[b,'YMD_date'  := df[1,'ymd']]
    ACHR[b,'NumVisits'  :=  length(unique(df[['Visit_ID']]))]
    ACHR[b,'Clinic' := df[1,'Clinic']]
    ACHR[b,'NumUniqueProcedures' := length(unique(df[['Proc']]))]
    ACHR[b,'NumUniqueDiagnosisGroups'  := length(unique(df[['Diagnosis_Group']]))]  
    ACHR[b,'NumPhysicians'  := length(unique(df[['Physician']]))]  
    ACHR[b,'Weekday'  := df[1,'Weekday']]
    ACHR[b,'Month'  := df[1,'Month']]
    
    
    
  } # loop thru buckets
  
  # return the table
  return(ACHR)
}

###################################################################
###################################################################



# ACHR stands for Antecedents of Complexity in Healthcare Routines.  
# This is function is set up to compute process parameters on thousands of patient visits.
ACHR_batch_visits <- function(occ,TN, CFs) {

DV= newColName(CFs)

# pick subsets -- one visit at a time in this version, but could be more
bucket_list <- make_buckets_1(occ, 'Visit_ID')

# get the size (number of buckets)
N = length(bucket_list)
print(bucket_list)


# pre-allocate the data.table.  Tables are supposed to be faster.
ACHR = data.table(bucket=integer(N),
                  NEvents = integer(N),
                  VisitStart= character(N),
                  VisitStartInt = integer(N),
                  VisitDuration=double(N),   
                  NetComplexity=double(N),
                  CompressRatio = double(N),
                  Visit_ID  = character(N), 
                  Subject_ID  = character(N), 
                  Clinic = character(N),
                  Procedure = character(N),
                  Diagnosis = character(N),   # diag in the raw data
                  Diagnosis_group  = character(N), 
                  Physician  = character(N), 
                  Weekday  = character(N), 
                  Month  = character(N)
                  )

# Now add columns for the IVs.  There will be three for each IV

# Add the IV columns
for (cf in CFs){

  ACHR[, paste0(cf,"_count"):= double(N)]
  ACHR[, paste0(cf,"_compression"):= double(N)]
  ACHR[, paste0(cf,"_entropy"):= double(N)]

}

# loop through the buckets. Result will be data frame with one row per bucket
for (i in 1:N){

  b = i #  as.integer(bucket_list[i])

  # print once every 100 visits
  if (b%%100==0) {print(b)}
  
  # select the threads that go in this bucket
    df = occ[occ[[TN]] ==bucket_list[i],]

    # bucket number
    ACHR[b,bucket := b]

    # length of the thread (number of rows)
    ACHR[b,NEvents := nrow(df)]

    # only do the computations if there are more than two occurrences
    if (nrow(df) > 2) {

      # compute the duration of the visit in hours
      ACHR[b,VisitDuration := difftime(max(lubridate::ymd_hms(df$tStamp)),  min(lubridate::ymd_hms(df$tStamp)), units='hours') ]
     
    # compressibility of DV
    ACHR[b,CompressRatio := compression_index(df,DV)]

    # NetComplexity of DV
    # First get the network
    n = threads_to_network_original(df,TN, DV)
    ACHR[b,NetComplexity := estimate_network_complexity( n )]

  # compute stuff on each context factor
  for (cf in CFs){

    # Count the unique elements in each cf
    ACHR[b, paste0(cf,"_count") :=  length(unique(df[[cf]])) ]

    # get the compression
    ACHR[b, paste0(cf,"_compression") := compression_index(df,cf) ]

    # get the entropy
    ACHR[b, paste0(cf,"_entropy") := compute_entropy(table(df[[cf]])[table(df[[cf]])>0]) ]

  }
} # kf nrows > 2

    # Now copy in the rest of data  
    # this only works because one visit is one bucket
    
 #   ACHR[b,'VisitStart'  := as.POSIXct( df[1,'tStamp']) ] 
    ACHR[b,'VisitStartInt'  :=   df[1,'tStamp'] ] 
    
    ACHR[b,'VisitStart'  :=   as.character( df[1,'tStamp'] ) ] 
    
    ACHR[b,'Visit_ID'  := df[1,'Visit_ID']]
    ACHR[b,'Subject_ID'  := df[1,'Subject_ID']]
    ACHR[b,'Clinic' := df[1,'Clinic']]
    ACHR[b,'Procedure' := df[1,'Proc']]
    ACHR[b,'Diagnosis' := df[1,'Diag']]
    ACHR[b,'Diagnosis_group'  := df[1,'Diagnosis_Group']]
    ACHR[b,'Physician'  := df[1,'Physician']]
    ACHR[b,'Weekday'  := df[1,'Weekday']]
    ACHR[b,'Month'  := df[1,'Month']]
    
    
    
} # loop thru buckets

# return the table
return(ACHR)
}


# Each bucket is a list of thread numbers that can be used to subset the list of occurrences
make_buckets_1 <- function(o, criteria){

  return( unique(o[[criteria]]) )

}





# 
# make_box_plots <- function(){
# ggboxplot(ACHR_test[NEvents>100 & Clinic=='DRH'], x = "VisitMonth", y = "NetComplexity",
#           color = "VisitDay",
#           ylab = "Complexity", xlab = "Month (DRH)")
# }
# 



# clean up the raw occurrence data
# Remove blanks for n-gram functionality
cleanOccBatch <- function(fileRows){
  
  # extract tStamp
  tStamp <- fileRows$tStamp
  
  # confirm all spaces are converted to underscores in non tStamp columns; set as factors
  cleanedCF <- data.frame(lapply(fileRows[2:ncol(fileRows)], function(x){ gsub(" ","_",x)}))
  
  # bind tStamp back to cleaned data
  complete <- cbind(tStamp,cleanedCF)
  
  # force tStamp into a "YMD_HMS" format
  # complete$tStamp <- as.character(complete$tStamp)
  # complete$tStamp <- parse_date_time(complete$tStamp, c("dmy HMS", "dmY HMS", "ymd HMS","dmy HM", "dmY HM", "ymd HM"))
  
  # add weekday and month
  # complete$weekday <- as.factor(weekdays(as.Date(complete$tStamp)))
  # complete$month   <- as.factor(months(as.Date(complete$tStamp)))
  
  return(complete)
}

get_timeScale <- function(){'hr'}

