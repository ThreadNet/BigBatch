##########################################################################################################
# THREADNET:  Batch processing for larger data sets
#
#  ACHR_Batch_PreProcess
#
# (c) 2018 Michigan State University. This software may be used according to the terms provided in the
# GNU General Public License (GPL-3.0) https://opensource.org/licenses/GPL-3.0?
# Brian Pentland
# Absolutely no warranty!
##########################################################################################################


# Separate out the two basicfunctions
# 1) reading/cleaning the occurrences
# 2) Threading the data by  adding ThreadNum and  SeqNum to  the  threads

redo_ACHR_data_from_scratch  <- function(fname){

  # Pick you point of view, or multiple POV...
  THREAD_CF = c('VISIT_ID')
  EVENT_CF = c('Action','Role','Workstation')
  
  #  first read the data
  o = read_ACHR_data( fname )
  
  # Thread occurrences for each POV  
  occ = thread_occurrences( o,  THREAD_CF, EVENT_CF ,fname )
  
  #  Now aggregate by clinic-day for  each POV
  clinic_days = ACHR_batch_clinic_days(occ, THREAD_CF, EVENT_CF, fname)
  
  # Aggregate by visit for each POV
  visits = ACHR_batch_visits(occ, THREAD_CF, EVENT_CF, fname)
  
}

# this function reads the raw data from URMC
# It returns the  data  frame  and  also saves it as an Rdata file with the  same name as the csv
read_ACHR_data <- function(fname){

  library(tidyr)
  library(data.table)
  library(dplyr)
  library(ThreadNet)
  library(ngram)
  library(lubridate)
  library(expss)
  
  
  # This code is tailored for reading in new data from URMC, October 2018
  # Assumes there  is a column  called 'Timestamps" in column 13 and  "Visit_ID"
  # read the file into data frame
  d <- fread( paste0(fname, '.csv') )
  
  # Sort by visit and timestamp
  d  <- arrange(d,desc(Visit_ID,asc(Timestamps)))
  
  # move and rename the columns, as needed
  # Mainly need  to get the  Timestamps into the first column and rename them tStamp
  cn=colnames(d)
  cnx = cn[c(13,1:12,14:15)]
  d=setcolorder(d,cnx)
  setnames(d,'Timestamps','tStamp')
  setnames(d,'X','seqn')
  
  # This converts numbers to char and replaces spaces with underscore
   d = cleanOccBatch(d)
  
  save(d, file=paste0(fname, '.rData'))
  
  return(d)
  
}

# clean up the raw occurrence data
# Remove blanks for n-gram functionality
cleanOccBatch <- function(fileRows){
  
  # extract tStamp
  tStamp <- fileRows$tStamp
  
  # First  get  rid of any underscores pre-existing in the data, then
  # confirm all spaces are converted to underscores in non tStamp columns; set as factors
  cleanedCF <- data.frame(lapply(fileRows[2:ncol(fileRows)], function(x){ gsub("_","",x)
                                                                          gsub(" ","_",x)}))
  
  # bind tStamp back to cleaned data
  complete <- cbind(tStamp,cleanedCF)
  
  
  # Old code  for forcing tStamp into a "YMD_HMS" format -- not  needed  for  URMC  EMR  data
  # complete$tStamp <- as.character(complete$tStamp)
  # complete$tStamp <- parse_date_time(complete$tStamp, c("dmy HMS", "dmY HMS", "ymd HMS","dmy HM", "dmY HM", "ymd HM"))
  
  # add weekday and month
  # complete$weekday <- as.factor(weekdays(as.Date(complete$tStamp)))
  # complete$month   <- as.factor(months(as.Date(complete$tStamp)))
  
  return(complete)
}
  
#########################################################################################
# This function adds columns  for  the the thread, as requested
# Then it sorts by thread and tStamp and  adds thread/sequence numbers to  each thread.
# 
# Typical column names from the UMRC file
# ALL_CF = c('Action','Workstation','Role','Clinic','Diagnosis_Group')
# thread_CF = c('Visit_ID', 'Role', 'Workstation')
# event_CF = c('Action','Role','Workstation')
# TN = 'Visit_ID'
# It returns  the threaded set of occurrences  and also saves it as  Rdata. 
thread_occurrences <- function(occ,THREAD_CF,EVENT_CF, fname='emr'){
  
  print(paste('Number of occurrences: ', nrow(occ)))
  
  # these will be  new column names
  new_thread_col = newColName(THREAD_CF)
  new_event_col = newColName(EVENT_CF)
  
  print(paste('new_thread_col: ',new_thread_col))
  print(paste('new_event_col: ', new_event_col))
  
  # Add names for the new columns, if necessary
  if  (!(new_thread_col  %in% colnames(occ))) {  occ = combineContextFactors(occ,THREAD_CF,new_thread_col) }
  if  (!(new_event_col  %in% colnames(occ)))  {  occ = combineContextFactors(occ,EVENT_CF,new_event_col) }
  
  # ThreadNet code assumes these columns will be there, so we need to add them
  new_occ_VR$threadNum = integer(nrow(new_occ_VR))
  new_occ_VR$seqNum =   integer(nrow(new_occ_VR))
  
  # Special case for  Visit_ID, which is already filled in URMC  EMR data
  if  (THREAD_CF == 'VISIT_ID') {
    occ$threadNum = occ$Visit_ID
    occ$seqNum = occ$seqn
    
    #  say how many there  are... 
    print(paste('Number of threads (visits) in this POV: ', length(unique( occ$Visit_ID ))))
    
    #  Now  sort the  data  set by the new  threadNum and tStamp
    occ = occ[order(occ[[new_thread_col]],occ$tStamp),]
    
  } 
  else {
  
  #  Now  sort the  data  set by the new  threadNum and tStamp
  occ = occ[order(occ[[new_thread_col]],occ$tStamp),]
  
  # get the list of unique identifies for the threads. The length of this list is the number of threads
  pov_list = unique(occ[[new_thread_col]])
  print(paste('Number of threads in this POV: ', length(pov_list)))
  
  start_row=1
  thrd=1
  for (p in pov_list){
    
    # get the length of the thread
    tlen = sum(occ[[new_thread_col]]==p)
    
    # guard against error
    if (length(tlen)==0) tlen=0
    if (tlen>0){
     
      #compute the index of the end row
      end_row = start_row+tlen-1
      
      # they all get the same thread number and incrementing seqNum
      occ[start_row:end_row, "threadNum"] <- as.matrix(rep(as.integer(thrd),tlen))
      occ[start_row:end_row, "seqNum"] <- seq(tlen)
      
      # increment the counters for the next thread
      start_row = end_row + 1
      thrd=thrd+1
      
      if (thrd %% 1000 ==0) print(paste0('Thread count =  ',thrd))
      
    }
  }
 }
  
  ### This is URMC specific:   For all POV combinations, add the clinic and clinic_day
  
  #  get the date only from the timestamp
  occ$ymd <- format(as.POSIXct(new_occ$tStamp),"%Y-%m-%d")
  
  # make new columns for clinic + day and events
   occ = unite(occ, 'Clinic_ymd', c('Clinic','ymd'),sep='_',remove = 'false')
  
    
   # Now save the result for  later and return it, as well. 
  save(occ, file=paste0(paste(fname,new_thread_col,new_event_col,sep='+'), '.Rdata'))
  print(paste('Saved threaded  occurrences: ', nrow(occ)))
  
  return(occ)
  

  ############   extra  stuff not used  #############
  # # split occ data frame by threadNum to find earliest time value for that thread
  # # then substract that from initiated relativeTime from above
  # occ_split = lapply(split(occ, occ$threadNum),
  #                    function(x) {x$relativeTime = difftime(x$relativeTime,  min(lubridate::ymd_hms(x$tStamp)), units=timescale ); x})
  # 
  # # # row bind data frame back together
  # occ= data.frame(do.call(rbind, occ_split))


}
  

  

###################################################################
###################################################################
ACHR_batch_clinic_days <- function(occ,TN, CFs) {
  
  # Name for column that has events -- three variations
  DV1= newColName(CFs[1])
  DV2= newColName(CFs[1:2])
  DV3= newColName(CFs[1:3])
  
  print(DV1)
  print(DV2)
  print(DV3)
  
  # name for column that defines the buckets
  BU = 'Clinic_ymd'
  
    # first  get the date only 
#    occ$ymd <- format(as.POSIXct(occ$tStamp),"%Y-%m-%d")
    
    # make new columns as needed for clinic + day and events
#    occ = unite(occ, 'Clinic_ymd', c('Clinic','ymd'),sep='_',remove = 'false')

    

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
                    A_NetComplexity=double(N),
                    AR_NetComplexity=double(N),
                    ARW_NetComplexity=double(N),
                    A_Nodes=double(N),
                    AR_Nodes=double(N),
                    ARW_Nodes=double(N),
                    A_Edges=double(N),
                    AR_Edges=double(N),
                    ARW_Edges=double(N),
                    A_CompressRatio = double(N),
                    AR_CompressRatio = double(N),
                    ARW_CompressRatio = double(N),
                    A_Entropy = double(N),
                    AR_Entropy = double(N),
                    ARW_Entropy = double(N),
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
      ACHR[b,A_CompressRatio := compression_index(df,DV1)]
      ACHR[b,AR_CompressRatio := compression_index(df,DV2)]
      ACHR[b,ARW_CompressRatio := compression_index(df,DV3)]  
      
      # NetComplexity of DV
      # First get the network
      n = threads_to_network_original(df,TN, DV1)
      ACHR[b,A_NetComplexity := estimate_network_complexity( n )]
      ACHR[b,A_Nodes := nrow(n$nodeDF) ]
      ACHR[b,A_Edges := nrow(n$edgeDF) ]
      
      n = threads_to_network_original(df,TN, DV2)
      ACHR[b,AR_NetComplexity := estimate_network_complexity( n )]
      ACHR[b,AR_Nodes := nrow(n$nodeDF) ]
      ACHR[b,AR_Edges := nrow(n$edgeDF) ]     
      
      n = threads_to_network_original(df,TN, DV3)
      ACHR[b,ARW_NetComplexity := estimate_network_complexity( n )]
      ACHR[b,ARW_Nodes := nrow(n$nodeDF) ]
      ACHR[b,ARW_Edges := nrow(n$edgeDF) ]  
      
      # get the entropy for AR and ARW
      ACHR[b, AR_Entropy  := compute_entropy(table(df[[DV2]])[table(df[[DV2]])>0]) ]
      ACHR[b, ARW_Entropy  := compute_entropy(table(df[[DV3]])[table(df[[DV3]])>0]) ]
      
      
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
  
  # copy this one for consistency
  ACHR[,'A_Entropy' := ACHR[,'Action_entropy'] ]
  
  # return the table
  return(ACHR)
}

###################################################################
###################################################################



# ACHR stands for Antecedents of Complexity in Healthcare Routines.  
# This is function is set up to compute process parameters on thousands of patient visits.
ACHR_batch_visits <- function(occ,TN, CFs) {

  DV1= newColName(CFs[1])
  DV2= newColName(CFs[1:2])
  DV3= newColName(CFs[1:3])
  
  print(DV1)
  print(DV2)
  print(DV3)

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
                  A_NetComplexity=double(N),
                  AR_NetComplexity=double(N),
                  ARW_NetComplexity=double(N),
                  A_Nodes=double(N),
                  AR_Nodes=double(N),
                  ARW_Nodes=double(N),
                  A_Edges=double(N),
                  AR_Edges=double(N),
                  ARW_Edges=double(N),
                  A_CompressRatio = double(N),
                  AR_CompressRatio = double(N),
                  ARW_CompressRatio = double(N),
                  A_Entropy = double(N),
                  AR_Entropy = double(N),
                  ARW_Entropy = double(N),
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
      ACHR[b,A_CompressRatio := compression_index(df,DV1)]
      ACHR[b,AR_CompressRatio := compression_index(df,DV2)]
      ACHR[b,ARW_CompressRatio := compression_index(df,DV3)]  
      
    # NetComplexity of DV
    # First get the network
      # NetComplexity of DV
      # First get the network
      n = threads_to_network_original(df,TN, DV1)
      ACHR[b,A_NetComplexity := estimate_network_complexity( n )]
      ACHR[b,A_Nodes := nrow(n$nodeDF) ]
      ACHR[b,A_Edges := nrow(n$edgeDF) ]
      
      n = threads_to_network_original(df,TN, DV2)
      ACHR[b,AR_NetComplexity := estimate_network_complexity( n )]
      ACHR[b,AR_Nodes := nrow(n$nodeDF) ]
      ACHR[b,AR_Edges := nrow(n$edgeDF) ]     
      
      n = threads_to_network_original(df,TN, DV3)
      ACHR[b,ARW_NetComplexity := estimate_network_complexity( n )]
      ACHR[b,ARW_Nodes := nrow(n$nodeDF) ]
      ACHR[b,ARW_Edges := nrow(n$edgeDF) ]  
      
      # get the entropy for AR and ARW
      ACHR[b, AR_Entropy  := compute_entropy(table(df[[DV2]])[table(df[[DV2]])>0]) ]
      ACHR[b, ARW_Entropy  := compute_entropy(table(df[[DV3]])[table(df[[DV3]])>0]) ]
      
    
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

# copy this one for consistency
ACHR[,'A_Entropy' := ACHR[,'Action_entropy'] ]

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





get_timeScale <- function(){'hr'}

