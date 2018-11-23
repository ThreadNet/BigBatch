##########################################################################################################
# THREADNET:  Batch processing for larger data sets
#
#  ACHR_Batch_PreProcess
#
# (c) 2018 Michigan State University. This software may be used according to the terms provided in the
# GNU General Public License (GPL-3.0) https://opensource.org/licenses/GPL-3.0?
# Brian Pentland
# Absolutely no warranty!
#
# Separate out the two basic functions
# 1) reading/cleaning the occurrences
# 2) Threading the data by  adding ThreadNum and  SeqNum to  the  threads
#
##########################################################################################################


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
# Tested Nov 16.
read_ACHR_data <- function(fname){

  library(tidyr)
  library(data.table)
  library(dplyr)
  library(ThreadNet)
  library(ngram)
  library(lubridate)
  library(anytime)
  
  
  # This code is tailored for reading in new data from URMC, October 2018
  # Assumes there  is a column  called 'Timestamps" in column 13 and  "Visit_ID"
  # read the file into data frame
   d <-   fread( paste0(fname, '.csv') )
  
  # Sort by visit and timestamp
   d  <- arrange(d,desc(Visit_ID,asc(Timestamps)))
  
  ##################################################################################
  # This  stuff  is URMC specific 
  #
  # move and rename the columns, as needed
  # Mainly need  to get the  Timestamps into the first column and rename them tStamp
  cn=colnames(d)
  cnx = cn[c(13,1:12,14:15)]
  d=setcolorder(d,cnx)
  setnames(d,'Timestamps','tStamp')
  setnames(d,'V1','seqn')

  # This converts numbers to char and replaces spaces with underscore
   d = cleanOccBatch(d)
  
    ### This is URMC specific: Add the clinic and clinic_day  and fix the code for the highland  clinic
   
   # get the date only from the timestamp
    d$ymd <- format(  as.POSIXct(d$tStamp, format="%Y-%m-%d"),  "%Y-%m-%d" ) 
   
   # make new columns for clinic + day and events
   d = unite(d, 'Clinic_ymd', c('Clinic','ymd'),sep='_',remove = 'false')
   
   # Fix the name of the Highland  clinic
   d$Clinic = gsub('HH_POB','HHPOB',d$Clinic)
   
  # Save the result   
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
  cleanedCF <- data.frame(lapply(fileRows[2:ncol(fileRows)], function(x){  gsub(" ","_",x)  }) )
  
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
# This function adds columns  for  the the thread, as requested ( Tested Nov 16 )
# Then it sorts by thread and tStamp and  adds thread/sequence numbers to  each thread.
# 
# Typical column names from the UMRC file
# ALL_CF = c('Action','Workstation','Role','Clinic','Diagnosis_Group')
# thread_CF = c('Visit_ID', 'Role', 'Workstation')
# event_CF = c('Action','Role','Workstation')

# It returns  the threaded set of occurrences  and also saves it as  Rdata. 
thread_occurrences <- function(occ, THREAD_CF, fname='emr'){
  
  print(paste('Number of occurrences: ', nrow(occ)))
  
  # these will be  new column names
  new_thread_col = newColName(THREAD_CF)
 # new_event_col = newColName(EVENT_CF)
  
  print(paste('new_thread_col: ',new_thread_col))
 # print(paste('new_event_col: ', new_event_col))
  
  # Add names for the new columns, if necessary
  if  (!(new_thread_col  %in% colnames(occ))) {  occ = combineContextFactors(occ,THREAD_CF,new_thread_col) }
# if  (!(new_event_col  %in% colnames(occ)))  {  occ = combineContextFactors(occ,EVENT_CF,new_event_col) }
  
  # ThreadNet code assumes these columns will be there, so we need to add them
  occ$threadNum = integer(nrow(occ))
  occ$seqNum =   integer(nrow(occ))
  
  #_#_#_#_#_#_ I think what we need to do to get the correct visit_role is sort by  visit and time stamp
  #_#_#_#_#_#_ But then SPLIT the visit according to visit_role or visit_workstation or whatever
  #  Now  sort the  data  set by the new  threadNum and tStamp
  occ = occ[order(occ[[new_thread_col]],occ$tStamp),]
  
  
  #_#_#_#_#_#_ This will be incorrect if  the same role is occurs more than once in a single visit. 
  #_#_#_#_#_#_ Need to break  it out according to the changes in role within a visit...  
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
  
   # Now save the result for  later and return it, as well. 
#  save(occ, file=paste0(paste(fname,new_thread_col,new_event_col,sep='+'), '.Rdata'))
  save(occ, file=paste0(paste(fname,new_thread_col,sep='+'), '.Rdata'))
  
  print(paste('Saved threaded  occurrences: ', nrow(occ)))
  
  return(occ)
}


  ############   extra  stuff not used  #############
  # split occ data frame by threadNum to find earliest time value for that thread
  # # then substract that from initiated relativeTime from above
  # occ_split = lapply(split(occ, occ$threadNum),
  #                    function(x) {x$relativeTime = difftime(x$relativeTime,  min(lubridate::ymd_hms(x$tStamp)), units=timescale ); x})
  # 
  # # # row bind data frame back together
  # occ= data.frame(do.call(rbind, occ_split))



  

  
