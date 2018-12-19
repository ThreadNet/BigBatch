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
library(tidyr)
library(data.table)
library(dplyr)
library(ThreadNet)
library(ngram)
library(lubridate)

redo_ACHR_data_from_scratch  <- function(fname){

  setwd("~/Documents/@NSF Healthcare Routines/Deidentified_Derm_data_Nov_2018")
  
  fname = 'audit_111818'
  
  # Pick you point of view, or multiple POV...
  THREAD_CF = c('Visit_ID')
  EVENT_CF = c('Action','Role','Workstation')
  ALL_CF = c('Action','Role','Workstation')
  
  #  first read the data
  o = read_ACHR_data( fname )  # nrow(o) = 7728927
  
  # fix the role IDs -- needs file called RoleChangeTable.csv in the local directory
  o = fix_derm_role_ID(o)
  
  # Thread occurrences  adds threadNum and seqNum to each thread, as defined by the the thread_CF.
  # Threads are always WITHIN VISITS in this code -- so they can be whole visits or chunks of visits (e.g., visit_ID_Role)
  # make two versions -- one for whole visits and one for visit_ID_Role
  ot = thread_occurrences( o, 'Visit_ID' ,fname )   # nrow(ot)= 57835
  otr = thread_occurrences( o, c('Visit_ID','Role') ,fname ) # nrow(otr) = 527666
  
  
  #  Now aggregate by clinic-day for  each POV
  clinic_days = ACHR_batch_clinic_days(occ, THREAD_CF, EVENT_CF, fname)
  
  # Aggregate by visit for each POV
  # Whole visits
  visits = ACHR_batch_visits(ot, THREAD_CF, EVENT_CF, fname)
  
  # And then by role
  VRThreads <- ACHR_batch_visit_role_threads(otr, c('Action'), c('Action','Role','Workstation') , visits)

    
  # compute  trajectories to see difference from a reference graph
  traj = graph_trajectory( ot, 'Clinic_ymd', EVENT_CF, 2, 1, 0, 'Clinic_trajectory')
  
}

# this function reads the raw data from URMC
# It returns the  data  frame  and  also saves it as an Rdata file with the  same name as the csv
# Tested Nov 16.
read_ACHR_data <- function(fname){


  
  
  # This code is tailored for reading in new data from URMC, October 2018
  # Assumes there  is a column  called 'Timestamps" in column 13 and  "Visit_ID"
  # read the file into data.table
   d <-   fread( paste0(fname, '.csv') )
  
  # Sort by visit and timestamps (will do this later with setkey)
   d  <- arrange(d,desc(Visit_ID,asc(Timestamps)))

  # This  stuff  is URMC specific 
  #
  # move and rename the columns, as needed
  # Mainly need  to get the  Timestamps into the first column and rename them tStamp
  cn=colnames(d)
  cnx = cn[c(13,1:12,14:15)]
  d=setcolorder(d,cnx)
  setnames(d,'Timestamps','tStamp')
  setnames(d,'V1','seqn')

   # Fix the name of the Highland  clinic
   d$Clinic = gsub('HH POB','HHPOB',d$Clinic)
  
  
  # This converts numbers to char and replaces spaces with underscore
    d = cleanOccBatch(d)
 
    ### This is URMC specific: Add the clinic and clinic_day  and fix the code for the highland  clinic
   
   # get the date only from the timestamp
    d$ymd <- format(  as.POSIXct(d$tStamp, format="%Y-%m-%d"),  "%Y-%m-%d" ) 
   
   # make new columns for clinic + day and events
   d = unite(d, 'Clinic_ymd', c('Clinic','ymd'),sep='_',remove = 'false')
   
   
  # Save the result   
  save_file_name = paste0(fname, '.rData')
  save(d, file=save_file_name) 
  print(paste("Saved ",nrow(d), " occurrences in ",save_file_name))
  
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
  complete <- as.data.table(cbind(tStamp,cleanedCF))
  
  # Old code  for forcing tStamp into a "YMD_HMS" format -- not  needed  for  URMC  EMR  data
   # complete$tStamp <- as.character(complete$tStamp)
   # complete$tStamp <- parse_date_time(complete$tStamp, c("dmy HMS", "dmY HMS", "ymd HMS","dmy HM", "dmY HM", "ymd HM"))
   
  # add weekday and month
  # complete$weekday <- as.factor(weekdays(as.Date(complete$tStamp)))
  # complete$month   <- as.factor(months(as.Date(complete$tStamp)))
  
  return(complete)
}
  
#######################################################################################
# In the derm data, the role "technician" was used for people who were not technicians
# Pass in the table of occurrences, fix it, and return it
# read in the list of changes
fix_derm_role_ID <- function(o){
  
  # read in the new role ID
  rc = read.csv('RoleChangeTable.csv', stringsAsFactors = FALSE)
  
  
  #  occurrences are a data.table so use :=
  for (r in 1:nrow(rc)) {
    
    old_rid = rc$OLD_Role_ID[r]
    new_rid = rc$NEW_Role_ID[r]
    new_role = rc$NEW_Role[r]
    
    # pick out rows with old Role_ID
    # old_rid_vector = o$Role_ID==old_rid
    # 
    # o[old_rid_vector, Role_ID := new_rid]
    # o[old_rid_vector, Role := new_role]
    
    # crazy data table  syntax... do it all in one shot
    o[o$Role_ID==old_rid, ':=' (Role = new_role, Role_ID = new_rid ) ]
    
  }
  
  return(o)
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
  
  # Assume the data are sorted by visit_ID and tStamp -- the natural flow of the visit.
  idx_list = which(occ[[new_thread_col]] !=dplyr::lag(occ[[new_thread_col]]))
  idx_list = c(idx_list, nrow(occ))
  
  # creates two lists of indices so they match up correctly. 
  start_idx = c(1, head(idx_list,-1))  # add 1 to the front and drop the last
  end_idx =   idx_list - 1    # shift back by one
  
  print(paste('Number of threads in this POV: ', length(idx_list)))
  
  # sapply(seq(length(idx_list)), function(x)
    for (x in seq(length(idx_list)))
    {   if (x %% 1000 ==0) print(paste0('Thread count =  ',x))
     occ[start_idx[x]:end_idx[x],'threadNum' := x ]
     occ[start_idx[x]:end_idx[x],'seqNum':= c(1:(end_idx[x]-start_idx[x]+1)) ]
    }
  
  # Now save the result for  later and return it, as well. 
  #  save(occ, file=paste0(paste(fname,new_thread_col,new_event_col,sep='+'), '.Rdata'))
  
  save_file_name = paste0(paste0(paste(fname,new_thread_col,sep='+'), '.rData'))
  save(occ, file=save_file_name)
  
  print(paste('Saved threaded occurrences: ', nrow(occ), " in ",save_file_name))

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



  

  
