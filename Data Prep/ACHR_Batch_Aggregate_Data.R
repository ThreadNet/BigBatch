##########################################################################################################
# THREADNET:  Batch processing for larger data sets
#
#  ACHR_Batch_Aggregate_Data
#
# (c) 2018 Michigan State University. This software may be used according to the terms provided in the
# GNU General Public License (GPL-3.0) https://opensource.org/licenses/GPL-3.0?
# Brian Pentland
# Absolutely no warranty!
##########################################################################################################


# Two  functions to aggregate the occurrences 
# 1) By thread (usually visit, but could be visit_role)
# 2) By Clinic_day
  

###################################################################
###################################################################
# This is function is set up to aggregate the occurrences within threads -- typically visits
# occ = pre-processed threaded occurrences
# TN = threadNum in most cases
# CFs can be chosen
# event_CFs define changes within threads
# ALL_CFs are used to computer the CF_alignment


#  need to bring Role_ID forward when looking at Visit_Role, so we can track residents over time.

ACHR_batch_threads <- function(occ,THREAD_CFs, EVENT_CFs, ALL_CFs) {
  
  library(tidyr)
  library(dplyr)
  library(ThreadNet)
  library(ngram)
  library(lubridate)
  library(stringr)
  library(data.table)
  
  # Add  columns for combinations of CFs if needed
  thread_col = newColName(THREAD_CFs)
  TN = thread_col
  if  (!(thread_col  %in% colnames(occ)))  {  occ = combineContextFactors(occ,THREAD_CFs,thread_col) }

  
  new_event_col = newColName(EVENT_CFs)
  if  (!(new_event_col  %in% colnames(occ)))  {  occ = combineContextFactors(occ,EVENT_CFs,new_event_col) }
 
  all_cf_col = newColName(ALL_CFs)
  if  (!(all_cf_col  %in% colnames(occ)))  {  occ = combineContextFactors(occ,ALL_CFs,all_cf_col) } 
  
  # set key on the data.table
  setkeyv(occ, TN)
  
  # get the list of buckets
  bucket_list <- unique(occ[[TN]])
  
  # print the number of buckets
  print(paste0('Number of buckets=', length(bucket_list)  ))
  
  
  # make data frame with results
 Thrds = data.table( t(sapply( bucket_list, 
                           function(b){
                             
                             # select a subset 
                             df= occ[get(TN)==b]
                    
                            # make sure it  is sorted by timestamp
                             df=df[order(df$tStamp),]
                             
                             # get the network -- only if there are enough rows...
                             
                              if (nrow(df)>2)   n = threads_to_network_original(df,TN, new_event_col) 
                              else  n = list(edgeDF=t(c(0)),nodeDF=t(c(0)))
                             
                             # compute each parameter and put them in a vector
                             c(
                               bucket=b,
                               Clinic = as.character(df$Clinic[1]),
                               Physician  = as.character(df$Physician[1]),
                               ymd = as.character(df$ymd[1]),
                               Clinic_ymd = as.character(df$Clinic_ymd[1]),
                               Weekday  = as.character(df$Weekday[1]),
                               Month  = as.character(df$Month[1]),
                               Phase =  compute_phase(df$tStamp[1]),
                               threadNum = as.numeric(df$threadNum[1]),
                               Visit_ID  = as.character(df$Visit_ID[1]),
                               Subject_ID  = as.character( df$Subject_ID[1] ),
                               Role_ID = as.character( df$Role_ID[1] ),
                               NEvents = nrow(df),
                               ThreadStart= as.character(df$tStamp[1] ),
                               ThreadDuration= compute_thread_duration(df), 
                               VisitDuration= compute_visit_duration(df), 
                               wait_time1 = compute_wait_time1(df),
                               wait_time2 = compute_wait_time2(df),
                               Visit_number = as.numeric(df[1,'Visitnum_Total']),
                               LOC_CPT =  as.character(df[1,'LOS_CPT']),
                               NetComplexity=estimate_network_complexity( n ),
                               Nodes=nrow(n$nodeDF),
                               Edges=nrow(n$edgeDF),
                               CompressRatio = compression_index(df,new_event_col),
                               Entropy = compute_entropy(table(df[[new_event_col]])[table(df[[new_event_col]])>0]),
                               NumProcedures = count_procedures(df$Proc[1]),
                               NumDiagnoses = count_diagnoses(df$Diag[1]),
                               Proc = as.character(df[1,'Proc']),
                               Diagnosis = as.character(df[1,'Diag']),    
                               Diagnosis_group  = as.character(df[1,'Diagnosis_Group']),
                               CF_Alignment = 1,    # make placeholder, but compute below
                               ALL_CF_count = length(unique(df[[all_cf_col]])),
                               ALL_CF_entropy = compute_entropy(table(df[[all_cf_col]])[table(df[[all_cf_col]])>0]), 
                             
                             # name these columns afterwards
                             sapply(ALL_CFs, function(cf){
                               c( length(unique(df[[cf]]))
                                  # ,
                                  # compute_entropy(table(df[[cf]])[table(df[[cf]])>0]) 
                                  ) })
                             
                             )
                      } )))
 
 
 # name the last columns -- code has to match above
 cn =  as.vector(sapply(ALL_CFs, function(cf){
   c( paste0(cf,"_count")
    #   ,  
    # paste0(cf,"_entropy") 
    ) }))
 
 # now assign them to the last columns
 last_col = ncol(Thrds)
 first_col= last_col-length(cn)+1
 setnames(Thrds, c(first_col:last_col), cn)
 
 # Compute the alignment of the context factors
 Thrds$CF_Alignment =  as.numeric( as.character(Thrds$Action_count)) / as.numeric( as.character(Thrds$ALL_CF_count ))
 
 save_file_name = paste0(paste('Thrds',TN,new_event_col,sep='+'), '.Rdata')
 save(Thrds, file=save_file_name)
 
 print(paste('Saved ', nrow(Thrds), " records in ",save_file_name))
 
 return(Thrds) 
  
}



get_timeScale <- function(){'hr'}



# Need to write all  these functions
compute_phase <- function(t){
  
  t=as.Date(t)
  
  if ( t >= '2016-01-01' & t < '2016-06-07')  return('one')
  if ( t >= '2016-06-07' & t < '2016-09-01')  return('two')
  if ( t >= '2016-09-01' & t < '2017-04-15')  return('three')
  if ( t >= '2017-04-15' & t < '2017-09-01')  return('four')
  if ( t >= '2017-09-01' & t < '2017-12-31')  return('five')
    
    # else return NA
    return(NA)

}
compute_thread_duration <- function(df){
  
  return(difftime(lubridate::ymd_hms(df$tStamp[nrow(df)]),  lubridate::ymd_hms(df$tStamp[1]), units='hours' ))
  
}

compute_visit_duration <- function(df){
  
  # get the checkin time.  Note that they sometimes look at patient data in advance of the visit, so you 
  # cannot use the first occurrences to mark the start of the thread
  w1= grep('CHECKIN_TIME',df$Action)[1]
  w2= grep('AVS_PRINT_TIME',df$Action)[1]
  
  return(difftime(lubridate::ymd_hms(df$tStamp[w2]),  lubridate::ymd_hms(df$tStamp[w1]), units='hours' ))
  
}

# First wait time... from checkin until they get their vitals taken
compute_wait_time1 <- function(df){
  
  #  find the first  occurrence of MR_VN_VITALS. Corresponds closely to wait time
 # w= grep('MR_VN_VITALS',df$Action)[1]
  w= grep('VITALS',df$Action)[1]
  
  return( difftime( lubridate::ymd_hms(df$tStamp[w]),lubridate::ymd_hms(df$tStamp[1]),  units='hours' ) )
   
}

# Second wait time... from  vitals until they see a doctor or other medical person
compute_wait_time2 <- function(df){
  
  #  find the first  occurrence of MR_VN_VITALS. Get the time t1
 # w1= grep('MR_VN_VITALS',df$Action)[1]
  w1= grep('VITALS',df$Action)[1]
  
  if (is.na(w1))
    return(NA)
  else
  {
  t1 = lubridate::ymd_hms(df$tStamp[w1])
  
  
  # get the workstation and role used for the vitals at t1
  work1 = df$Workstation[w1]
  role1 = df$Role[w1]
  
  # get set of actions at the same workstation or different workstations that happen later in the visit
  # Truncate the visit up to w1
  restdf = df[w1:nrow(df),]

  # then get the set of action at the same workstation or different workstations (in that visit)
  sameWS =  restdf[restdf$Workstation == work1,]
  differentWS =  restdf[restdf$Workstation != work1,]
  
  # get list of other roles
  otherRoles = paste(setdiff(unique(restdf$Role),c('Technician','Technologist','Staff','Unknown')),
                    collapse = '|')
  
  # Now find the marker for when they see a physician
  # procedure is differenrt in BRKPT and HHPOB clinics
  if (df$Clinic[1]  %in% c('BRKPT','HHPOB') )
    { # LPN brings a different workstation on wheels at these clinics  
    firstLPN =  grep('LPN',differentWS$Role)[1] 
    t2 = lubridate::ymd_hms(differentWS$tStamp[firstLPN])
    
    # Other roles use the same workstation
    firstDoc =  grep('Physician|Resident|Registered', sameWS$Role)[1] 
    t2 = lubridate::ymd_hms(sameWS$tStamp[firstDoc])
    }
  else 
  {  # a different role, using the same workstation 
    w2 = max(w1, grep(otherRoles, sameWS$Role) )[1]
    t2 = lubridate::ymd_hms(sameWS$tStamp[w2])
    }
    
  # return the difference in time from t1 to t2
  return( difftime( t2, t1,  units='hours' ) )
  }
}


# THIS ONE IS GOOD.
count_procedures <- function( p ){
  
  # if no procedures, return zero
  if (is.na( p )) return( 0 )
  
  # get the overall number of items
  total_num = str_count( p, '#@#')
  
  # adjust for the visit codes that are not actual procedures
  v=c('99211','99212','99213','99214','99215','99201','99202','99203','99204','99205')
  
  # just look at the first occurrence in the df
  num_visit_codes  = sum(str_count( p,  v ))
  
  # limit  to non-negative 
  total_num = max(0, total_num - num_visit_codes) 
  
  return( total_num  )
}

# This one  is easy
count_diagnoses <- function(d){
  
  # Sometimes the marker is missing, so set the floor to one.  There is always at least one. 
  return( max(1, str_count( d, '#@#' )))
}






###################################################################
###################################################################
# This is function is set up to aggregate the occurrences among collections of visit -- typically clinic_days
# occ = pre-processed threaded occurrences
# TN = threadNum in most cases
# CFs can be chosen  
ACHR_batch_clinic_days <- function(occ,TN='Clinic_ymd', EVENT_CFs, ALL_CFs) {
  
  library(tidyr)
  library(dplyr)
  library(ThreadNet)
  library(ngram)
  library(lubridate)
  library(stringr)
  library(data.table)
  
  # make a list of unique buckets.   
  bucket_list <- unique(occ[[TN]])
  
  # print the number of buckets
  print(paste0('Number of buckets=', length(bucket_list)  ))
  
  # Add  columns for combinations of CFs if needed
  new_event_col = newColName(EVENT_CFs)
  if  (!(new_event_col  %in% colnames(occ)))  {  occ = combineContextFactors(occ,EVENT_CFs,new_event_col) }
  
  all_cf_col = newColName(ALL_CFs)
  if  (!(all_cf_col  %in% colnames(occ)))  {  occ = combineContextFactors(occ,ALL_CFs,all_cf_col) } 
  
  
  # make data frame with results
 cds = data.table( t(sapply( bucket_list, 
                                function(b){
                                  
                                  # select a subset of occurrences for the bucket
                                  df = occ[ occ[[TN]] ==b , ] 
                                  
                                  
                                  # get the network -- only if there are enough rows...
                                  if (nrow(df)>2)   n = threads_to_network_original(df,TN, new_event_col) 
                                  else  n = list(edgeDF=t(c(0)),nodeDF=t(c(0)))
                                  
                                  # compute each parameter and put them in a vector
                                  c(
                                    bucket=b,
                                    Clinic = df[1,'Clinic'],
                                    Physician  = as.character(df[1,'Physician']),
                                    ymd = df[1,'ymd'],
                                    Clinic_ymd = df[1,'Clinic_ymd'],
                                    Weekday  = df[1,'Weekday'],
                                    Month  = df[1,'Month'],
                                    Phase =  compute_phase(df$tStamp[1]),
                                    threadNum = df[1,'threadNum'],
                                    Visit_ID  = df[1,'Visit_ID'],
                                    Subject_ID  = df[1,'Subject_ID'],
                                    NEvents = nrow(df),
                                    bucketStart= as.character( df[1,'tStamp'] ),
                                    bucketStartInt = df[1,'tStamp'],
                                    bucketDuration= difftime(max(lubridate::ymd_hms(df$tStamp)),  min(lubridate::ymd_hms(df$tStamp)), units='hours' ),
                                    #   wait_time = compute_wait_time(df),
                                    #   Visit_number = integer(N),
                                    #   LOC_CPT =  character(N),
                                    #  NumUniqueProcedures = length(unique(df[['Proc']])
                                    NumVisits =  length(unique(df[['Visit_ID']])),
                                    NumUniqueDiagnosisGroups = length(unique(df[['Diagnosis_Group']])), 
                                    NumPhysicians = length(unique(df[['Physician']])),
                                    TotalStaff = length(unique(df[['Role_ID']])),
                                    NetComplexity=estimate_network_complexity( n ),
                                    Nodes=nrow(n$nodeDF),
                                    Edges=nrow(n$edgeDF),
                                    CompressRatio = compression_index(df,new_event_col),
                                    Entropy = compute_entropy(table(df[[new_event_col]])[table(df[[new_event_col]])>0]),
                                    # NumProcedures = count_procedures(df$Proc[1]),
                                    # NumDiagnoses = count_diagnoses(df$Diag[1]),
                                    # Diagnosis_group  = df[1,'Diagnosis_Group'],
                                    CF_Alignment = 1,    # make placeholder, but compute below
                                    ALL_CF_count = length(unique(df[[all_cf_col]])),
                                    ALL_CF_entropy = compute_entropy(table(df[[all_cf_col]])[table(df[[all_cf_col]])>0]), 
                                    
                                    # name these columns afterwards
                                    sapply(ALL_CFs, function(cf){
                                      c( length(unique(df[[cf]]))
                                         # ,
                                         # compute_entropy(table(df[[cf]])[table(df[[cf]])>0]) 
                                         ) })
                                    
                                  )
                                } )))
  
  
  # name the last columns -- code has to match above
  cn =  as.vector(sapply(ALL_CFs, function(cf){
    c( paste0(cf,"_count")
       # ,  
       # paste0(cf,"_entropy") 
       ) }))
  
  # now assign them to the last columns
  last_col = ncol(cds)
  first_col= last_col-length(cn)+1
  setnames(cds, c(first_col:last_col), cn)
  
  # Compute the alignment of the context factors
  cds$CF_Alignment =  as.numeric( as.character(cds$Action_count)) / as.numeric( as.character(cds$ALL_CF_count ))
  
  save(cds, file=paste0(paste('cds',TN,new_event_col,sep='+'), '.Rdata'))
  
  return(cds) 
  
}


###################################################################
###################################################################
##  OLD   ###
##  OLD   ###
##  OLD   ###
# ACHR_batch_threads_old <- function(occ,TN, EVENT_CFs, ALL_CFs) {
# 
#   library(tidyr)
#   library(dplyr)
#   library(ThreadNet)
#   library(ngram)
#   library(lubridate)
#   library(stringr)
#   
#   new_event_col = newColName(EVENT_CFs)
#   if  (!(new_event_col  %in% colnames(occ)))  {  occ = combineContextFactors(occ,EVENT_CFs,new_event_col) }
#   
# # pick subsets --here we aggregate one thread at a time  
# bucket_list <- make_buckets_1(occ, TN)
# 
# # get the size (number of buckets)
# N = length(bucket_list)
# c
# # pre-allocate the result
# ACHR = data.table(bucket=integer(N),
#                   NEvents = integer(N),
#                   ThreadStart= character(N),
#                   ThreadStartInt = integer(N),
#                   ThreadDuration =double(N),   
#                   NetComplexity=double(N),
#                   Nodes=double(N),
#                   Edges=double(N),
#                   CompressRatio = double(N),
#                   Entropy = double(N),
#                   NumProcedures = double(N),
#                   NumDiagnoses = double(N),
#                   wait_time = double(N),
#                   Visit_number = integer(N),
#                   LOC_CPT =  character(N), 
#                   Visit_ID  = character(N), 
#                   Subject_ID  = character(N), 
#                   Clinic = character(N),
#                   ymd = character(N),
#                   Clinic_ymd = character(N),
#                   Procedure = character(N),
#                   Diagnosis = character(N),   # diag in the raw data
#                   Diagnosis_group  = character(N), 
#                   Physician  = character(N), 
#                   Weekday  = character(N), 
#                   Month  = character(N),
#                   Phase =  character(N),
#                   CF_Alignment = double(N) 
#                   )
# 
# # Now add columns for the CFs.  There will be two for each CF
# for (cf in ALL_CFs){
#   ACHR[, paste0(cf,"_count"):= double(N)]
#   ACHR[, paste0(cf,"_entropy"):= double(N)]
# }
# 
# # loop through the buckets. Result will be data frame with one row per bucket
# for (i in 1:N){
# 
#   b = i #  as.integer(bucket_list[i])
# 
#   # print once every 100 visits
#   if (b%%100==0) {print(b)}
#   
#   # select the threads that go in this bucket
#     df = occ[occ[[TN]] ==bucket_list[i],]
# 
#     # bucket number
#     ACHR[b,bucket := b]
# 
#     # length of the thread (number of rows)
#     ACHR[b,NEvents := nrow(df)]
# 
#     # only do the computations if there are more than two occurrences
#     if (nrow(df) > 2) {
# 
#       # compute the duration of the visit in hours
#       ACHR[b,ThreadDuration := difftime(max(lubridate::ymd_hms(df$tStamp)),  min(lubridate::ymd_hms(df$tStamp)), units='hours') ]
#      
#     # compressibility of DV
#       ACHR[b,CompressRatio := compression_index(df,new_event_col)]
#       # ACHR[b,AR_CompressRatio := compression_index(df,DV2)]
#       # ACHR[b,ARW_CompressRatio := compression_index(df,DV3)]  
#       
#     # NetComplexity of DV
#     # First get the network
#       # NetComplexity of DV
#       # First get the network
#       n = threads_to_network_original(df,TN, new_event_col)
#       ACHR[b,NetComplexity := estimate_network_complexity( n )]
#       ACHR[b,Nodes := nrow(n$nodeDF) ]
#       ACHR[b,Edges := nrow(n$edgeDF) ]
#       
#     
#   # compute stuff on each context factor
#   for (cf in ALL_CFs){
# 
#     # Count the unique elements in each cf
#     ACHR[b, paste0(cf,"_count") :=  length(unique(df[[cf]])) ]
# 
#     # get the entropy
#     ACHR[b, paste0(cf,"_entropy") := compute_entropy(table(df[[cf]])[table(df[[cf]])>0]) ]
# 
#   }
# } # kf nrows > 2
# 
#     # Now copy in the rest of data  
#     # this only works because one visit is one bucket
#     
#  #   ACHR[b,'VisitStart'  := as.POSIXct( df[1,'tStamp']) ] 
#     ACHR[b,'ThreadStartInt'  :=   df[1,'tStamp'] ] 
#     
#     ACHR[b,'ThreadStart'  :=   as.character( df[1,'tStamp'] ) ] 
#     
#     ACHR[b,'Visit_ID'  := df[1,'Visit_ID']]
#     ACHR[b,'Subject_ID'  := df[1,'Subject_ID']]
#     ACHR[b,'Clinic' := df[1,'Clinic']]
#     ACHR[b,'Procedure' := df[1,'Proc']]
#     ACHR[b,'Diagnosis' := df[1,'Diag']]
#     ACHR[b,'Diagnosis_group'  := df[1,'Diagnosis_Group']]
#     ACHR[b,'Physician'  := df[1,'Physician']]
#     ACHR[b,'Weekday'  := df[1,'Weekday']]
#     ACHR[b,'Month'  := df[1,'Month']]
#     
#     
#     
# } # loop thru buckets
# 
# 
#   save(ACHR, file=paste0(paste('Threads',TN,new_event_col,sep='+'), '.Rdata'))
# 
# 
# # return the table
# return(ACHR)
# }


###################################################################
###################################################################
# ACHR_batch_clinic_days <- function(occ,TN, EVENT_CFs, ALL_CFs) {
#   
#   # Name for column that has events -- three variations
#   DV1= newColName(CFs[1])
#   DV2= newColName(CFs[1:2])
#   DV3= newColName(CFs[1:3])
#   
#   print(DV1)
#   print(DV2)
#   print(DV3)
#   
#   # name for column that defines the buckets
#   BU = 'Clinic_ymd'
#   
#   # first  get the date only 
#   #    occ$ymd <- format(as.POSIXct(occ$tStamp),"%Y-%m-%d")
#   
#   # make new columns as needed for clinic + day and events
#   #    occ = unite(occ, 'Clinic_ymd', c('Clinic','ymd'),sep='_',remove = 'false')
#   
#   
#   
#   # pick subsets -- one visit at a time in this version, but could be more
#   bucket_list <- unique(occ[['Clinic_ymd']])
#   
#   # get the size (number of buckets)
#   N = length(bucket_list)
#   print(N)
#   
#   # pre-allocate the data.table.  Tables are supposed to be faster.
#   ACHR = data.table(bucket=integer(N),
#                     Clinic_date =  character(N),
#                     YMD_date = character(N),
#                     NEvents = integer(N),
#                     NumVisits = numeric(N),
#                     A_NetComplexity=double(N),
#                     AR_NetComplexity=double(N),
#                     ARW_NetComplexity=double(N),
#                     A_Nodes=double(N),
#                     AR_Nodes=double(N),
#                     ARW_Nodes=double(N),
#                     A_Edges=double(N),
#                     AR_Edges=double(N),
#                     ARW_Edges=double(N),
#                     A_CompressRatio = double(N),
#                     AR_CompressRatio = double(N),
#                     ARW_CompressRatio = double(N),
#                     A_Entropy = double(N),
#                     AR_Entropy = double(N),
#                     ARW_Entropy = double(N),
#                     Clinic = character(N),
#                     NumUniqueProcedures = numeric(N),
#                     NumUniqueDiagnosisGroups  = numeric(N), 
#                     NumPhysicians  = numeric(N), 
#                     Weekday  = character(N), 
#                     Month  = character(N)
#   )
#   
#   # Now add columns for the IVs.  There will be three for each IV
#   
#   # Add the IV columns
#   for (cf in CFs){
#     
#     ACHR[, paste0(cf,"_count"):= double(N)]
#     ACHR[, paste0(cf,"_compression"):= double(N)]
#     ACHR[, paste0(cf,"_entropy"):= double(N)]
#     
#   }
#   
#   # loop through the buckets. Result will be data frame with one row per bucket
#   for (i in 1:N){
#     
#     b = i #  as.integer(bucket_list[i])
#     
#     # print once every 10 buckets
#     if (b%%10==0) {print(b)}
#     
#     # select the threads that go in this bucket
#     df = occ[occ[[BU]] ==bucket_list[i],]
#     
#     # bucket number
#     ACHR[b,bucket := b]
#     
#     # length of the thread (number of rows)
#     ACHR[b,NEvents := nrow(df)]
#     
#     # only do the computations if there are more than two occurrences
#     if (nrow(df) > 2) {
#       
#       # compute the duration of the visit in hours
#       #     ACHR[b,VisitDuration := difftime(max(lubridate::ymd_hms(df$tStamp)),  min(lubridate::ymd_hms(df$tStamp)), units='hours') ]
#       
#       # compressibility of DV
#       ACHR[b,A_CompressRatio := compression_index(df,DV1)]
#       ACHR[b,AR_CompressRatio := compression_index(df,DV2)]
#       ACHR[b,ARW_CompressRatio := compression_index(df,DV3)]  
#       
#       # NetComplexity of DV
#       # First get the network
#       n = threads_to_network_original(df,TN, DV1)
#       ACHR[b,A_NetComplexity := estimate_network_complexity( n )]
#       ACHR[b,A_Nodes := nrow(n$nodeDF) ]
#       ACHR[b,A_Edges := nrow(n$edgeDF) ]
#       
#       n = threads_to_network_original(df,TN, DV2)
#       ACHR[b,AR_NetComplexity := estimate_network_complexity( n )]
#       ACHR[b,AR_Nodes := nrow(n$nodeDF) ]
#       ACHR[b,AR_Edges := nrow(n$edgeDF) ]     
#       
#       n = threads_to_network_original(df,TN, DV3)
#       ACHR[b,ARW_NetComplexity := estimate_network_complexity( n )]
#       ACHR[b,ARW_Nodes := nrow(n$nodeDF) ]
#       ACHR[b,ARW_Edges := nrow(n$edgeDF) ]  
#       
#       # get the entropy for AR and ARW
#       ACHR[b, AR_Entropy  := compute_entropy(table(df[[DV2]])[table(df[[DV2]])>0]) ]
#       ACHR[b, ARW_Entropy  := compute_entropy(table(df[[DV3]])[table(df[[DV3]])>0]) ]
#       
#       
#       # compute stuff on each context factor
#       for (cf in CFs){
#         
#         # Count the unique elements in each cf
#         ACHR[b, paste0(cf,"_count") :=  length(unique(df[[cf]])) ]
#         
#         # get the compression
#         ACHR[b, paste0(cf,"_compression") := compression_index(df,cf) ]
#         
#         # get the entropy
#         ACHR[b, paste0(cf,"_entropy") := compute_entropy(table(df[[cf]])[table(df[[cf]])>0]) ]
#         
#       }
#     } # kf nrows > 2
#     
#     # Now copy in the rest of data  
#     # this works because one visit is one bucket
#     # count the number of different diagnoses --> typical indicators of complexity
#     
#     
#     ACHR[b,'Clinic_date'  := df[1,'Clinic_ymd']]
#     ACHR[b,'YMD_date'  := df[1,'ymd']]
#     ACHR[b,'NumVisits'  :=  length(unique(df[['Visit_ID']]))]
#     ACHR[b,'Clinic' := df[1,'Clinic']]
#     ACHR[b,'NumUniqueProcedures' := length(unique(df[['Proc']]))]
#     ACHR[b,'NumUniqueDiagnosisGroups'  := length(unique(df[['Diagnosis_Group']]))]  
#     ACHR[b,'NumPhysicians'  := length(unique(df[['Physician']]))]  
#     ACHR[b,'Weekday'  := df[1,'Weekday']]
#     ACHR[b,'Month'  := df[1,'Month']]
#     
#     
#     
#   } # loop thru buckets
#   
#   # copy this one for consistency
#   ACHR[,'A_Entropy' := ACHR[,'Action_entropy'] ]
#   
#   # return the table
#   return(ACHR)
# }

