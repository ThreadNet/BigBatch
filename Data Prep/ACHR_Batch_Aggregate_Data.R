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
  
# Rewrite for speed: 
# a) go back to for loop, but just for the stuff that is computed based on the df
# b) grab first row of each thread seqNum ==1 from the table
# c) merge them back  together
# d) Probably need to do different versions for visit vs. visit_role because a lot of the parameters make no sense...

# April 2019 Added  compute_graph_entropy

###################################################################
###################################################################
# This is function is set up to aggregate the occurrences within threads -- typically visits
# occ = pre-processed threaded occurrences
# TN = threadNum in most cases
# CFs can be chosen
# event_CFs define changes within threads
# ALL_CFs are used to compute the CF_alignment



ACHR_batch_threads <- function(occ, EVENT_CFs, ALL_CFs) {
  
  library(tidyr)
  library(dplyr)
  library(ThreadNet)
  library(ngram)
  library(lubridate)
  library(stringr)
  library(data.table)
  
  
  # get the first occurrence from each thread -- this will be used for the stuff that never changes, like clinic
  firstOcc=occ[occ$seqNum==1,]
  
  # always used threadNum for the TN
  TN = 'threadNum'
  
  # Add  columns for combinations of CFs if needed
  # don't need to add new column for the thread_CFs.  That one has to be here.
  # thread_col = newColName(THREAD_CFs)
  # TN = thread_col
  # if  (!(thread_col  %in% colnames(occ)))  {  occ = combineContextFactors(occ,THREAD_CFs,thread_col) }
  
  new_event_col = newColName(EVENT_CFs)
  if  (!(new_event_col  %in% colnames(occ)))  {  occ = combineContextFactors(occ,EVENT_CFs,new_event_col) }
  
  all_cf_col = newColName(ALL_CFs)
  if  (!(all_cf_col  %in% colnames(occ)))  {  occ = combineContextFactors(occ,ALL_CFs,all_cf_col) } 
  
  # set key on the data.table for the threadNum
  occ=as.data.table(occ)
  setkeyv(occ, TN)
  
  # get the list of buckets
  bucket_list <- sort( unique(occ[[TN]]) )
  
  # get the number of buckets
  N = length(bucket_list)
  print(paste0('Number of buckets=', N  ))
  
  # allocate the data table with N rows, then for speed, we use := to update it.
  ACHR = data.table(threadNum =  integer(N),
                    Phase  = character(N), 
                    NEvents = integer(N),
                    ThreadDuration =double(N),
                    VisitDuration= double(N),
                    wait_time1 = double(N),
                    wait_time2 = double(N),
                    LOS =double(N),
                    NetComplexity=double(N),
                    Nodes=double(N),
                    Edges=double(N),
                    CompressRatio = double(N),
                    Entropy = double(N),
                    NumProcedures = double(N),
                    NumDiagnoses = double(N),
                    ALL_CF_count =  integer(N),
                    ALL_CF_entropy = double(N),
                    CF_Alignment = double(N) )

    # Now add columns for the CF counts.  
    for (cf in ALL_CFs){
      ACHR[, paste0(cf,"_count"):= double(N)]
      # ACHR[, paste0(cf,"_compression"):= double(N)]
      # ACHR[, paste0(cf,"_entropy"):= double(N)]
    }
  
  # update data table with results - only compute the ones that require looking at the whole thread
  for (b in 1:N){
    
    if (b %% 1000 ==0) print(paste0('Thread count =  ',b))
    
    # select a subset 
    df= occ[get(TN)==bucket_list[b]]
    
    # make sure it is sorted by timestamp
    df=df[order(df$tStamp),]
    
    # get the network -- only if there are enough rows...
    if (nrow(df)>2)   n = threads_to_network_original(df,TN, new_event_col) 
    else  n = list(edgeDF=t(c(0)),nodeDF=t(c(0)))
    
    # compute each parameter and update the table
    ACHR[b,threadNum := bucket_list[b] ]
    ACHR[b,Phase :=  compute_phase(df$tStamp[1])]
    ACHR[b,NEvents := nrow(df) ]
    ACHR[b,ThreadDuration:= compute_thread_duration(df) ] 
    ACHR[b,VisitDuration:= compute_visit_duration(df) ] 
    ACHR[b, wait_time1 := compute_wait_time1(df) ]
    ACHR[b, wait_time2 := compute_wait_time2(df) ]
    ACHR[b, NetComplexity:=estimate_network_complexity( n ) ]
    ACHR[b, Nodes:=nrow(n$nodeDF) ]
    ACHR[b,Edges:=nrow(n$edgeDF) ]
    ACHR[b,CompressRatio := compression_index(df,new_event_col) ]
  #  ACHR[b,Entropy := compute_entropy(table(df[[new_event_col]])[table(df[[new_event_col]])>0]) ]
    ACHR[b,Entropy := compute_graph_entropy( df[[new_event_col]]) ]
    ACHR[b,NumProcedures := count_procedures(df$Proc[1]) ]
    ACHR[b,NumDiagnoses := count_diagnoses(df$Diag[1]) ]
    ACHR[b,ALL_CF_count := length(unique(df[[all_cf_col]])) ]
    # ACHR[b,ALL_CF_entropy := compute_entropy(table(df[[all_cf_col]])[table(df[[all_cf_col]])>0]) ]
    ACHR[b,ALL_CF_entropy := compute_graph_entropy( df[[all_cf_col]])  ] 
    
    # Count the unique elements in each cf  
    for (cf in ALL_CFs){  ACHR[b, paste0(cf,"_count") :=  length(unique(df[[cf]])) ] }
  } 
  
  # Compute the alignment of the context factors
  ACHR$CF_Alignment =  as.numeric( as.character(ACHR$Action_count)) / as.numeric( as.character(ACHR$ALL_CF_count ))
  
  # convert level of service to 1-5 integer 
  ACHR$LOS = convert_LOS( firstOcc$LOS_CPT )
  
  # Merge the results with the first row from each thread
  print('Merging results...')
  Thrds=merge(x=ACHR, y=firstOcc, by.x='threadNum', by.y ='threadNum' ,all=TRUE) 
  
  save_file_name = paste0(paste('Thrds',TN,new_event_col,sep='+'), '.Rdata')
  save(Thrds, file=save_file_name)
  
  print(paste('Saved ', nrow(Thrds), " records in ",save_file_name))
  
  return(Thrds) 
}

#################################################################################
# streamlined version for threads that are PART OF A VISIT.  So you can merge with visits
# to get the visit context.  Just compute the minimum for speed.
ACHR_batch_visit_role_threads <- function(occ, EVENT_CFs, ALL_CFs, visits) {
  
  library(tidyr)
  library(dplyr)
  library(ThreadNet)
  library(ngram)
  library(lubridate)
  library(stringr)
  library(data.table)
  
  
  # get the first occurrence from each thread -- this will be used for the stuff that never changes, like clinic
 # firstOcc=occ[occ$seqNum==1,]
  
  # always used threadNum for the TN
  TN = 'threadNum'
  
  # Add  columns for combinations of CFs if needed
  # don't need to add new column for the thread_CFs.  That one has to be here.
  # thread_col = newColName(THREAD_CFs)
  # TN = thread_col
  # if  (!(thread_col  %in% colnames(occ)))  {  occ = combineContextFactors(occ,THREAD_CFs,thread_col) }
  
  new_event_col = newColName(EVENT_CFs)
  if  (!(new_event_col  %in% colnames(occ)))  {  occ = combineContextFactors(occ,EVENT_CFs,new_event_col) }
  
  all_cf_col = newColName(ALL_CFs)
  if  (!(all_cf_col  %in% colnames(occ)))  {  occ = combineContextFactors(occ,ALL_CFs,all_cf_col) } 
  
  # set key on the data.table for the threadNum
  occ=as.data.table(occ)
  setkeyv(occ, TN)
  
  # get the list of buckets
  bucket_list <- sort( unique(occ[[TN]]) )
  
  # get the number of buckets
  N = length(bucket_list)
  print(paste0('Number of buckets=', N  ))
  
  # allocate the data table with N rows, then for speed, we use := to update it.
  ACHR = data.table(Visit_ID  = character(N), 
                    Role_VR = character(N), 
                    Role_ID_VR =  character(N), 
                    threadNumVR =  integer(N),
                    NEventsVR = integer(N),
                    threadStartVR =  character(N), 
                    ThreadDurationVR =double(N),
                    LOS =double(N),
                    NetComplexityVR=double(N),
                    NodesVR=double(N),
                    EdgesVR=double(N),
                    CompressRatioVR = double(N),
                    EntropyVR = double(N),
                    ALL_CF_countVR =  integer(N),
                    ALL_CF_entropyVR = double(N) )
  
  # Now add columns for the CF counts.  
  for (cf in ALL_CFs){
    ACHR[, paste0(cf,"_countVR"):= double(N)]
    # ACHR[, paste0(cf,"_compression"):= double(N)]
    # ACHR[, paste0(cf,"_entropy"):= double(N)]
  }
  
  # update data table with results - only compute the ones that require looking at the whole thread
  for (b in 1:N){
    
    if (b %% 1000 ==0) print(paste0('Thread count =  ',b))
    
    # select a subset 
    df= occ[get(TN)==bucket_list[b]]
    
    # get the network -- only if there are enough rows...
    if (nrow(df)>2)   n = threads_to_network_original(df,TN, new_event_col) 
    else  n = list(edgeDF=t(c(0)),nodeDF=t(c(0)))
    
    # compute each parameter and update the table
    ACHR[b,Visit_ID := as.character(df$Visit_ID[1]) ]
    ACHR[b,Role_VR:= as.character(df$Role[1]) ]
    ACHR[b,Role_ID_VR:= as.character(df$Role_ID[1]) ]
    ACHR[b,threadNumVR := bucket_list[b] ]
    ACHR[b,NEventsVR := nrow(df) ]
    ACHR[b,threadStartVR := as.character(min(lubridate::ymd_hms(df$tStamp))) ]
    ACHR[b,ThreadDurationVR:= difftime(max(lubridate::ymd_hms(df$tStamp)),  min(lubridate::ymd_hms(df$tStamp)), units='hours' ) ] 
    ACHR[b, NetComplexityVR:=estimate_network_complexity( n ) ]
    ACHR[b, NodesVR:=nrow(n$nodeDF) ]
    ACHR[b,EdgesVR:=nrow(n$edgeDF) ]
    ACHR[b,CompressRatioVR := compression_index(df,new_event_col) ]
    # ACHR[b,EntropyVR := compute_entropy(table(df[[new_event_col]])[table(df[[new_event_col]])>0]) ]
    ACHR[b,EntropyVR := compute_graph_entropy( df[[new_event_col]])  ]
    
    ACHR[b,ALL_CF_countVR := length(unique(df[[all_cf_col]])) ]
    # ACHR[b,ALL_CF_entropyVR := compute_entropy(table(df[[all_cf_col]])[table(df[[all_cf_col]])>0]) ] 
    ACHR[b,ALL_CF_entropyVR := compute_graph_entropy( df[[all_cf_col]])  ] 
    
    # Count the unique elements in each cf  
    for (cf in ALL_CFs){  ACHR[b, paste0(cf,"_countVR") :=  length(unique(df[[cf]])) ] }
  } 
  
  
  # Merge the results with the first row from each thread
  print('Merging results...')
  VRThrds=merge(x=ACHR, y=visits, by.x='Visit_ID', by.y ='Visit_ID' ,all=TRUE) 
  
  save_file_name = paste0(paste('VRThrds',TN,new_event_col,sep='+'), '.Rdata')
  save(VRThrds, file=save_file_name)
  
  print(paste('Saved ', nrow(VRThrds), " records in ",save_file_name))
  
  return(VRThrds) 
}

#################################################################################
# streamlined version for adding or recomputing columns that describe a visit.  
# So you can merge with visits to get rest of the columns.  Just compute the minimum for speed.
# Can re-write this as needed... 
ACHR_batch_visit_add_columns <- function(occ, EVENT_CFs, visits) {
  
  library(tidyr)
  library(dplyr)
  library(ThreadNet)
  library(ngram)
  library(lubridate)
  library(stringr)
  library(data.table)
  
  
  # get the first occurrence from each thread -- this will be used for the stuff that never changes, like clinic
  # firstOcc=occ[occ$seqNum==1,]
  
  # always used threadNum for the TN
  TN = 'threadNum'
  
  # Add  columns for combinations of CFs if needed
  # don't need to add new column for the thread_CFs.  That one has to be here.
  # thread_col = newColName(THREAD_CFs)
  # TN = thread_col
  # if  (!(thread_col  %in% colnames(occ)))  {  occ = combineContextFactors(occ,THREAD_CFs,thread_col) }
  
  new_event_col = newColName(EVENT_CFs)
  if  (!(new_event_col  %in% colnames(occ)))  {  occ = combineContextFactors(occ,EVENT_CFs,new_event_col) }
  
 
  # set key on the data.table for the threadNum
  occ=as.data.table(occ)
  setkeyv(occ, TN)
  
  # get the list of buckets
  bucket_list <- sort( unique(occ[[TN]]) )
  
  # get the number of buckets
  N = length(bucket_list)
  print(paste0('Number of buckets=', N  ))
  
  # allocate the data table with N rows, then for speed, we use := to update it.
  ACHR = data.table(Visit_ID  = character(N), 
                    new_wait_time1 = double(N),
                    new_wait_time2 = double(N) )
                    
                    # Role_VR = character(N), 
                    # Role_ID_VR =  character(N), 
                    # threadNumVR =  integer(N),
                    # NEventsVR = integer(N),
                    # threadStartVR =  character(N), 
                    # ThreadDurationVR =double(N),
                    # NetComplexityVR=double(N),
                    # NodesVR=double(N),
                    # EdgesVR=double(N),
                    # CompressRatioVR = double(N),
                    # EntropyVR = double(N),
                    # ALL_CF_countVR =  integer(N),
                    # ALL_CF_entropyVR = double(N) )
  
  
  # update data table with results - only compute the ones that require looking at the whole thread
  for (b in 1:N){
    
    if (b %% 1000 ==0) print(paste0('Thread count =  ',b))
    
    # select a subset 
    df= occ[get(TN)==bucket_list[b]]
    
    # get the network -- only if there are enough rows...
    # if (nrow(df)>2)   n = threads_to_network_original(df,TN, new_event_col) 
    # else  n = list(edgeDF=t(c(0)),nodeDF=t(c(0)))
    # 
    # compute each parameter and update the table
    ACHR[b,Visit_ID := as.character(df$Visit_ID[1]) ]
    ACHR[b, new_wait_time1 := compute_wait_time1(df) ]
    ACHR[b, new_wait_time2 := compute_wait_time2(df) ]
  } 
  
  
  # Compute the alignment of the context factors
  #  ACHR$CF_AlignmentVR =  as.numeric( as.character(ACHR$Action_countVR)) / as.numeric( as.character(ACHR$ALL_CF_countVR ))
  
  # Merge the results with the first row from each thread
  print('Merging results...')
 # VRThrds=merge(x=ACHR, y=visits, by.x='Visit_ID', by.y ='Visit_ID' ,all=TRUE) 
  VRThrds=merge(x=ACHR, y=visits, by.x='Visit_ID', by.y ='Visit_ID') 
  
  save_file_name = paste0(paste('VRThrds',TN,new_event_col,sep='+'), '.Rdata')
  save(VRThrds, file=save_file_name)
  
  print(paste('Saved ', nrow(VRThrds), " records in ",save_file_name))
  
  return(VRThrds) 
}
##################################################################################
# Helper functions
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
  
  w1= grep('CHECKIN_TIME',df$Action)[1]
  w2= grep('VITALS',df$Action)[1]
  
  return( max(0, as.numeric(difftime( lubridate::ymd_hms(df$tStamp[w2]),lubridate::ymd_hms(df$tStamp[w1]),  units='hours' ) )))
   
}

# Second wait time... from  vitals until they see a doctor or other medical person
# NEW VERSION -- look for "AC_VISIT_NAVIGATOR" by a different role than "VITALS"
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
    # work1 = df$Workstation[w1]
    role1 = df$Role[w1]
    
    # get set of actions at the same workstation or different workstations that happen later in the visit
    # Truncate the visit up to w1
    restdf = df[(w1+1):nrow(df),]
    
    # get the actions performed by other roles in the rest of the visit
    differentRole =  restdf[restdf$Role != role1,]
    
    # Now look  for first instance of  AC_VISIT_NAVIGATOR that is by one of the other roles
       # first_AC_NAV =  grep('AC_VISIT_NAVIGATOR',differentRole$Action)[1]
       # t2 = lubridate::ymd_hms(differentRole$tStamp[first_AC_NAV])
       
     # Try first instance of different role. period.  
      t2 = lubridate::ymd_hms(differentRole$tStamp[1])
      
      
    # return the difference in time from t1 to t2
    return( max(0, as.numeric(difftime( t2, t1,  units='hours' ) )))
  }
}

# Second wait time... from  vitals until they see a doctor or other medical person
# compute_wait_time2 <- function(df){
#   
#   #  find the first  occurrence of MR_VN_VITALS. Get the time t1
#  # w1= grep('MR_VN_VITALS',df$Action)[1]
#   w1= grep('VITALS',df$Action)[1]
#   
#   if (is.na(w1))
#     return(NA)
#   else
#   {
#   t1 = lubridate::ymd_hms(df$tStamp[w1])
#   
#   
#   # get the workstation and role used for the vitals at t1
#   work1 = df$Workstation[w1]
#   role1 = df$Role[w1]
#   
#   # get set of actions at the same workstation or different workstations that happen later in the visit
#   # Truncate the visit up to w1
#   restdf = df[w1:nrow(df),]
# 
#   # then get the set of action at the same workstation or different workstations (in that visit)
#   sameWS =  restdf[restdf$Workstation == work1,]
#   differentWS =  restdf[restdf$Workstation != work1,]
#   
#   # get list of other roles
#   otherRoles = paste(setdiff(unique(restdf$Role),c('Technician','Technologist','Staff','Unknown')),
#                     collapse = '|')
#   
#   # Now find the marker for when they see a physician
#   # procedure is differenrt in BRKPT and HHPOB clinics
#   if (df$Clinic[1]  %in% c('BRKPT','HHPOB') )
#     { # LPN brings a different workstation on wheels at these clinics  
#     firstLPN =  grep('LPN',differentWS$Role)[1] 
#     t2 = lubridate::ymd_hms(differentWS$tStamp[firstLPN])
#     
#     # Other roles use the same workstation
#     firstDoc =  grep('Physician|Resident|Registered', sameWS$Role)[1] 
#     t2 = lubridate::ymd_hms(sameWS$tStamp[firstDoc])
#     }
#   else 
#   {  # a different role, using the same workstation 
#     w2 = max(w1, grep(otherRoles, sameWS$Role) )[1]
#     t2 = lubridate::ymd_hms(sameWS$tStamp[w2])
#     }
#     
#   # return the difference in time from t1 to t2
#   return( difftime( t2, t1,  units='hours' ) )
#   }
# }


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

count_daily_procedures <- function( df ) {
  
  # get the string of procedures for that day
  allproc  = unlist(df[df$seqNum==1, 'Proc'])
  
  if (length(allproc)>0)  
    return( sum(sapply(allproc, function(p) {count_procedures(p)} )) )
  else return(0)
  
  
}

# convert the LOS_CPT codes into a five point numeric scale
# operate on whole column to create new column
convert_LOS <- function(los_cpt){
  
  # get the 5th character
  s = sapply(los_cpt, function(cpt)  { substring(as.character(cpt),5,5) } )
  
  los = as.numeric(s)
  
  # convert 9 and NA to 1
  los[los==9]=NA


return(los)  
  
}

###################################################################
###################################################################
# This is function is set up to aggregate the occurrences among collections of visit -- typically clinic_days
# occ = pre-processed threaded occurrences
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
  N =  length(bucket_list) 
  
  # print the number of buckets
  print(paste0('Number of buckets=', N  ))
  
  # Add  columns for combinations of CFs if needed
  new_event_col = newColName(EVENT_CFs)
  if  (!(new_event_col  %in% colnames(occ)))  {  occ = combineContextFactors(occ,EVENT_CFs,new_event_col) }
  
  all_cf_col = newColName(ALL_CFs)
  if  (!(all_cf_col  %in% colnames(occ)))  {  occ = combineContextFactors(occ,ALL_CFs,all_cf_col) } 
  
  # this will speed up retrieving the subsets
  setkeyv(occ, TN)
  
  # allocate the data table with N rows, then for speed, we use := to update it.
  cds = data.table(Clinic_ymd  = character(N), 
                    Clinic = character(N),
                    ymd = character(N),
                    Phase =  character(N), 
                    NEvents = integer(N),
                    ClinicDayStart =  character(N), 
                    ClinicDayDuration =double(N),
                    NumVisits =  integer(N),
                    NumUniqueDiagnosisGroups = integer(N), 
                    NumPhysicians = integer(N),
                    TotalStaff = integer(N),
                    NetComplexity=double(N),
                    Nodes=double(N),
                    Edges=double(N),
                    CompressRatio = double(N),
                    Entropy = double(N),
                    ALL_CF_count =  integer(N),
                    ALL_CF_entropy = double(N),
                    CF_Alignment  = double(N))
  
  # Now add columns for the CF counts.  
  for (cf in ALL_CFs){
    cds[, paste0(cf,"_countVR"):= double(N)]
    # cds[, paste0(cf,"_compression"):= double(N)]
    # cds[, paste0(cf,"_entropy"):= double(N)]
  }
  
  
  # make data frame with results
 for (b in 1:N){
   
   
   if (b %% 100 ==0) print(paste0('Thread count =  ',b))
   
   
   # select a subset of occurrences for the bucket
   df = occ[ occ[[TN]] == bucket_list[b]  ] 
  
  
  # get the network -- only if there are enough rows...
  if (nrow(df)>2)   n = threads_to_network_original(df,TN, new_event_col) 
  else  n = list(edgeDF=t(c(0)),nodeDF=t(c(0)))
  
  # compute each parameter and put them in a vector
  
  cds[b,Clinic_ymd :=  bucket_list[b] ]
  cds[b,Clinic := as.character( df$Clinic[1] ) ]
  cds[b,ymd := as.character( df$ymd[1] ) ]
  cds[b,Phase  :=   compute_phase(df$tStamp[1]) ]
  cds[b,threadNum  :=  as.numeric(df[1,'threadNum']) ]
  cds[b,NEvents  :=  nrow(df) ]
  cds[b,ClinicDayStart :=  as.character(lubridate::ymd_hms(df$tStamp[1]))  ]
  cds[b,ClinicDayDuration :=  difftime(max(lubridate::ymd_hms(df$tStamp)),  min(lubridate::ymd_hms(df$tStamp)), units := 'hours' ) ]
  cds[b,NumVisits  :=   length(unique( df$Visit_ID )) ]
  cds[b,NumUniqueDiagnosisGroups  :=  length(unique( df$Diagnosis_Group )) ] 
  cds[b,NumPhysicians  :=  length(unique( df$Physician )) ]
  cds[b,TotalStaff  :=  length(unique( df$Role_ID )) ]
  cds[b, NetComplexity := estimate_network_complexity( n ) ]
  cds[b,Nodes := nrow(n$nodeDF) ]
  cds[b,Edges := nrow(n$edgeDF) ]
  cds[b,CompressRatio  :=  compression_index(df,new_event_col) ]
  # cds[b,Entropy  :=  compute_entropy(table(df[[new_event_col]])[table(df[[new_event_col]])>0]) ]
  cds[b,Entropy  :=  compute_graph_entropy( df[[new_event_col]])  ]
  
  cds[b,NumProceduresPerDay  :=  count_daily_procedures(df) ]
  cds[b,CF_Alignment  :=  1 ]    # make placeholder, but compute below
  cds[b,ALL_CF_count  :=  length(unique(df[[all_cf_col]])) ]
  # cds[b,ALL_CF_entropy  :=  compute_entropy(table(df[[all_cf_col]])[table(df[[all_cf_col]])>0]) ] 
  cds[b,ALL_CF_entropy  :=  compute_graph_entropy( df[[all_cf_col]])  ] 
  
  # Count the unique elements in each cf  
  for (cf in ALL_CFs){  cds[b, paste0(cf,"_countVR") :=  length(unique(df[[cf]])) ] }
  
  
} 
  
  # Compute the alignment of the context factors
  cds$CF_Alignment =  as.numeric( as.character(cds$Action_count)) / as.numeric( as.character(cds$ALL_CF_count ))
  
  save(cds, file=paste0(paste('cds',TN,new_event_col,sep='+'), '.Rdata'))
  
  return(cds) 
  
}

compute_graph_entropy <- function(s){
  
  # guard against insufficient data...
  if (length(s) <2) return(0)
  
  # convert s into text vector
  text_vector =  concatenate(s) 
  
  # get the 2-grams in the sequence s. 'prop' is the proportion of each edge. It sums to 1. 
  p = get.phrasetable(ngram(text_vector,2))[['prop']]
  
  # return Shannon entropy
  return(-sum(p*log(p)))
  
}

