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


# Two  ways to aggregate the occurrences 
# 1) By thread (usually visit)
# 2) By Clinic_day
  



###################################################################
###################################################################
# This is function is set up to aggregate the occurrences within threads -- typically visits
# occ = pre-processed threaded occurrences
# TN = threadNum in most cases
# CFs can be chosen -- they don't have to match the POV

ACHR_batch_threads <- function(occ,TN, EVENT_CFs, ALL_CFs) {
  
  library(tidyr)
  library(dplyr)
  library(ThreadNet)
  library(ngram)
  library(lubridate)
  library(stringr)
  
  # don't remember why  this is here...
  new_event_col = newColName(EVENT_CFs)
  if  (!(new_event_col  %in% colnames(occ)))  {  occ = combineContextFactors(occ,EVENT_CFs,new_event_col) }
  
  # get the list of buckets
  bucket_list <- make_buckets_1(occ, TN)
  
  # make data frame with results
 Thrds = data.frame( t(sapply( bucket_list, 
                           function(b){
                             
                             # select a subset 
                             df = occ[ occ[[TN]] ==b , ] 
                             
                            # Only run for visits with more than two occurrences
                        
                             # get the network
                             n = threads_to_network_original(df,TN, new_event_col)
                             
                             # compute each parameter and put them in a vector
                             c(
                               bucket=b,
                               NEvents = nrow(df),
                               ThreadStart= as.character( df[1,'tStamp'] ),
                               ThreadStartInt = df[1,'tStamp'],
                               ThreadDuration= difftime(max(lubridate::ymd_hms(df$tStamp)),  min(lubridate::ymd_hms(df$tStamp)), units='hours' ),
                               NetComplexity=estimate_network_complexity( n ),
                               Nodes=nrow(n$nodeDF),
                               Edges=nrow(n$edgeDF),
                               CompressRatio = compression_index(df,new_event_col),
                               Entropy = compute_entropy(table(df[[new_event_col]])[table(df[[new_event_col]])>0]),
                               NumProcedures = count_procedures(df$Proc[1]),
                               NumDiagnoses = count_diagnoses(df$Diag[1]),
                               wait_time = compute_wait_time(df),
                             #   Visit_number = integer(N),
                             #   LOC_CPT =  character(N),
                                Visit_ID  = df[1,'Visit_ID'],
                                Subject_ID  = df[1,'Subject_ID'],
                                Clinic = df[1,'Clinic'],
                                ymd = df[1,'ymd'],
                                Clinic_ymd = df[1,'Clinic_ymd'],
                                Proc = df[1,'Proc'],
                                Diagnosis = df[1,'Diag'],   # diag in the raw data
                                Diagnosis_group  = df[1,'Diagnosis_Group'],
                                Physician  = df[1,'Physician'],
                                Weekday  = df[1,'Weekday'],
                                Month  = df[1,'Month'],
                                Phase =  compute_phase(df$tStamp[1]),
                                CF_Alignment = compute_alignment(df,TN, EVENT_CFs, ALL_CFs ),
                             
                             # need to figure out how to name these
                             sapply(ALL_CFs, function(cf){
                               c(
                                 assign( paste0(cf,"_count"), length(unique(df[[cf]]))),
                                 assign( paste0(cf,"_entropy"), compute_entropy(table(df[[cf]])[table(df[[cf]])>0])) 
                               ) })
                             
                             )
                             
  } )))
 
 save(Thrds, file=paste0(paste('Thrds',TN,new_event_col,sep='+'), '.Rdata'))
 
 return(Thrds) 
  
}

# Each bucket is a list of thread numbers that can be used to subset the list of occurrences
make_buckets_1 <- function(o, criteria){
  
  return( unique(o[[criteria]]) )
  
}

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

# alignment is a measure of how much the number of nodes is increased by adding context factors
# always base this on the "Action" -- how many more "actions" are there? 
# could also do this based on entropy, or write the function to use the values that were already computed
compute_alignment <- function(df,TN, EVENT_CFs, ALL_CFs ){
  
  # first get numerator
  n= length(unique(df$Action))
  
  # Now get the total number with all  of the CFs used to define the events
  m= length(unique( df[[newColName(EVENT_CFs)]]))

  return(n/m)
}


compute_wait_time <- function(df){
  
  #  find the first  occurrence of the chief complain.  Seems to correspond closely to wait time
  #  Make sure it's not in the first 5 occurrences, and make sure to  return a  value of at least 1 if it never occurs
  w= grep('MR_VN_CHIEF_COMPLAINT',df$Action)[1]
  
  return( difftime( lubridate::ymd_hms(df$tStamp[w]),lubridate::ymd_hms(df$tStamp[1]),  units='hours' ) )
   
}

count_procedures <- function( p ){
  
  # get the overall number of items
  total_num=length(grep('#@#', p ))
  
  # adjust for the visit codes that are not actual procedures
  v=c('99211','99212','99213','99214','99215','99201','99202','99203','99204','99205')
  
  # just look at the first occurrence in the df
  num_visit_codes  = sum(sapply(v,function(x){length(grep( x, p ))}))
  
  # limit  to non-negative 
  total_num = max(0, total_num - num_visit_codes) 
  
  return( total_num  )
}

# This one  is easy
count_diagnoses <- function(d){
  return(length(grep('#@#', d )))
}


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
ACHR_batch_clinic_days <- function(occ,TN, EVENT_CFs, ALL_CFs) {
  
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




# 
# make_box_plots <- function(){
# ggboxplot(ACHR_test[NEvents>100 & Clinic=='DRH'], x = "VisitMonth", y = "NetComplexity",
#           color = "VisitDay",
#           ylab = "Complexity", xlab = "Month (DRH)")
# }
# 





get_timeScale <- function(){'hr'}

