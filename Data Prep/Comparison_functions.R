##########################################################################################################
# THREADNET:  Batch processing for larger data sets

# (c) 2017 Michigan State University. This software may be used according to the terms provided in the
# GNU General Public License (GPL-3.0) https://opensource.org/licenses/GPL-3.0?
# Absolutely no warranty!
##########################################################################################################


library(dplyr)

get_date_subsets <-function(c) {
b1 <<- new_occ_VR %>% filter(Clinic==eval(c) & as.Date(tStamp) > '2016-01-01' & as.Date(tStamp) < '2016-06-06')
b2 <<- new_occ_VR %>% filter(Clinic==eval(c) & as.Date(tStamp) > '2016-06-8' & as.Date(tStamp) < '2016-09-01')
b3 <<- new_occ_VR %>% filter(Clinic==eval(c) & as.Date(tStamp) > '2016-09-02' & as.Date(tStamp) < '2017-04-15')
b4 <<- new_occ_VR %>% filter(Clinic==eval(c) & as.Date(tStamp) > '2017-04-15' & as.Date(tStamp) < '2017-09-01')
b5 <<- new_occ_VR %>% filter(Clinic==eval(c) & as.Date(tStamp) > '2017-09-02' & as.Date(tStamp) < '2017-12-31')

rc1 <<- new_occ_VR %>% filter(Clinic=='REDCK' & as.Date(tStamp) > '2016-01-01' & as.Date(tStamp) < '2016-06-06')
rc2 <<- new_occ_VR %>% filter(Clinic=='REDCK' & as.Date(tStamp) > '2016-06-8' & as.Date(tStamp) < '2016-09-01')
rc3 <<- new_occ_VR %>% filter(Clinic=='REDCK' & as.Date(tStamp) > '2016-09-02' & as.Date(tStamp) < '2017-04-15')
rc4 <<- new_occ_VR %>% filter(Clinic=='REDCK' & as.Date(tStamp) > '2017-04-15' & as.Date(tStamp) < '2017-09-01')
rc5 <<- new_occ_VR %>% filter(Clinic=='REDCK' & as.Date(tStamp) > '2017-09-02' & as.Date(tStamp) < '2017-12-31')

smh1 <<- new_occ_VR %>% filter(Clinic=='SMH' & as.Date(tStamp) > '2016-01-01' & as.Date(tStamp) < '2016-06-06')
smh2 <<- new_occ_VR %>% filter(Clinic=='SMH' & as.Date(tStamp) > '2016-06-8' & as.Date(tStamp) < '2016-09-01')
smh3 <<- new_occ_VR %>% filter(Clinic=='SMH' & as.Date(tStamp) > '2016-09-02' & as.Date(tStamp) < '2017-04-15')
smh4 <<- new_occ_VR %>% filter(Clinic=='SMH' & as.Date(tStamp) > '2017-04-15' & as.Date(tStamp) < '2017-09-01')
smh5 <<- new_occ_VR %>% filter(Clinic=='SMH' & as.Date(tStamp) > '2017-09-02' & as.Date(tStamp) < '2017-12-31')

hh1 <<- new_occ_VR %>% filter(Clinic=='HHPOB' & as.Date(tStamp) > '2016-01-01' & as.Date(tStamp) < '2016-06-06')
hh2 <<- new_occ_VR %>% filter(Clinic=='HHPOB' & as.Date(tStamp) > '2016-06-8' & as.Date(tStamp) < '2016-09-01')
hh3 <<- new_occ_VR %>% filter(Clinic=='HHPOB' & as.Date(tStamp) > '2016-09-02' & as.Date(tStamp) < '2017-04-15')
hh4 <<- new_occ_VR %>% filter(Clinic=='HHPOB' & as.Date(tStamp) > '2017-04-15' & as.Date(tStamp) < '2017-09-01')
hh5 <<- new_occ_VR %>% filter(Clinic=='HHPOB' & as.Date(tStamp) > '2017-09-02' & as.Date(tStamp) < '2017-12-31')

}

# c is for the clinic 
# m is the metric
compare_periods <-function(c,m){

  # b1 <- new_occ_VR %>% filter(Clinic==c & as.Date(tStamp) > '2016-01-01' & as.Date(tStamp) < '2016-06-06')
  # b2 <- new_occ_VR %>% filter(Clinic==c & as.Date(tStamp) > '2016-06-8' & as.Date(tStamp) < '2016-09-01')
  # b3 <- new_occ_VR %>% filter(Clinic==c & as.Date(tStamp) > '2016-09-02' & as.Date(tStamp) < '2017-04-15')
  # b4 <- new_occ_VR %>% filter(Clinic==c & as.Date(tStamp) > '2017-04-15' & as.Date(tStamp) < '2017-09-01')
  # b5 <- new_occ_VR %>% filter(Clinic==c & as.Date(tStamp) > '2017-09-02' & as.Date(tStamp) < '2017-12-31')
  
  
  a1=routine_metric(b1,'Action',m )
  a2=routine_metric(b2,'Action',m )
  a3=routine_metric(b3,'Action',m )
  a4=routine_metric(b4,'Action',m )
  a5=routine_metric(b5,'Action',m )
  
  r1=routine_metric(b1,'Role',m )
  r2=routine_metric(b2,'Role',m )
  r3=routine_metric(b3,'Role',m )
  r4=routine_metric(b4,'Role',m )
  r5=routine_metric(b5,'Role',m )
  
  ws1=routine_metric(b1,'Workstation',m )
  ws2=routine_metric(b2,'Workstation',m )
  ws3=routine_metric(b3,'Workstation',m )
  ws4=routine_metric(b4,'Workstation',m )
  ws5=routine_metric(b5,'Workstation',m )
  
  rc_a1=routine_metric(rc1,'Action',m )
  rc_a2=routine_metric(rc2,'Action',m )
  rc_a3=routine_metric(rc3,'Action',m )
  rc_a4=routine_metric(rc4,'Action',m )
  rc_a5=routine_metric(rc5,'Action',m )
  
  smh_a1=routine_metric(smh1,'Action',m )
  smh_a2=routine_metric(smh2,'Action',m )
  smh_a3=routine_metric(smh3,'Action',m )
  smh_a4=routine_metric(smh4,'Action',m )
  smh_a5=routine_metric(smh5,'Action',m )
  
  hh_a1=routine_metric(hh1,'Action',m )
  hh_a2=routine_metric(hh2,'Action',m )
  hh_a3=routine_metric(hh3,'Action',m )
  hh_a4=routine_metric(hh4,'Action',m )
  hh_a5=routine_metric(hh5,'Action',m )
  
  lex=data.frame(c(a1,a2,a3,a4,a5),
                 c(hh_a1,hh_a2,hh_a3,hh_a4,hh_a5),
                 c(rc_a1,rc_a2,rc_a3,rc_a4,rc_a5),
                 c(smh_a1,smh_a2,smh_a3,smh_a4,smh_a5)
                 )
  
  ws1=ws1[ws1>0]/nrow(b1)
  ws2=ws2[ws2>0]/nrow(b2)
  ws3=ws3[ws3>0]/nrow(b3)

}


# o is the list of occurrences
# m is the metric we want to see
routine_metric <- function(o, cf, m){
  
  switch(m,
         len = nrow(o),
         lex_size = length(unique(o[[cf]])),
         first_order = table(o[[cf]]),
         first_order_sort = sort( table(o[[cf]]), decreasing = T)
         
         )
  
}

get_range <- function(o, start_date, end_date){
  
  o %>% filter(as.Date(tStamp) > start_date & as.Date(tStamp) < end_date)
  
}
