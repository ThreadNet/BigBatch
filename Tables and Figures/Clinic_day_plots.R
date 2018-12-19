##########################################################################################################
# THREADNET:  Batch processing for larger data sets

# (c) 2017 Michigan State University. This software may be used according to the terms provided in the
# GNU General Public License (GPL-3.0) https://opensource.org/licenses/GPL-3.0?
# Absolutely no warranty!
##########################################################################################################

# These are just examples -- need to  work out tables for particular purposes.

# here are some example plots 
# much better to use ggplot

# cds is for clinic_days

more_plots <- function(){
  
  # More plots  -- these operate on the merged clinic-day trajetories
  plot(cds$ymd[cds$Clinic.x=='BATAVIA'],cds$cosine_dist_from_first[cds$Clinic.x=='BATAVIA'])
  plot(cds$ymd[cds$Clinic.x=='BATAVIA'],cds$NumVisits[cds$Clinic.x=='BATAVIA'])
  plot(cds$ymd[cds$Clinic.x=='BATAVIA'],cds$NumUniqueProcedures[cds$Clinic.x=='BATAVIA'])
  plot(cds$ymd[cds$Clinic.x=='BATAVIA'],cds$NumUniqueDiagnosisGroups[cds$Clinic.x=='BATAVIA'])
  plot(cds$ymd[cds$Clinic.x=='BATAVIA'],cds$cosine_dist_from_first[cds$Clinic.x=='BATAVIA'])
  
}


# This is the basic plot to compare all the clinic trajectories over time
plots_for_papers <- function(){
  
  library(ggplot2)
  
  # compare the complexity of the five clinics
  plot(cds$Clinic.x,cds$A_NetComplexity)
  
  # Clinic change over time, with and without smoothing
  ggplot(data = cds, aes(x = ymd, y = Dist_from_reference, group=Clinic)) + geom_line(aes(color=Clinic))
  ggplot(data = cds, aes(x = ymd, y = rollmean(Dist_from_reference,5,na.pad=TRUE), group=Clinic)) + geom_line(aes(color=Clinic))
  
  # Dec 6th version 
  ggplot(data = cdt, aes(x = ymd.x, y = Dist_from_reference, group=Clinic.x)) + 
    geom_line(aes(color=Clinic.x)) + theme(axis.text.x=element_blank(), axis.text.y=element_blank())

    ggplot(data = cdt, aes(x = ymd.x, y = rollmean( Dist_from_reference, 5,na.pad=TRUE), group=Clinic.x)) + 
    geom_line(aes(color=Clinic.x)) + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
    
    ggplot(data = cdt_435[,1:6], aes(x = ymd, y = rollmean( Dist_from_reference, 10,na.pad=TRUE), group=Clinic)) + 
      geom_line(aes(color=Clinic)) + theme(axis.text.x=element_blank())
    
  
  
  # ROle change over time
  ggplot(data = cds, aes(x = ymd, y = Dist_from_reference, group=Role)) + geom_line(aes(color=Role))
  
  # get the wait times and visit duration for each clinic
  # convert to minutes. 
 a= visits  %>% 
    group_by(Clinic) %>%
    summarize(n = n(),
              NEvents = mean(NEvents, na.rm = TRUE),
              wait1 = mean(wait1*60, na.rm = TRUE),
              sdwait1 = sd(wait1*60, na.rm = TRUE),
              wait2 = mean(wait2*60, na.rm = TRUE),
              sdwait2 = sd(wait2*60, na.rm = TRUE),
              Duration  = mean(VisitDuration, na.rm = TRUE)
              )
  
 # For some reason, the sd function was not working correctly above...  do it again here. 
 visits  %>% 
    group_by(Clinic) %>%
    summarize(n = n(),
              sdwait1 = sd(wait1*60, na.rm = TRUE),
              sdwait2 = sd(wait2*60, na.rm = TRUE)
    )
  
  # get the precentage of actions by role... based on occurrences
 Clinic_role_pct = otr %>%
    group_by(Clinic, Role)  %>%
    summarize(n=n()) %>%
    mutate(RolePct = n/sum(n))  %>%
   select(Clinic, Role, RolePct) %>%
   spread(Role, RolePct)
    
 # get LOS by clinic based on visits
 Clinic_LOS_pct = visits %>%
   group_by(Clinic, LOS_CPT)  %>%
   summarize(n=n()) %>%
   mutate(LOSPct = n/sum(n)) %>%
   spread(Clinic, Role)
 
 # Make contingency table Clinic x  LOS
 C_LOS = visits %>%
   group_by(Clinic, LOS_CPT)  %>%
   summarize(n=n()) %>%
   spread( LOS_CPT,n)
 
 # Convert NA to zero and then convert to table. Carefully remove  rows/cols that make no sense
 C_LOS= as.data.frame(C_LOS)
 C_LOS[is.na(C_LOS)] =0
 C_LOS = matrix(unlist(C_LOS), nrow=6, ncol=20)[1:5,2:19]
 chisq.test(C_LOS)
 

# Make contingency table Clinic x  Role
C_R = otr %>%
  group_by(Clinic, Role)  %>%
  summarize(n=n()) %>%
  spread( Role, n)

# Convert NA to zero
C_R = as.data.frame(C_R)
C_R[is.na(C_R)] =0
C_R = matrix(unlist(C_R), nrow=5, ncol=10)[1:5,2:10]
chisq.test(C_R)

 Clinic_staffing = otr %>%
  group_by(Clinic,Role)  %>%
  summarize(total_staff=n_distinct(Role_ID)) %>%
  spread( Role, total_staff)
 
 
 # Get the handoffs per visit.  Use the distinct chunks in the Visit-Role threads (VRThrds)
 Role_Handoffs = VRThreads %>%
   group_by(Clinic,Visit_ID) %>%
   summarize(ho = n_distinct(threadNumVR))

 Role_Handoffs_by_clinic = Role_Handoffs %>%
   group_by(Clinic) %>%
   summarize(avg_ho = mean(ho), 
             median_ho = median(ho),
             max_ho = max(ho),
             sd_ho = sd(ho))
 
 # get the average staffing levels per day
 Clinic_daily_staffing = otr %>%
   group_by(Clinic,ymd,Role)  %>%
   summarize(total_staff=n_distinct(Role_ID)) 
 
 Clinic_AVG_daily_staffing = Clinic_daily_staffing %>%
   group_by(Clinic,Role)  %>%
   summarize(AVG_staff=mean(total_staff, na.rm = TRUE))  %>%
   spread( Role, AVG_staff)
 
 Clinic_AVG_daily_staffing = Clinic_daily_staffing %>%
   group_by(Clinic,Role)  %>%
   summarize(AVG_staff=median(total_staff, na.rm = TRUE))  %>%
   spread( Role, AVG_staff)
 
 # let's look at the actions by clinic_role, using the visit_role threads
actions_by_role =  VRThreads %>%
   group_by(Clinic,Role_VR) %>%
   summarize(avgActions = mean(Action_countVR))    %>%
   spread( Role_VR, avgActions)

# Look at alignment by clinic
visits %>%
  group_by(Clinic,Phase)  %>%
  summarize(w=max(ThreadDuration,na.rm = TRUE)  )   %>%
  spread( Phase, w)
 

# Get LOS by clinic-day
LOS_by_clinic_day = visits %>%
  group_by(Clinic_ymd)  %>%
  summarize(los=mean(LOS,na.rm = TRUE)  )


# See if diagnoses differ by  time period
Diag_by_phase = visits %>%
  group_by(Phase,Diagnosis_Group)  %>%
  summarize(n=n() )  %>%
  spread( Diagnosis_Group, n)

# See if actions differ by time period
#  cds_435_1 is a dataframe with one for for each clinic_day and ~300 columns for the actions
Action_by_phase = cds_435_1 %>%
  group_by(Phase)  %>%
  summarize_at(8:307,median ) 

}


remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

# converts visit_ID into character string with  factor  levels for actions
make_visit_string <- function(occ, vid){
  
  return(paste(as.character(as.numeric(occ$Action[occ$Visit_ID==as.character(vid)])),collapse = '-'))
}

# Converts list of actions into regex of factor levels  
make_process_pattern <- function(lvls,  p){
  
  return( paste( sapply(p, function(x) {as.character(which(lvls==x))}) , collapse = '.*') )
}

# Find process patterns
find_process_pattern <- function(occ,  p, vid){
  
  x = make_process_pattern( levels(ot$Action),  p)
  print(x)
  
   m = sapply(vid, function(v) {
                              (grepl( x, make_visit_string(occ, v))) } )
  # 
  

#  m = grep( x, make_visit_string(occ, vid))
  
  return(m)
  
}

# 
# make_box_plots <- function(){
# ggboxplot(ACHR_test[NEvents>100 & Clinic=='DRH'], x = "VisitMonth", y = "NetComplexity",
#           color = "VisitDay",
#           ylab = "Complexity", xlab = "Month (DRH)")
# }
# 
