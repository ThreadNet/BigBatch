##########################################################################################################
# THREADNET:  Batch processing for larger data sets

# (c) 2017 Michigan State University. This software may be used according to the terms provided in the
# GNU General Public License (GPL-3.0) https://opensource.org/licenses/GPL-3.0?
# Absolutely no warranty!
##########################################################################################################

# These are just examples -- need to  work out tables for particular purposes.

# These tables are based on cd=clinic_days and  v=visits or other tables with a similar structure
# The names of specific columns probably need to be adjusted

make_ACHR_tables <- function(cd,v) {
  
  library(dplyr)
  library(expss)
  
  # cd = apply_labels(cd,
  #                   Clinic_date =  'ClinicDate',
  #                   YMD_date = 'Date',
  #                   NEvents = 'Num Steps',
  #                   NumVisits = 'Num Visits',
  #                   NetComplexity= 'Complexity',
  #                   CompressRatio = 'Compressibility',
  #                   Clinic = 'Clinic',
  #                   NumUniqueProcedures = 'Num Unique Procedures',
  #                   NumUniqueDiagnosisGroups  = 'Num Unique Diagnosis Groups', 
  #                   NumPhysicians  = 'Num Physicians', 
  #                   Weekday  = 'Weekday', 
  #                   Month  = 'Month')
  # 
  expss_output_viewer()
  
  cd %>% 
    tab_cells(A_NetComplexity,NumVisits,NEvents/NumVisits,NumUniqueProcedures,NumUniqueDiagnosisGroups,
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
    tab_cells( NEvents, VisitDuration, A_NetComplexity, CompressRatio, Action_count, Workstation_count, Role_count ) %>% 
    tab_cols(total(), Clinic) %>% 
    tab_stat_mean() %>% 
    tab_pivot()
  
  v %>% 
    tab_cells( NEvents, VisitDuration, A_NetComplexity,  Action_count, Workstation_count, Role_count ) %>% 
    tab_cols(total(), Diagnosis_group) %>% 
    tab_stat_mean() %>% 
    tab_pivot() %>% 
    tab_transpose()
  
  v %>% 
    tab_cells(  NEvents, VisitDuration, A_NetComplexity,  Action_count ) %>% 
    tab_cols(total(), Diagnosis_group) %>% 
    tab_stat_mean_sd_n() %>% 
    tab_pivot() %>% 
    tab_transpose()
  
  v %>% 
    tab_cells( NEvents, VisitDuration, A_NetComplexity,   Action_count  ) %>% 
    tab_cols(total(), Physician) %>% 
    tab_stat_mean_sd_n() %>% 
    tab_pivot() %>% 
    tab_transpose()
  
   
  
}


