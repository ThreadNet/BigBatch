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
  ggplot(data = cds, aes(x = ymd, y = corr_with_first, group=Clinic)) + geom_line(aes(color=Clinic))
  ggplot(data = cds, aes(x = ymd, y = rollmean(corr_with_first,5,na.pad=TRUE), group=Clinic)) + geom_line(aes(color=Clinic))
  
}

# 
# make_box_plots <- function(){
# ggboxplot(ACHR_test[NEvents>100 & Clinic=='DRH'], x = "VisitMonth", y = "NetComplexity",
#           color = "VisitDay",
#           ylab = "Complexity", xlab = "Month (DRH)")
# }
# 
