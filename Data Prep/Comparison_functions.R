##########################################################################################################
# THREADNET:  Batch processing for larger data sets

# (c) 2017 Michigan State University. This software may be used according to the terms provided in the
# GNU General Public License (GPL-3.0) https://opensource.org/licenses/GPL-3.0?
# Absolutely no warranty!
##########################################################################################################


library(dplyr)




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
