##########################################################################################################
# THREADNET:  Batch processing for larger data sets

# (c) 2017 Michigan State University. This software may be used according to the terms provided in the
# GNU General Public License (GPL-3.0) https://opensource.org/licenses/GPL-3.0?
# Absolutely no warranty!
##########################################################################################################

library(parallel)

ptest <-function(o) {
  x<<-0
  oo<<-o
  
  l=1:nrow(o)
   # return( sapply(l,f1 ) )

  new_o = cbind(o,  sapply(l,f1 ) )
   return(new_o)
}

f1 <- function(r){
  x<<-x+1
  return(as.character(oo[r,5]))
}

