
# these functions will take clinic data and produce trajectories (distance of the graph from a reference graph)

# start with auditfinal_10022018.csv --> use ACHR_Batch_v2 to make new_occ.  
# The code here starts with new_occ (~ 7.72m occurrences that we break down into 1857 clinic_days and 57835 visits)

library(tidyr)
library(data.table)
library(dplyr)
library(ThreadNet)
library(ngram)
library(lubridate)
library(expss)
library(philentropy)

# Use this to put the trajectory data with the clinic_day data  
merge_clinic_days_and_trajectories <- function(){
  
  clinic_days$Clinic_date = gsub('HH_POB','HHPOB',clinic_days$Clinic_date)
  cds=merge(x=test0,y=clinic_days,by.x='Clinic_ymd',by.y ='Clinic_date' ,all=TRUE) 
  
}

# these functions support the moving window

# Simplified...  you just pick the unique TN
get_threadList <- function(e,TN){ return(unique(e[[TN]])) }

# always just one clinic day -- so use the "bucket" definition from ACHR batch... 
# e = new_occ
# b = bucket
get_bucket <- function(o, b ){
  
  # get the list of threads
 return( filter(o,Clinic_ymd==b) )
  
}



#' @param e data frame for POV
#' @param blist is the bucket list.  Each bucket is a "window"
#' @param cf is the column that defines events (e.g. 'Action_Role')
#' @export
graph_trajectory  <- function(e, cf){
  
  # get total size of possible matrix
  Max_Order = length(unique(e[[cf]]))^2
  
  
  # make data frame for results
  vt=data.frame( ngrams=character(), freq=integer(), bid=integer() )
  
  # here is the list of buckets
  blist = unique(e[['Clinic_ymd']])
  
  # take the underscore  out so we can use it to strsplit below... and sort it
  blist=gsub('HH_POB','HHPOB',blist)
  blist = sort(blist)
  
  # now many buckets?
  nb = length(blist)
  print(blist)
  print(paste('nb =',nb))
  
  complexity_idx=list(nb)
  
  bcount=0
  # scan through the data
  for (b in blist){
    bcount= bcount +1
   # print(paste('b =',b))
    
    # get text vector for the whole data set.  Bucket needs at least 2 threads
      th = get_bucket(e, b)
      if (nrow(th)>2) { ngdf = count_ngrams(th, 'threadNum', cf, 2)[1:2]  }
      
      nodes = length(unique(th[[cf]]))
      edges = nrow(ngdf)
      complexity_idx[bcount] = estimate_task_complexity_index( nodes ,edges )
     
      
     # print(paste('nrow ngdf =',nrow(ngdf)))
    
    # add the bucket number and name
    ngdf$bid = bcount
    # ngdf$bname = b
    
    # append the columns to the end
    # vt is the whole set of all ngrams in all the windows
    vt=rbind(vt,ngdf)
  }
   print(head(complexity_idx))

   
  # convert to factor
  vt$ngrams = factor(vt$ngrams)
  
  # nWindows = length(unique(vt$bid))
  nWindows = nb
  #print(paste('nWindows =',nWindows))
 
  # get the set of unique ngrams for the whole data set
  vt_unique = data.frame(ngrams=unique(vt$ngrams))
  #print(vt_unique)
  
  # put the results here
  windowFreqMatrix = matrix(0,nrow=nWindows, ncol=nrow(vt_unique))
  
  for (i in 1:nWindows){
    
    # get the merged list
    vtmerge = merge(x=vt_unique, y=vt[vt$bid==i,], by='ngrams', all.x = TRUE)
    
    # use the bid.y to get the whole vector, but replace the NA with zeros
    bb=vtmerge[vtmerge$bid==i,'freq']
    bb[is.na(bb)] <- 0
    
    windowFreqMatrix[i,]=bb
    
    print( paste('num non-zero=', sum(windowFreqMatrix[i,] > 0)))
    
  }
  
 # return(windowFreqMatrix)
  
  #  correlate each row with the first one stick it in a dataframe
  df =data.frame(window=2:nWindows,
                 Clinic_ymd = blist[2:nWindows],
                 Clinic= unlist(lapply(2:nWindows,
                                       function(i){ unlist(strsplit(blist[i],'_'))[1] })),
                 ymd=unlist(lapply(2:nWindows,
                                   function(i){ unlist(strsplit(blist[i],'_'))[2] })),
                 pct_retained = unlist(lapply(2:nWindows,
                                              function(i){sum((windowFreqMatrix[i-1,]>0)+(windowFreqMatrix[i,]>0)==2)/
                                                  sum(windowFreqMatrix[i-1,]>0)
                                                  })) ,
                 pct_possible = unlist(lapply(2:nWindows,
                                              function(i){sum(windowFreqMatrix[i,] > 0)/Max_Order })) ,
                 complexity = unlist( complexity_idx[2:nWindows]) ,
                 corr_with_first= unlist(lapply(2:nWindows,
                                            function(i){abs( cor(windowFreqMatrix[i,],windowFreqMatrix[1,])) })),
                 cosine_dist_from_first= unlist(lapply(2:nWindows,
                                                function(i){distance(rbind(windowFreqMatrix[i,],windowFreqMatrix[1,]),
                                                                     method='cosine' ) })),
                 corr_with_next= unlist(lapply(2:nWindows,
                                                function(i){abs( cor(windowFreqMatrix[i-1,],windowFreqMatrix[i,])) })),
                 cosine_dist_from_next= unlist(lapply(2:nWindows,
                                                function(i){distance(rbind(windowFreqMatrix[i-1,],windowFreqMatrix[i,]),
                                                                     method='cosine' ) }))
                 )
  

return(df)
  
  # # get the ngram data and labels
  # b_df=as.data.frame(windowFreqMatrix[2:nWindows,])
  # colnames(b_df)=vt_unique$ngrams
  # 
  # # stick the ngram frequencies on the end for good measure
  # return(cbind(df,b_df))

}


plots_for_papers <- function(){
  
  # compare the complexity of the five clinics
  plot(cds$Clinic.x,cds$A_NetComplexity)
  
  
  
}
 