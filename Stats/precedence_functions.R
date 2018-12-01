
# need to test if this is working correctly
# It should be counting the number of times that each element in the lexicon precedes each other element.

# make some data
S=lapply(1:10,function(x) {c('a','b','a','c',sample(letters[5:6],1+round(runif(1)*10),replace=TRUE)) })

# take list of sequences, S
precedes <- function(S){
  
  # Get lexicon
  lex = sort(unique(unlist(S)))
  
  # for each element in the lexicon...
  result = matrix(data = sapply(lex, function(a){
    
      # Find its position(s) in each sequence. pos is a list of vectors
      pos=sapply(S, function(s) {grep(a,s)})
    
      # now for the last position in each string, find the lexical elements that precede it 
      # get the substrings
      subS = mapply(function(p,s){if (max(p)>1) s[1:max(p)]},pos, S)
 
      # Now find all of the other lexical elements that precede a in each subS
      colSums(sapply(lex,function(x){str_detect(subS,x)})) }
      ) , 
      
      ncol=length(lex))

  rownames(result) = lex
  colnames(result) = lex
  
  # normalize by the number of sequences,so "always" = 1
  result = result/length(S)
  
  return(result)
  
}