# Define function to create a multistate m-array from capture history data
marray <- function(ch, unobs = 0){
  ns <- length(table(ch)) - 1 + unobs
  no <- ncol(ch)
  out <- matrix(0, ncol = ns*(no-1)+1, nrow = ns*(no-1))
  
  # Remove capture histories of individuals that are marked at last occasion
  get.first <- function(x) min(which(x!=0))
  first <- apply(ch, 1, get.first)
  last.only <- which(first==no)
  if (length(last.only) > 0) ch <- ch[-last.only,]
  
  # Compute m-array
  for (i in 1:nrow(ch)){
    cap.occ <- which(ch[i,]!=0)
    state <- ch[i,cap.occ]
    if (length(state) == 1) {
      out[state[1]+ns*(cap.occ[1]-1), ns*(no-1)+1] <- out[state[1]+
                                                      ns*(cap.occ[1]-1), ns*(no-1)+1] + 1
    }
    
    if (length(state) > 1) {
      for (t in 2:length(cap.occ)){
        out[(cap.occ[t-1]-1)*ns+state[t-1], (cap.occ[t]-2)*
              ns+state[t]] <- out[(cap.occ[t-1]-1)*ns+state[t-1], (cap.occ[t]-2)*
                                    ns+state[t]] + 1
      } # t
      
      if (max(cap.occ) < no){
        out[(cap.occ[t]-1)*ns+state[t], ns*(no-1)+1] <- out[(cap.occ[t]-1)*
                                                              ns+state[t], ns*(no-1)+1] + 1
      } # i
      
    } # if
    
  } # i
  
  return(out)
  
}    

# Excecute function to create m-array from capture history data

marr <- marray(CH)
