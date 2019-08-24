# This function returns:
# the separation vector and distance between a pair of atoms
# The logical parameter periodic determines whether or not PBC is enforced 
# Raymond Atta-Fynn
# NYC Data Academy, August 2, 2019
min_vec_dist <- function(x,y,L,periodic){
  z <- y-x
  if(periodic==TRUE){
    z[1] <- z[1]-L*anint(z[1]/L)
    z[2] <- z[2]-L*anint(z[2]/L)
    z[3] <- z[3]-L*anint(z[3]/L)
  }
  return(list(z,sqrt(sum(z^2))))
}
