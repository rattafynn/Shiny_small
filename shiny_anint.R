# anint function for periodic boundary conditions (PBC)
# It ensures that an atom never leaves the simulation box
# x is a real number
# Raymond Atta-Fynn
# NYC Data Academy, August 2, 2019
anint <- function(x){
  if(x>0){
    return(as.integer(x+0.5))
  } else{
    return(as.integer(x-0.5))
  }
}
