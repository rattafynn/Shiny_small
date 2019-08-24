
# Generates data for RDF plot
# Raymond Atta-Fynn
# NYC Data Academy, August 2, 2019
rdfplot <- function(pos,L,periodic){
  rmin<-1.5
  rmax<-L/2.0
  nbin<-51
  dr<-(rmax-rmin)/(nbin-1)
  N<-nrow(pos)
  myhist<-rep(0,nbin)
  
  
  for (i in 1:(N-1)){
    if(pos[i,1]=="O"){
      tmp1<-pos[i,2:4]
      for (j in (i+1):N){
        if(pos[j,1]=="O"){
          tmp2<-pos[j,2:4]
          dist <- min_vec_dist(tmp1,tmp2,L,periodic)[[2]]
          bin=as.integer((dist-rmin)/dr) + 1
          if(bin<nbin){
            myhist[bin] <- myhist[bin] + 1 
          }
        }
      }
    }
  }
  
  rdfdata <- data.frame(r0=numeric(),rdfoutput=numeric())
  rho <- N/L^3
  for (i in 1:nbin){
    rl <- rmin + (i-1)*dr
    ru <- rl + dr
    norm_factor <- (4*pi*rho/3)*(ru^3 - rl^3)
    myhist[i] <- myhist[i]/norm_factor
    rdfdata <- rbind(rdfdata,data.frame(r0=0.5*(ru+rl),rdfoutput=myhist[i]))
  }
  return(rdfdata)
}
