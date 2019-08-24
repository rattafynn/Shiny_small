#************************************************************************
#               MOLECULAR MODEL BUILDING ROUTINES                       *
#                     Raymond Atta-Fynn                                 *
#                      NYC Data Academy                                 *
#           Routine generating initial configuration                    *
#                       August 2, 2019                                  *
#************************************************************************


# This function checks whether a pair of atoms are closer than within rcut
# iflag =0 is returned if atoms are separated by a distance greater than rcut
# iflag =1 is returned if atoms are separated by a distance less than or equal to rcut
# Raymond Atta-Fynn
# NYC Data Academy, August 2, 2019
closest_approach <- function(j,pos,y,L,rcut,periodic){
  iflag <- 0
  for (i in 1:j){
    atom <-pos[i,1]
    if(atom=="O"){
          z <- pos[i,2:4]
         dist <- min_vec_dist(z,y,L,periodic)[[2]]
         if (dist < rcut){
             iflag <- 1
             break  
         }

      }
    }
  return(iflag)
}

closest_approach_ion <- function(z,y,L,rcut,periodic){
  iflag <- 0
  dist <- min_vec_dist(z,y,L,periodic)[[2]]
  if(dist < rcut)iflag <- 1
  return(iflag)
}

# This function generates the random x,y,z positions of a given atom in a cubic box of length L 
# -L/2 <= x <= L/2
# -L/2 <= y <= L/2
# -L/2 <= z <= L/2
# Raymond Atta-Fynn
# NYC Data Academy, August 2, 2019
gen_random_position <- function(L){
  return(runif(3,-L/2.0,L/2.0))
}

# This function generates the random vector along which an OH bond will be aligned
# Raymond Atta-Fynn
# NYC Data Academy, August 2, 2019
gen_H_vector <- function(x){
  vmod =0
  while(vmod<1E-4){
    tmp <- runif(3,-x/2,x/2)
    vmod = sqrt(sum(tmp^2))
  }
  return(tmp/vmod)
}

# This function computes the angle formed by the water bonds
# Periodic boundary may or may not be enforced
# Raymond Atta-Fynn
# NYC Data Academy, August 2, 2019
compute_angle <-function(x,y,z,L,periodic){
  tmp1 <- min_vec_dist(x,y,L,periodic)
  tmp2 <- min_vec_dist(x,z,L,periodic)
  v1=tmp1[[1]]
  v2=tmp2[[1]]
  
  d1=tmp1[[2]]
  d2=tmp2[[2]]
  
  return(acos(sum(v1*v2)/(d1*d2))*180/pi)
}

# This function attaches two hydrogen atoms attached to oxygen to form water 
# -L/2 <= x <= L/2, -L/2 <= y <= L/2, -L/2 <= z <= L/2
# OH_dist is the O-H bond distance
# Raymond Atta-Fynn
# NYC Data Academy, August 2, 2019
gen_H_atoms <- function(x,L,OH_dist,periodic){
  y <- x + OH_dist*gen_H_vector(2)
  angle <- 104.5
  while (abs(angle) > 3.0){
    z <- x + OH_dist*gen_H_vector(2)
    angle <- compute_angle(x,y,z,L,periodic) - 104.5
  } 
  return(list(y,z))
}


# This are the main function which actually generates the model
# It takes only three inputs: 
# the total of water molecules N, the system type, and periodicity
# Raymond Atta-Fynn
# NYC Data Academy, August 2, 2019
generate_initial_configuration <- function(N,system_type,periodic,metal,model){
  avogadro <- 6.02214E23
  O_mass <- atomic_mass("O")
  H_mass <- atomic_mass("H")
  water_mass <- O_mass + 2*H_mass
  total_mass <- N*water_mass
  density <- 1.0
  rcut <- 2.2 
  OH_dist<-1.0
  L <- 1.0E8*(total_mass/(density*avogadro))^(1.0/3.0)  

#  if(system_type=="pure water"){
#    L = 1.0E8*(total_mass/(density*avogadro))^(1.0/3.0)  
#  } else if(system_type=="Cm(III) in Water"){
#    total_mass = total_mass + 247
#    L = 1.0E8*(total_mass/(density*avogadro))^(1.0/3.0)
#  }  else if(system_type=="Al(III) in Water"){
#    total_mass = total_mass + 27
#    L = 1.0E8*(total_mass/(density*avogadro))^(1.0/3.0)
#  }  else if(system_type=="Th(IV) in Water"){
#    total_mass = total_mass + 232
#    L = 1.0E8*(total_mass/(density*avogadro))^(1.0/3.0)
#  }

# Initialize atom counter and empty frame  
  i <- 0
  pos <- data.frame(ATOM=character(),X=numeric(),Y=numeric(),Z=numeric())

# This is the initial position of the ion
# It would be assigned last
  tmp0 <-rep(0.0,3) 
  
# Now begin to build the model  
  if(system_type=="pure water"){
    tmp = gen_random_position(L)
    while( i < 3*N){
      tmp = gen_random_position(L)
      if(i==0){
        iflag  <- 0
      } else {
        iflag  <- closest_approach(i,pos,tmp,L,rcut,periodic)
      }

      if(iflag==0){
        pos0 <- gen_H_atoms(tmp,L,OH_dist,periodic)
        H1   <- pos0[[1]]
        H2   <- pos0[[2]]
        pos <- rbind(pos,data.frame(ATOM='O',X=tmp[1],Y=tmp[2],Z=tmp[3]))
        pos <- rbind(pos,data.frame(ATOM='H',X=H1[1],Y=H1[2],Z=H1[3]))
        pos <- rbind(pos,data.frame(ATOM='H',X=H2[1],Y=H2[2],Z=H2[3]))
        i <- i + 3
        print(paste("Generated water molecule",i/3))
      }
    }
  } else {
    tmp = gen_random_position(L)
    while( i < 3*N){
      tmp = gen_random_position(L)
      if(i==0){
        iflag  <- 0
      } else {
        iflag  <- closest_approach(i,pos,tmp,L,rcut,periodic)
      }
      iflag0 <- closest_approach_ion(tmp0,tmp,L,rcut,periodic)
  
      if(iflag0==0 & iflag==0){
        pos0 <- gen_H_atoms(tmp,L,OH_dist,periodic)
        H1   <- pos0[[1]]
        H2   <- pos0[[2]]
        pos <- rbind(pos,data.frame(ATOM='O',X=tmp[1],Y=tmp[2],Z=tmp[3]))
        pos <- rbind(pos,data.frame(ATOM='H',X=H1[1],Y=H1[2],Z=H1[3]))
        pos <- rbind(pos,data.frame(ATOM='H',X=H2[1],Y=H2[2],Z=H2[3]))
        i <- i + 3
        print(paste("Generated water molecule",i/3))
      }
    }
    pos <- rbind(pos,data.frame(ATOM=metal,X=tmp0[1],Y=tmp0[2],Z=tmp0[3]))
    print(paste("Generated position of metal atom",metal))
  }

# Atomic masses  
  mass<-c()
  q<-c()
# conv_factor converts the atomic mass from a.m.u. to (eV*fs^2)/A^2
  conv_factor=103.6426915168282606
  for(i in 1:nrow(pos)){
    mass[i]<-atomic_mass(pos[i,1])*conv_factor
    q[i]<-ion_water_12_6_4_parameters(pos[i,1],model)[[4]]
  }

  energy_param<-c()
  if(system_type=="pure water"){
    tmp <-ion_water_12_6_4_parameters("O",model)
    energy_param[1]<-tmp[[1]]
    energy_param[2]<-tmp[[2]]
  } else {
    tmp <-ion_water_12_6_4_parameters("O",model)
    energy_param[1]<-tmp[[1]]
    energy_param[2]<-tmp[[2]]
    tmp <-ion_water_12_6_4_parameters(metal,model)
    energy_param[3]<-tmp[[1]]
    energy_param[4]<-tmp[[2]]
    energy_param[5]<-tmp[[3]]
  }

  return(list(L,pos,mass,q,energy_param))
}
