#Raymond Atta-Fynn
#NYC Data Academy
#Routine for building a model for water molecules for atomistic modeling
#August 2, 2019

library(rgl)

#**************************************************************************************
#               RGL VISUALIZATION ROUTINES                                            *
#The routines from lines 16-85 were adopted from the following website:               *
#        http://www.sthda.com/english/wiki/print.php?id=211                           *
#                                                                                     *
# I added the routine from lines 87-102 to build the bonds between atoms              *
#**************************************************************************************

#The function rgl_init() will create a new RGL device if requested or if there is no opened device:
#' @param new.device a logical value. If TRUE, creates a new device
#' @param bg the background color of the device
#' @param width the width of the device
#' This routine was adopted from: http://www.sthda.com/english/wiki/print.php?id=211
rgl_init <- function(new.device = FALSE, bg = "white", width = 800) { 
  if( new.device | rgl.cur() == 0 ) {
    rgl.open()
    par3d(windowRect = 50 + c( 0, 0, width, width ) )
    rgl.bg(color = bg )
  }
  rgl.clear(type = c("shapes", "bboxdeco"))
  rgl.viewpoint(theta = 15, phi = 20, zoom =0.7)
}

get_colors <- function(groups, group.col = palette()){
  groups <- as.factor(groups)
  ngrps <- length(levels(groups))
  if(ngrps > length(group.col)) 
    group.col <- rep(group.col, ngrps)
  color <- group.col[as.numeric(groups)]
  names(color) <- as.vector(groups)
  return(color)
}

#rgl_add_axes(): A custom function to add x, y and z axes
# x, y, z : numeric vectors corresponding to
#  the coordinates of points
# axis.col : axis colors
# xlab, ylab, zlab: axis labels
# show.plane : add axis planes
# show.bbox : add the bounding box decoration
# bbox.col: the bounding box colors. The first color is the
# the background color; the second color is the color of tick marks
#' This routine was adopted from: http://www.sthda.com/english/wiki/print.php?id=211
rgl_add_axes <- function(x, y, z, L,axis.col = "grey",
                         xlab = "X", ylab="Y", zlab="Z", show.plane = FALSE, 
                         show.bbox = FALSE, bbox.col = c("#333377","black"),lwd=6,size=10)
{ 
  
  lim <- function(x){c(-x/2,x/2)}
  # Add axes
  xlim <- lim(L); ylim <- lim(L); zlim <- lim(L)
  rgl.lines(xlim, c(0, 0), c(0, 0), color = "black",lwd=4)
  rgl.lines(c(0, 0), ylim, c(0, 0), color = "black",lwd=4)
  rgl.lines(c(0, 0), c(0, 0), zlim, color = "black",lwd=4)
  
  # Add a point at the end of each axes to specify the direction
  axes <- rbind(c(xlim[2], 0, 0), c(0, ylim[2], 0), 
                c(0, 0, zlim[2]))
  rgl.points(axes, color = "black", size = 3)
  #  shapelist3d(axes)
  # Add axis labels
  rgl.texts(axes, text = c(xlab, ylab, zlab), color = "black",
            adj = c(0.5, -0.8), size = 15)
  #
  # Add plane
  #  if(show.plane) 
  #    xlim <- xlim/1.1; zlim <- zlim /1.1
  #  rgl.quads( x = rep(xlim, each = 2), y = c(0, 0, 0, 0),
  #             z = c(zlim[1], zlim[2], zlim[2], zlim[1]))
  
  # Add bounding box decoration
  if(show.bbox){
    rgl.bbox(color=c(bbox.col[1],bbox.col[1]), alpha = 0.5,
             emission=bbox.col[1], specular=bbox.col[1], shininess=5, 
             xlen = 0, ylen = 0, zlen = 0,draw_front = FALSE,lwd=6) 
  }
  
}

# Build bond between atoms based on cut-off distances
rgl_bond <- function(x,y,z,nr,rcut,bond_color,bond_thickness){
  if(nr%%3==0){
    for(i in 1:nr-1){
      c1 <- c(x[i],y[i],z[i])
      for(j in (i+1):nr){
        c2 <- c(x[j],y[j],z[j])
        r <- sqrt(sum((c2-c1)^2))
        if (r<=rcut){
          lines3d(c(x[i], x[j]), c(y[i], y[j]), c(z[i], z[j]),
                  lwd=bond_thickness,color =bond_color)
        }
      }
    }
  } else {
    for(i in 1:(nr-1)){
      c1 <- c(x[i],y[i],z[i])
      for(j in (i+1):nr){
        c2 <- c(x[j],y[j],z[j])
        r <- sqrt(sum((c2-c1)^2))
        if (r<=rcut){
          lines3d(c(x[i], x[j]), c(y[i], y[j]), c(z[i], z[j]),
                  lwd=bond_thickness,color =bond_color)
        }
      }
    }
    c2 <- c(x[nr],y[nr],z[nr])
    for(i in 1:(nr-1)){
      c1 <- c(x[i],y[i],z[i])
      r <- sqrt(sum((c2-c1)^2))
      if (r<=3.0 & ((i-1)%%3==0)){
        lines3d(c(x[i], x[nr]), c(y[i], y[nr]), c(z[i], z[nr]),
                lwd=bond_thickness,color ="green")
      } 
    }
  }
  
}