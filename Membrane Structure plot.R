# Load necessary packages:

library(rgl)

# Load membrane data:

setwd("C:/Users/ferry/OneDrive/Dokumente/__UNI/9. Semester/Stallfudien 2/Projekt 2")
load("Finite_Element_Mesh.RData")
summary(membrane)
summary(edge)
summary(support)
summary(truss)

######################
##  Zeichne Membran ##
######################
open3d()
axes3d()
for(i in 1:nrow(membrane)){
  triangle <- coordinates[membrane[i, 2:4], 2:4]
  triangles3d(triangle, alpha = 0.5)
  lines3d(rbind(triangle, triangle[1,])) # to distinguish each segment
}
# draw support
for(i in 1:nrow(support)){
  lines3d(coordinates[support[i, 2:3], 2:4], color = "black")}
# draw truss
for(i in 1:nrow(truss)){
  lines3d(coordinates[truss[i, 2:3], 2:4], lwd = 2)}

rgl.snapshot("sunsail.png")


