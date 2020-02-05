dist1=get.distance(ac.coord=ac.coord,grid.coord=grid.coord,n.grid=n.grid,n.ac=n.ac)
dist2=GetDistance(AcCoord=data.matrix(ac.coord),GridCoord=data.matrix(grid.coord), 
                  Ngrid=n.grid, Nac=n.ac)
dist1[1:5,1:5]
dist2[1:5,1:5]

plot(dist1,dist2)
hist(dist1-dist2)
x2=(ac.coord$x[2]-grid.coord$x[1])^2
y2=(ac.coord$y[2]-grid.coord$y[1])^2
sqrt(x2+y2)

x2=(grid.coord$x[2]-ac.coord$x[1])^2
y2=(grid.coord$y[2]-ac.coord$y[1])^2
sqrt(x2+y2)