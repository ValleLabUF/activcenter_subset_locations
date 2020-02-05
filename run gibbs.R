rm(list=ls(all=TRUE))
set.seed(10)

#read important functions
library('Rcpp')
sourceCpp('aux1.cpp')
source('gibbs sampler.R')
source('gibbs functions.R')

#get data
dat=read.csv('fake data.csv',as.is=T) #frequency of visitation in each location (column) for each time segment (row)
grid.coord=read.csv('fake data grid.csv',as.is=T) #geographical coordinates of locations

#basic setup
ngibbs=1000
nburn=ngibbs/2
n.ac=20
gamma1=0.1

#initial coordinates of activity centers (define this based on data instead of uninformative as below)
ind=sample(nrow(grid.coord),size=n.ac)
ac.coord.init=grid.coord[ind,]

#potential locations for activity centers (AC)
possib.ac=grid.coord #these don't have to be identical (i.e., we can define AC's on a coarser grid)

#run gibbs sampler
options(warn=2)
dat=data.matrix(dat)
res=gibbs.activity.center(dat=dat,grid.coord=grid.coord,n.ac=n.ac,
                          ac.coord.init=ac.coord.init,gamma1=gamma1,
                          possib.ac=possib.ac)

