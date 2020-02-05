create.grid=function(dat,crs,extent,res,buffer) {
  grid<- raster(extent + buffer)
  res(grid)<- res
  proj4string(grid)<- crs
  grid[]<- 0
  
  grid
}
#------------------------------------------------
grid.summary.table=function(dat,crs,extent,res,buffer){  #dat must already have tseg assigned
  
  #create grid and extract coords per cell
  grid<- create.grid(dat=dat, crs=crs, extent=extent, res=res, buffer=buffer)
  grid.cell.locs<- coordinates(grid) %>% data.frame()
  names(grid.cell.locs)<- c("x", "y")
  grid.cell.locs$grid.cell<- 1:length(grid)
  grid.coord<- grid.cell.locs[grid.cell.locs$grid.cell %in% dat$grid.cell,]
  
  
  grid.coord
}
#------------------------------------------------
df.to.list=function(dat) {  #only for id as col in dat
  id<- unique(dat$id)
  n=length(id)
  dat.list<- vector("list", n)
  names(dat.list)<- id
  
  for (i in 1:length(id)) {
    dat.list[[i]]<- dat[dat$id==id[i],]
  }
  dat.list
}
#------------------------------------------------
get.summary.stats_obs=function(dat){  #dat must have tseg assigned; for all IDs
  
  #change values of grid cells for easy manipulation
  dat$grid.cell<- as.factor(dat$grid.cell)
  levels(dat$grid.cell)<- 1:length(levels(dat$grid.cell))
  dat$grid.cell<- as.numeric(dat$grid.cell)
  
  #create list of input and to store output
  dat.list<- df.to.list(dat = dat)
  id<- unique(dat$id)
  n<- length(id)
  obs.list<- vector("list", n)
  names(obs.list)<- id
  
  
  #calculate # of obs in each grid.cell by tseg
  for (i in 1:length(dat.list)) {
    ntseg=max(dat.list[[i]]$tseg)
    nloc=length(unique(dat$grid.cell))
    res=matrix(0, ntseg, nloc)
    colnames(res)=1:nloc
  
    for (j in 1:ntseg){
      ind=dat.list[[i]] %>% filter(tseg==j) %>% group_by(grid.cell) %>% tally()
      res[j,ind$grid.cell]=ind$n #takes count of each cluster within given time segment
    }
    id<- rep(unique(dat.list[[i]]$id), ntseg)
    res=cbind(id, res)
    obs.list[[i]]=res
  }
  obs<- do.call(rbind.data.frame, obs.list) %>% data.matrix()
  obs
}
