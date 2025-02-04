create.grid=function(dat,crs,extent,res,buffer) {
  grid<- raster(extent + buffer)
  res(grid)<- res
  proj4string(grid)<- crs
  grid[]<- 0
  
  grid
}
#------------------------------------------------
grid.summary.table=function(dat,crs,extent,res,buffer,method){  #dat must already have tseg and grid.cell assigned
  
  #method can be any one of "centroid", "mean", or "median", which represent the centroid of occupied grid cells or the mean/median of location of observations per occupied grid cell
  
  if (method == "centroid") {
    #create grid and extract coords per cell
    grid<- create.grid(dat=dat, crs=crs, extent=extent, res=res, buffer=buffer)
    grid.cell.locs<- coordinates(grid) %>% data.frame()
    names(grid.cell.locs)<- c("x", "y")
    grid.cell.locs$grid.cell<- 1:length(grid)
    grid.coord<- grid.cell.locs[grid.cell.locs$grid.cell %in% dat$grid.cell,]
    
  } else if(method == "mean") {
    grid.coord<- dat %>%
      group_by(grid.cell) %>%
      summarise(x = round(mean(x), 0), y = round(mean(y), 0)) %>%
      ungroup() %>%
      dplyr::select(x, y, grid.cell) %>%
      data.frame()
    
  } else if(method == "median") {
    grid.coord<- dat %>%
      group_by(grid.cell) %>%
      summarise(x = round(median(x), 0), y = round(median(y), 0)) %>%
      ungroup() %>%
      dplyr::select(x, y, grid.cell) %>%
      data.frame() 
  
  } else {
    stop("Method not recognized. Select either 'centroid', 'mean', or 'median'.")
  }
  
  
  
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
  nloc=length(unique(dat$grid.cell))
  
  
  #calculate # of obs in each grid.cell by tseg
  for (i in 1:length(dat.list)) {
    ntseg=max(dat.list[[i]]$tseg)
    res=matrix(0, ntseg, nloc)
    colnames(res)=1:nloc
  
    for (j in 1:ntseg){
      ind=dat.list[[i]] %>% filter(tseg == j) %>% group_by(grid.cell) %>% tally()
      res[j,ind$grid.cell]=ind$n #takes count of each cluster within given time segment
    }
    id1<- rep(unique(dat.list[[i]]$id), ntseg) %>% as.character()
    res=cbind(id = id1, res)
    obs.list[[i]]=res
  }
  obs<- do.call(rbind.data.frame, obs.list)
  obs<- mutate_at(obs, 2:ncol(obs), c("as.character")) %>% mutate_at(2:ncol(obs),
                                                                     c("as.numeric"))
  obs
}
#------------------------------------------------
plot.heatmap.ac=function(data, ac.coords, units) {
  #plots monthly, weekly, or daily proportions of AC use
  
  data<- data %>%
    mutate(day = lubridate::yday(date), week = lubridate::week(date),
           month = lubridate::month(date), year = lubridate::year(date))
  
  dat.sum<- data %>%
    group_by(year, !!units, ac) %>%
    tally() %>%
    group_by(year, !!units) %>%
    mutate(N=sum(n)) %>%
    mutate(prop = n/N)
  
  #create date var based on time unit of interest
  if (names(dat.sum)[2] == "month") {
    dat.sum$date<- dat.sum %>% 
      as_tibble() %>% 
      mutate(newyears = ymd(paste0(year,"-01-01"))) %>% 
      transmute(date = newyears + months(month-1)) %>% 
      pull()
    
  } else if (names(dat.sum)[2] == "week") {
    dat.sum$date<- dat.sum %>% 
      as_tibble() %>% 
      mutate(newyears = ymd(paste0(year,"-01-01"))) %>% 
      transmute(date = newyears + weeks(week-1)) %>% 
      pull()
    
  } else if (names(dat.sum)[2] == "day") {
    dat.sum$date<- dat.sum %>% 
      as_tibble() %>% 
      mutate(newyears = ymd(paste0(year,"-01-01"))) %>% 
      transmute(date = newyears + days(day-1)) %>% 
      pull()
  }
  
  ac.prop.df<- matrix(0, nrow(ac.coords)*length(unique(dat.sum$date)), 3)
  colnames(ac.prop.df)<- c("date","ac","prop")
  ac.prop.df[,1]<- rep(unique(dat.sum$date), each=nrow(ac.coords))
  ac.prop.df[,2]<- rep(1:nrow(ac.coords), length(unique(dat.sum$date)))
  for (i in 1:nrow(dat.sum)) {
    ind<- which(dat.sum$date[i] == ac.prop.df[,1] & dat.sum$ac[i] == ac.prop.df[,2])
    ac.prop.df[ind,3]<- dat.sum$prop[i]
  }
  ac.prop.df<- data.frame(ac.prop.df)
  
  
  print(
    ggplot(data = ac.prop.df, aes(x=lubridate::as_date(date), y=ac)) +
      geom_tile(aes(fill=prop), width = ceiling(365/max(dat.sum[,2]))) +
      scale_fill_viridis_c(paste("Proportion of\nobservations\nper", names(dat.sum)[2]),
                           limits = c(0,1)) +
      scale_y_continuous(trans = "reverse", breaks = 1:nrow(ac.coords), expand = c(0,0)) +
      scale_x_date(date_labels = "%b %Y", expand = c(0,0)) +
      labs(x="Time", y="Activity Center", title = paste("ID", unique(data$id))) +
      theme_bw() +
      theme(panel.grid = element_blank(), axis.title = element_text(size = 16),
            axis.text = element_text(size = 12), plot.title = element_text(size = 20))
  )
  
  
}
#------------------------------------------------
plot.ac.spatemp=function(data, ac.coords, units) {
  
  units<- enquo(units)
  
    par(ask = TRUE)
    map(data, ~plot.heatmap.ac(., ac.coords = ac.coords, units = units))
    par(ask = FALSE)
  
}
