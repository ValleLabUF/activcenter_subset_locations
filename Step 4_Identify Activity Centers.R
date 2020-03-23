set.seed(10)

#load libraries and read important functions
library('Rcpp')
library(raster)
library(sf)
library(tidyverse)
library(lubridate)
library(viridis)
library(progress)
library(rnaturalearth)
library(rnaturalearthdata)


sourceCpp('aux1.cpp')
source('gibbs sampler.R')
source('gibbs functions.R')
source('helper functions.R')



#################
### Load data ###
#################


#load data
dat<- read.csv("Snail Kite Gridded Data_TOHO15km.csv", header = T, sep = ",")
obs<- get.summary.stats_obs(dat)  #frequency of visitation in each location (column) for each time segment (row)
obs1<- as.matrix(obs[,-1])  #for proper use by model

#geographical coordinates of locations
utm.crs<- CRS('+init=epsg:32617')
extent<- extent(min(dat$x), max(dat$x), min(dat$y), max(dat$y))
res<- 15000
buffer<- 2*res
grid.coord<- grid.summary.table(dat=dat, crs=utm.crs, extent=extent, res=res, buffer=buffer,
                                method = "median")

#Load world map data
usa<- ne_states(country = "United States of America", returnclass = "sf")
fl<- usa %>% filter(name == "Florida") %>% sf::st_transform(fl, crs = "+init=epsg:32617")

# lakes
lakes10 <- ne_download(scale = 10, type = 'lakes', category = 'physical', returnclass = "sf")
lakes10<- sf::st_transform(lakes10, crs = "+init=epsg:32617") %>%
  sf::st_crop(xmin = min(dat$x-20000), xmax = max(dat$x+20000), ymin = min(dat$y-20000),
              ymax = max(dat$y+20000))


#Define initial activity centers (top 50 by # of obs)
tmp<- colSums(obs[,-1]) %>% data.frame(grid.cell = colnames(obs[,-1]), nobs = .) %>%
      arrange(desc(nobs)) %>% slice(n=1:50) %>% dplyr::select(grid.cell)
tmp<- tmp$grid.cell %>% as.character() %>% as.numeric()
ind<- sample(tmp, size = 20, replace = F)
ac.coord.init<- grid.coord[ind,]

#top 50
ac.coord.init2<- grid.coord[tmp,]

#potential locations for activity centers (AC)
possib.ac=grid.coord #these don't have to be identical (i.e., we can define AC's on a coarser grid)


#Viz possible AC locations
ggplot() +
  geom_sf(data = fl) +
  geom_sf(data = lakes10, fill = "lightblue", alpha = 0.65) +
  coord_sf(xlim = c(min(dat$x-10000), max(dat$x+10000)),
           ylim = c(min(dat$y-10000), max(dat$y+10000)), expand = FALSE) +
  geom_point(data = dat, aes(x, y, fill = "Observations"), shape = 21, size = 1,
             alpha = 0.2) +
  geom_point(data = grid.coord, aes(x, y, fill = "Possible AC Locations"), shape = 21,
             size = 2) +
  labs(x="Longitude", y="Latitude") +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  scale_fill_manual("", values = c("grey55","red")) +
  theme_bw() +
  theme(legend.position = "top", axis.title = element_text(size = 16),
        axis.text = element_text(size = 12), legend.text = element_text(size = 12))


#########################
### Run Gibbs sampler ###
#########################

#basic setup
ngibbs=5000
nburn=ngibbs/2
n.ac=20
gamma1=0.1

#run gibbs sampler
options(warn=2)

dat.res=gibbs.activity.center(dat=obs1,grid.coord=grid.coord[,-3],n.ac=n.ac,
                          ac.coord.init=ac.coord.init[,-3],gamma1=gamma1,
                          possib.ac=possib.ac[,-3])


#plot output and look at frequency of AC visitation
plot(dat.res$logl,type='l')
plot(dat.res$phi,type='l')

##############################################
### Extract AC Coordinates and Assignments ###
##############################################

##use ACs from iteration with max log likelihood (after burn-in)
ML<- dat.res$logl %>% order(decreasing = T)
ML<- ML[ML > 500][1]
ac<- dat.res$z[ML,]
phi<- dat.res$phi[ML,]
ac.coords<- matrix(NA, length(unique(ac)), 2)
colnames(ac.coords)<- c("x","y")
tmp<- dat.res$coord[ML,]

for (i in 1:length(unique(ac))) {
  ac.coords[i,]<- round(c(tmp[i], tmp[i+length(unique(ac))]), 0)
}

#reorder ACs from N -> S
ac.coords<- data.frame(ac.coords, ac=1:length(unique(ac))) %>% .[order(.$y, decreasing = TRUE),]
ac.coords<- cbind(ac.coords, ac.ns = 1:nrow(ac.coords))

ac.ns<- list()
for (i in 1:length(ac)) {
  ac.ns[i]<- which(ac[i] == ac.coords$ac)
}
ac.ns<- unlist(ac.ns)  #ac assignments order N -> S

table(ac)
table(ac.ns)





############################
### Add ACs to Dataframe ###
############################

tseg.length<- dat %>% group_by(id, tseg) %>% tally()
tseg.length<- tseg.length$n
ac.aug<- rep(ac.ns, times = tseg.length)

dat$ac<- ac.aug

#Calculate number of obs per AC
dat %>% group_by(id) %>% dplyr::select(ac) %>% table()



## Map

# ac.coords<- read.csv("Activity Center Coordinates_TOHO.csv", header = T, sep = ',')
# dat<- read.csv("Snail Kite Gridded Data_AC_TOHO.csv", header = T, sep = ',')



nests<- dat %>% group_by(id) %>% dplyr::select(c(id, x, y)) %>% slice(n=1)


# ACs and initial values
ggplot() +
  geom_sf(data = fl) +
  geom_sf(data = lakes10, fill = "lightblue", alpha = 0.65) +
  coord_sf(xlim = c(min(dat$x-20000), max(dat$x+20000)),
           ylim = c(min(dat$y-20000), max(dat$y+20000)), expand = FALSE) +
  geom_point(data = dat, aes(x, y, fill = "Raw"), shape = 21, size = 1,
             alpha = 0.2) +
  geom_point(data = nests, aes(x, y, fill = "Nests"), shape = 24, size = 3.5, alpha = 0.7) +
  geom_point(data = ac.coords, aes(x, y, fill = "ACs"), shape = 21, size = 3, alpha = 0.8) +
  labs(x="Longitude", y="Latitude") +
  scale_fill_manual("", values = c(viridis(n=5)[4],"red","grey55")) +
  guides(fill = guide_legend(override.aes = list(shape = c(21,24,21)))) +
  theme_bw()

# ACs and snail kite locs
ggplot() +
  geom_sf(data = fl) +
  geom_sf(data = lakes10, fill = "lightblue", alpha = 0.65) +
  coord_sf(xlim = c(min(dat$x-20000), max(dat$x+20000)),
           ylim = c(min(dat$y-20000), max(dat$y+20000)), expand = FALSE) +
  geom_point(data = dat, aes(x, y, color=ac), size=1, alpha = 0.6) +
  geom_point(data = ac.coords, aes(x, y, fill = ac.ns), color = "black", size = 3, pch = 21,
             stroke = 1, alpha = 0.8) +
  scale_color_viridis_c("Activity Center") +
  scale_fill_viridis_c(guide = F) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  facet_wrap(~id)


# Distance decay surface
grid<- create.grid(dat=dat, crs=utm.crs, extent=extent, res=res, buffer=buffer)
grid[]<- 0

#pre-calculate distances between each potential AC location (possib.ac) and each actual location in our data (grid.coord)
dist.mat=GetDistance(AcCoord=data.matrix(ac.coords[,-3]),GridCoord=data.matrix(grid.coord[,-3]), 
                     Ngrid=nrow(grid.coord), Nac=nrow(ac.coords))

prob=exp(-phi*dist.mat)
prob=prob/rowSums(prob)

#plot distance decay surface per AC
par(mfrow = c(4,5))
for (i in 1:nrow(ac.coords)) {
  grid[grid.coord.med$grid.cell]<- prob[,i]  #change prob column to explore each AC
  plot(grid, axes = F, useRaster = F)
  points(ac.coords$x, ac.coords$y, pch = 16)
  par(mfrow = c(1,1))


## AC Heatmap by month/week/day and year
dat.list<- df.to.list(dat)
plot.ac.spatemp(data = dat.list, ac.coords = ac.coords, units = day)




###################
### Save Output ###
###################

# write.csv(ac.coords, "Activity Center Coordinates_TOHO.csv", row.names = F)
# write.csv(dat, "Snail Kite Gridded Data_AC_TOHO.csv", row.names = F)
