gibbs.activity.center=function(dat,grid.coord,n.ac,ac.coord.init,gamma1,possib.ac){
  #basic setup
  n.tsegm=nrow(dat)
  n.grid=nrow(grid.coord)
  grid.coord=data.matrix(grid.coord)
  n.possib.ac=nrow(possib.ac)
  
  #initial values
  ac.ind=sample(n.grid,size=n.ac)
  z=sample(1:n.ac,size=n.tsegm,replace=T) #cluster membership
  phi=0.0001 #distance decay parameter
  theta=rep(1/n.ac,n.ac)
  
  #matrices to store results
  store.coord=matrix(NA,ngibbs,n.ac*2)
  store.z=matrix(NA,ngibbs,n.tsegm)
  store.param=matrix(NA,ngibbs,1) #to store phi
  store.logl=matrix(NA,ngibbs,1)
  store.theta=matrix(NA,ngibbs,n.ac)
  
  #MH stuff
  adaptMH=50
  jump1=list(phi=0.2,ac.ind=rep(1,n.ac))
  accept1=list(phi=0,ac.ind=rep(0,n.ac))
  
  #pre-calculate distances between each potential AC location (possib.ac) and each actual location in our data (grid.coord)
  dist.mat=GetDistance(AcCoord=data.matrix(possib.ac),GridCoord=data.matrix(grid.coord), 
                       Ngrid=nrow(grid.coord), Nac=nrow(possib.ac))
  
  #progress bar
  pb <- progress_bar$new(
    format = " iteration (:current/:total) [:bar] :percent [Elapsed: :elapsed, Remaining: :eta]",
    total = ngibbs, clear = FALSE, width= 100)
  
  #gibbs sampler
  for (i in 1:ngibbs){
    pb$tick()  #create progress bar
    
    #sample AC
    tmp=sample.ac(ac.ind=ac.ind,dat=dat,theta=theta,n.ac=n.ac,n.grid=n.grid,phi=phi,
                  dist.mat=dist.mat,n.tsegm=n.tsegm,n.possib.ac=n.possib.ac)
    ac.ind=tmp$ac.ind
    accept1$ac.ind=accept1$ac.ind+tmp$accept
    # ac.coord=ac.coord.true
    
    #sample phi
    tmp=sample.phi(ac.ind=ac.ind,dist.mat=dist.mat,n.grid=n.grid,
                   n.ac=n.ac,phi=phi,jump=jump1$phi,dat=dat,n.tsegm=n.tsegm,theta=theta)
    phi=tmp$phi
    accept1$phi=accept1$phi+tmp$accept
    # phi=phi.true
    
    #sample z
    z=sample.z(ac.ind=ac.ind,dist.mat=dist.mat,n.grid=n.grid,n.ac=n.ac,n.tsegm=n.tsegm,
               dat=dat,phi=phi,log.theta=log(theta))
    # z=z.true
    
    #sample theta
    v=sample.v(z=z,n.ac=n.ac,gamma1=gamma1)
    # theta=rep(NA,n.ac)
    # theta[1]=v[1]
    # tmp=(1-v[1])
    # for (j in 2:n.ac){
    #   theta[j]=v[j]*tmp
    #   tmp=tmp*(1-v[j])
    # } 
    theta=GetTheta(v=v,nac=n.ac)
    
    #get loglikel
    logl=get.loglikel(z=z,dist.mat=dist.mat,phi=phi,dat=dat,ac.ind=ac.ind,n.grid=n.grid)
    
    if (i<nburn & i%%adaptMH==0){
      #adapt MH
      tmp=print.adapt(accept1z=accept1,jump1z=jump1,accept.output=adaptMH)
      jump1=tmp$jump1
      accept1=tmp$accept1
      
      #re-order data from time to time according to theta (largest to smallest)
      ordem=order(theta,decreasing=T)
      znew=rep(NA,n.tsegm)
      ac.ind=ac.ind[ordem]
      for (j in 1:n.ac){
        cond=z==ordem[j]
        if (sum(cond)>0) znew[cond]=j
      }
      z=znew 
      theta=theta[ordem]
    }
    
    #store results
    store.coord[i,]=unlist(grid.coord[ac.ind,])
    store.z[i,]=z
    store.param[i,]=phi
    store.logl[i,]=logl
    store.theta[i,]=theta
  }
  list(coord=store.coord,z=store.z,phi=store.param,logl=store.logl,theta=store.theta)  
}