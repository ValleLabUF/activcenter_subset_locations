library(microbenchmark)

microbenchmark(
  ac.ind.tmp=sample.ac(ac.ind=ac.ind,dat=dat,theta=theta,n.ac=n.ac,n.grid=n.grid,phi=phi,
                       dist.mat=dist.mat,n.tsegm=n.tsegm,n.possib.ac=n.possib.ac),
  phi.tmp=sample.phi(ac.ind=ac.ind,dist.mat=dist.mat,n.grid=n.grid, n.ac=n.ac,phi=phi,jump=jump1$phi,
                     dat=dat,n.tsegm=n.tsegm,theta=theta),
  z=sample.z(ac.ind=ac.ind,dist.mat=dist.mat,n.grid=n.grid,n.ac=n.ac,n.tsegm=n.tsegm,
             dat=dat,phi=phi,log.theta=log(theta)),
  v=sample.v(z=z,n.ac=n.ac,gamma1=gamma1),  
  logl=get.loglikel(z=z,dist.mat=dist.mat,phi=phi,dat=dat,ac.ind=ac.ind,n.grid=n.grid)
)