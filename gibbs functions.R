#calculates part of the marginal loglikelihood after integrating our the z's
get.calc.mloglik=function(dist.mat.sel,phi,dat,n.tsegm,n.ac,n.grid,log.theta){
  prob=exp(-phi*dist.mat.sel)
  lprob=log(prob/rowSums(prob))
  
  # res=matrix(NA,n.tsegm,n.ac)
  # for (i in 1:n.ac){
  #   lprob.mat=matrix(lprob[i,],n.tsegm,n.grid,byrow=T)
  #   res[,i]=rowSums(dat*lprob.mat)+log.theta[i]
  # }
  # res
  
  res1=mloglikel(ntsegm=n.tsegm, nac=n.ac, ngrid=n.grid, lprob=lprob, dat=dat,
                 LogTheta=log.theta)
  res1
}
#-------------------------------------------------------------
sample.ac=function(ac.ind,dat,theta,n.ac,n.grid,phi,dist.mat,n.tsegm,n.possib.ac){
  ac.ind.orig=ac.ind.old=ac.ind
  
  for (i in 1:n.ac){
    ac.ind.new=ac.ind.old
    ind=sample(1:n.possib.ac,size=1)
    ac.ind.new[i]=ind
    
    #get marginal loglikel
    tmp.old=get.calc.mloglik(dist.mat.sel=dist.mat[ac.ind.old,],phi=phi,dat=dat,
                             n.tsegm=n.tsegm,n.ac=n.ac,n.grid=n.grid,log.theta=log(theta))
    tmp.new=get.calc.mloglik(dist.mat.sel=dist.mat[ac.ind.new,],phi=phi,dat=dat,
                             n.tsegm=n.tsegm,n.ac=n.ac,n.grid=n.grid,log.theta=log(theta))
    # max1=-apply(cbind(tmp.old,tmp.new),1,max) 
    max1=-GetMaxRows(mat=cbind(tmp.old,tmp.new))
    tmp.old1=tmp.old+max1
    tmp.new1=tmp.new+max1
    pold=sum(log(rowSums(exp(tmp.old1))))
    pnew=sum(log(rowSums(exp(tmp.new1))))
    
    #accept or reject MH
    k=acceptMH(p0=pold,p1=pnew,x0=ac.ind.old[i],x1=ac.ind.new[i],BLOCK=F)
    ac.ind.old[i]=k$x
  }
  list(ac.ind=ac.ind.old,accept=ac.ind.orig!=ac.ind.old)
}
#-----------------------------------
sample.phi=function(ac.ind,dist.mat,n.grid,n.ac,phi,jump,dat,n.tsegm,theta){
  old=phi
  new=abs(rnorm(1,mean=old,sd=jump)) #reflection proposal around zero
  
  #get marginal loglikel
  dist.mat1=dist.mat[ac.ind,]
  log.theta=log(theta)
  tmp.old=get.calc.mloglik(dist.mat.sel=dist.mat1,phi=old,dat=dat,
                           n.tsegm=n.tsegm,n.ac=n.ac,n.grid=n.grid,log.theta=log.theta)
  tmp.new=get.calc.mloglik(dist.mat.sel=dist.mat1,phi=new,dat=dat,
                           n.tsegm=n.tsegm,n.ac=n.ac,n.grid=n.grid,log.theta=log.theta)
  max1=GetMaxRows(mat=cbind(tmp.old,tmp.new))
  tmp.old1=tmp.old-max1
  tmp.new1=tmp.new-max1
  pold=sum(log(rowSums(exp(tmp.old1))))
  pnew=sum(log(rowSums(exp(tmp.new1))))
  
  #accept or reject MH
  k=acceptMH(p0=pold,p1=pnew,x0=old,x1=new,BLOCK=F)  
  logl=ifelse(k$accept==1,pnew,pold)
  list(phi=k$x,accept=k$accept)
}
#-----------------------------------
sample.z=function(ac.ind,dist.mat,n.grid,n.ac,n.tsegm,dat,phi,log.theta){
  #get distance
  dist1=dist.mat[ac.ind,]
  prob=exp(-phi*dist1)
  prob=prob/rowSums(prob)
  lprob=log(prob)
  
  #get loglikel
  logl=mloglikel(ntsegm=n.tsegm, nac=n.ac, ngrid=n.grid, lprob=lprob, dat=dat,
                 LogTheta=log.theta)
  
  # logl=matrix(NA,n.tsegm,n.ac)
  # for (i in 1:n.ac){
  #   lprob1=matrix(lprob[i,],n.tsegm,n.grid,byrow=T)
  #   logl[,i]=rowSums(dat*lprob1)+log.theta[i]    
  # }
  maximo=GetMaxRows(mat=logl)
  logl=logl-maximo
  tmp=exp(logl)
  prob=tmp/rowSums(tmp)
  
  #sample from multinomial
  z=rmultinom1(prob=prob,randu=runif(n.tsegm))
  z+1
}
#-----------------------------------
acceptMH <- function(p0,p1,x0,x1,BLOCK){   #accept for M, M-H
  # if BLOCK, then accept as a block,
  # otherwise, accept individually
  
  nz           <- length(x0)  #no. to accept
  if(BLOCK) nz <- 1
  
  a    <- exp(p1 - p0)       #acceptance PR
  z    <- runif(nz,0,1)
  keep <- which(z < a)
  
  if(BLOCK & length(keep) > 0) x0 <- x1
  if(!BLOCK)                   x0[keep] <- x1[keep]           
  accept <- length(keep)        
  
  list(x = x0, accept = accept)
}
#----------------------------
print.adapt = function(accept1z,jump1z,accept.output){
  accept1=accept1z; jump1=jump1z; 
  
  for (k in 1:length(accept1)){
    z=accept1[[k]]/accept.output
    print(names(accept1)[k])
    print(mean(z)); print(mean(jump1[[k]]))
  }
  
  for (k in 1:length(jump1)){
    cond=(accept1[[k]]/accept.output)>0.4 & jump1[[k]]<100
    jump1[[k]][cond] = jump1[[k]][cond]*2       
    cond=(accept1[[k]]/accept.output)<0.2 & jump1[[k]]>0.01
    jump1[[k]][cond] = jump1[[k]][cond]*0.5
    accept1[[k]][]=0
  }
  
  return(list(jump1=jump1,accept1=accept1))
}
#----------------------------
sample.v=function(z,n.ac,gamma1){
  tmp=table(z)
  n=rep(0,n.ac)
  n[as.numeric(names(tmp))]=tmp
  
  seq1=n.ac:1
  tmp=cumsum(n[seq1])
  n.ge=tmp[seq1]
  n.ge1=n.ge[-1]
  v=rbeta(n.ac-1,n[-n.ac]+1,n.ge1+gamma1)
  c(v,1)
}
#-----------------------------
get.loglikel=function(z,dist.mat,phi,dat,ac.ind,n.grid){
  #get distance
  dist1=dist.mat[ac.ind,]
  
  #get log-probabilities for multinomial
  prob=exp(-phi*dist1)
  prob=prob/rowSums(prob)
  lprob=log(prob)
  
  #calculate conditional loglikel
  uni.z=unique(z)
  res=0
  for(i in uni.z){
    cond=z==i
    dat1=dat[cond,]
    n=sum(cond)
    lprob1=matrix(lprob[i,],n,n.grid,byrow=T)
    res=res+sum(dat1*lprob1)
  }
  res
}