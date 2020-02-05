plot(res$logl,type='l')

plot(res$theta[ngibbs,],type='h')

plot(res$phi,type='l')
abline(h=phi.true,col='red')

k=data.frame(estim=res$z[ngibbs,],true1=z.true)
z=table(k);z

ordem=numeric()
for (i in 1:ncol(z)){
  ind=which(z[,i]==max(z[,i]))
  ordem=c(ordem,ind)
}
table(k)[ordem,]

n.ac=20
ac.coord=matrix(res$coord[ngibbs-1,],n.ac,2)
ac.coord[ordem,]
ac.coord.true

rango=range(c(ac.coord.true),ac.coord[ordem,])
plot(ac.coord.true$x,ac.coord[ordem,1],xlim=rango,ylim=rango)
lines(rango,rango,col='red')

plot(ac.coord.true$y,ac.coord[ordem,2],xlim=rango,ylim=rango)
lines(rango,rango,col='red')
