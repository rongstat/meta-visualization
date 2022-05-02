###################################
#     May 1, 2022  by Rong Ma
#           Simulations
###################################

library(cluster)
library(ggplot2)
library(uwot)
library(ggrepel)
library(Rtsne)
library(phateR)
library(Rfast)


angle <- function(x,y){
  return(sum(x*y)/sqrt(sum(x^2))/sqrt(sum(y^2)))
}

#################### Gaussian mixture setting
r=6
n=900
p=500
c.prop = c(0.1,0.15,0.1,0.1,0.25,0.3) #class proportion
rho=1 

theta.v=seq(2.5,5,length.out = 20)
s.index=matrix(ncol=length(theta.v), nrow=18)
concord=matrix(ncol=length(theta.v), nrow=18)
score.concord.cor=c()
for(i in 1:length(theta.v)){
  theta=theta.v[i] #signal strength
  
  #set mean vectors
  mu=rep(0,p)
  mu1=mu
  mu1[1]=rho
  mu2=mu
  mu2[2]=rho
  mu3=mu
  mu3[3]=rho
  mu4=mu
  mu4[4]=rho
  mu5=mu
  mu5[5]=rho
  mu6=mu
  mu6[6]=rho
  
  #generate 6 clusters
  data1=matrix(rep(mu1, n*c.prop[1]), nrow=n*c.prop[1], byrow=TRUE)
  data2=matrix(rep(mu2, n*c.prop[2]), nrow=n*c.prop[2], byrow=TRUE)
  data3=matrix(rep(mu3, n*c.prop[3]), nrow=n*c.prop[3], byrow=TRUE)
  data4=matrix(rep(mu4, n*c.prop[4]), nrow=n*c.prop[4], byrow=TRUE)
  data5=matrix(rep(mu5, n*c.prop[5]), nrow=n*c.prop[5], byrow=TRUE)
  data6=matrix(rep(mu6, n*c.prop[6]), nrow=n*c.prop[6], byrow=TRUE)
  # data set X
  data.X=rbind(data1,data2,data3,data4,data5,data6)
  data.Z=rnorm(p*n)
  data.Y = data.Z+data.X*theta
  info=rep(1:6, times = n*c.prop)
  
  dist.ori = as.matrix(dist(data.X))
  
  
  candidate.out = candidate.visual(data.Y, method=c("PCA", "MDS", "iMDS", "Sammon", "LLE","HLLE", "Isomap", 
                                                    "kPCA", "LEIM", "UMAP", "tSNE", "PHATE"))
  ensemble.out = ensemble.v.local(data.list=candidate.out[[1]], name.method = candidate.out[[2]])
  
  
  

  ###########
  
  K=dim(ensemble.out[[4]])[1]
  for(k in 1:K){
    concord[k,i] = angle(c(as.matrix(ensemble.out[[4]][k,,])), c(dist.ori))
    s.index[k,i] = summary(silhouette(as.numeric(factor(info)), dist = as.matrix(ensemble.out[[4]][k,,]), FUN=mean))$avg.width
    
  }
  
  
  concord[K+1,i] = angle(as.matrix(ensemble.out[[1]]), dist.ori)
  ensemble.data= umap(as.dist(ensemble.out[[1]]),  n_neighbors =30, min_dist = 0.1)
  s.index[K+1,i] = summary(silhouette(as.numeric(factor(info)), dist = as.matrix(dist(ensemble.data)), FUN=mean))$avg.width
  
  concord[K+2,i] = angle(as.matrix(apply(ensemble.out[[4]], c(2,3),mean)), dist.ori)
  naive.data=umap(as.dist(apply(ensemble.out[[4]], c(2,3),mean)),  n_neighbors =10)
  s.index[K+2,i] = summary(silhouette(as.numeric(factor(info)), dist = as.matrix(dist(naive.data)), FUN=mean))$avg.width
  score.concord.cor[i]=cor((concord[1:16,i]), apply(ensemble.out[[2]],2, mean))
}

setwd("~/Dropbox/Eric-Rong-James/ensemble-visual/Figures/simu")
data.plot = data.frame(concord = as.vector(t(concord)),
                  method = rep(c(candidate.out[[2]], "meta-spec", "meta-aver"), 
                               each=length(theta.v)),
                  theta = rep(theta.v, times = 18))
data.plot=data.plot[data.plot$method!="HLLE",]
data.plot$method <- with(data.plot, reorder(method, -concord, FUN=mean))
pdf(file = "simu-cluster-concord.pdf", width=4,height=4) 
#pdf(file = "simu-cluster-concord2.pdf", width=6,height=6) 
ggplot(data.plot, aes(x=theta, y=concord, group=factor(method))) +
  geom_line(aes(color=method), size=rep(c(0.8,0.8, rep(0.2,7),   rep(0.2,8) ), each=length(theta.v))) + 
  geom_point(data = subset(data.plot, method == "meta-spec" | method == "meta-aver"), aes(shape=method))+ 
  ylab("concordance")
dev.off()


####making plots: only use for theta=6

k=0
k=k+1
concord=c()
for(i in 1:n){
  concord[i]=angle(c(as.matrix(ensemble.out[[4]][k,i,])), c(dist.ori[i,]))
}
data.plot = data.frame(cluster=info, concord= concord)
pdf(file = paste0("cluster","-",candidate.out[[2]][k],".pdf"), width=3,height=3)
ggplot(data.plot, aes(x=factor(info), y=concord)) + geom_boxplot() + ggtitle(candidate.out[[2]][k])+
  ylab("concordance") + xlab("cluster") + ylim(0.84,0.99)
dev.off()

#visualize eigenscores
k=0
k=k+1
concord=c()
for(i in 1:n){
  concord[i]=angle(c(as.matrix(ensemble.out[[4]][k,i,])), c(dist.ori[i,]))
}
data.plot = data.frame(dim1=candidate.out[[1]][[k]][,1], dim2=candidate.out[[1]][[k]][,2], 
                       cluster=factor(info), eigenscore = c(ensemble.out[[2]][,k]), concordance=concord)
pdf(file = paste0("cluster","-",candidate.out[[2]][k],"-score.pdf"), width=5,height=4) 
ggplot(data.plot, aes(x=dim1, y=dim2)) + geom_point(size=1.5, aes(shape=cluster, color=concordance)) +
  scale_shape_manual(values=1:nlevels(data.plot$cluster)) + ggtitle(title.name[k])
dev.off()
#pdf(file = paste0("cluster","-",candidate.out[[2]][k],"-eigenscore-box.pdf"), width=3,height=3) 
#ggplot(data.plot, aes(x=cluster, y=eigenscore)) + geom_boxplot() + ggtitle(title.name[k])+ 
#  ylab("eigenscore") + xlab("cluster")
#dev.off()

#meta-spec
concord=c()
for(i in 1:n){
  concord[i]=angle(c(as.matrix(ensemble.out[[1]])[i,]), c(dist.ori[i,]))
}
data.plot = data.frame(cluster=info, concordance= concord)
pdf(file = paste0("cluster-meta.pdf"), width=3,height=3)
ggplot(data.plot, aes(x=factor(info), y=concord)) + geom_boxplot() + ggtitle("meta-spec") + 
  ylim(0.85,0.99)+ylab("concordance") + xlab("cluster")
dev.off()
data.plot = data.frame(dim1=ensemble.data[,1], dim2=ensemble.data[,2], 
                                cluster=factor(info), concordance=concord)
concord[1]=0.85
pdf(file = paste0("cluster","-meta","-score.pdf"), width=5,height=4) 
ggplot(data.plot, aes(x=dim1, y=dim2)) + geom_point(size=1.5, aes(shape=cluster, color=concordance)) +
  scale_shape_manual(values=1:nlevels(data.plot$cluster)) + ggtitle("meta-spec") 
dev.off()





######## load smiley 3d
setwd("~/Dropbox/Eric-Rong-James/data/3D datasets/Smiley")
data=read.csv("Smiley_2d.csv")
dim(data)
r=2
plot(data.X, xlab="x", ylab="y")

######## load mammoth
setwd("~/Dropbox/Eric-Rong-James/data/3D datasets/Mammoth")
library("rjson")
r=3
data <- fromJSON(file = "mammoth_3d.json")
data = matrix(unlist(data), ncol=3, byrow=T)

library(scatterplot3d)
scatterplot3d(x=data.X[,3], y=data.X[,1], z=data.X[,2], angle=20)
scatterplot3d(x=data.Y[,3], y=data.Y[,1], z=data.Y[,2], angle=20)

#################### run simulations
angle <- function(x,y){
  return(sum(x*y)/sqrt(sum(x^2))/sqrt(sum(y^2)))
}

n=500
p=300

set.seed(35) 
data.X=data[sample(dim(data)[1],n),]
dim(data.X)
dist.ori = as.matrix(dist(data.X))
theta.v=seq(1,3, length.out = 20) #face
theta.v=seq(0.007,0.012, length.out = 20) #mammoth
concord=matrix(ncol=length(theta.v), nrow=18)
concord.cor=matrix(ncol=length(theta.v), nrow=18)
score.concord.cor=c()
for(i in 1:length(theta.v)){
  theta=theta.v[i] #signal strength
  
  data.Z=rmvnorm(n,rep(0,p),diag(rep(1,p)))
  dim(data.Z)
  data.Y = data.Z
  data.Y[,1:r] = data.Y[,1:r]+as.matrix(data.X*theta)
  
  candidate.out = candidate.visual(data.Y, method=c("PCA", "MDS", "iMDS", "Sammon", "LLE", "HLLE", "Isomap", 
                                            "kPCA", "LEIM", "UMAP", "tSNE", "PHATE"), dim=2)
  ensemble.out = ensemble.v.local(data.list=candidate.out[[1]], name.method = candidate.out[[2]])
  
  K=dim(ensemble.out[[4]])[1]
  for(k in 1:K){
    concord[k,i] = angle(c(as.matrix(ensemble.out[[4]][k,,])), c(dist.ori))
    concord.cor[k,i] = angle(c(as.matrix(ensemble.out[[4]][k,,])), c(dist.ori))
  }
  
  concord[K+1,i] = angle(c(as.matrix(ensemble.out[[1]])),c( dist.ori))
  K.mat = exp(-as.matrix(ensemble.out[[1]])^2/quantile(ensemble.out[[1]],0.5)^2)
  eigen.K.mat = eigen(K.mat)
  U = eigen.K.mat$vectors[,2:3] %*% diag(eigen.K.mat$values[2:3])
  concord.cor[K+1,i] = angle(c(as.matrix(dist(U))),c( dist.ori))
  
  concord[K+2,i] = angle(c(apply(ensemble.out[[4]], c(2,3),mean)), c(dist.ori))
  K.mat = exp(-as.matrix(apply(ensemble.out[[4]], c(2,3),mean))^2/quantile(apply(ensemble.out[[4]], c(2,3),mean),0.5)^2)
  eigen.K.mat = eigen(K.mat)
  U = eigen.K.mat$vectors[,2:3] %*% diag(eigen.K.mat$values[2:3])
  concord.cor[K+2,i] = angle(c(as.matrix(dist(U))), c(dist.ori))
}


setwd("~/Dropbox/Eric-Rong-James/ensemble-visual/Figures/simu")
data.plot = data.frame(concord = as.vector(t(concord)),
                       method = rep(c(candidate.out[[2]], "meta-spec", "meta-aver"), 
                                    each=length(theta.v)),
                       theta = rep(theta.v, times = 18))
data.plot=data.plot[data.plot$method!="HLLE",]
data.plot$method <- with(data.plot, reorder(method, -concord, FUN=mean))
pdf(file = "simu-face-concord.pdf", width=4,height=4) 
#pdf(file = "simu-face-concord2.pdf", width=6,height=6) 
ggplot(data.plot, aes(x=theta, y=concord, group=factor(method))) +
  geom_line(aes(color=method), size=rep(c(0.8,0.8, rep(0.2,8),   rep(0.2,8) ), each=length(theta.v))) + 
  geom_point(data = subset(data.plot, method == "meta-spec" | method == "meta-aver"), aes(shape=method))+ 
  ylab("concordance")
dev.off()


setwd("~/Dropbox/Eric-Rong-James/ensemble-visual/Figures/simu")
data.plot = data.frame(concord = as.vector(t(concord)),
                       method = rep(c(candidate.out[[2]], "meta-spec", "meta-aver"), 
                                    each=length(theta.v)),
                       theta = rep(theta.v, times = 18))
data.plot=data.plot[data.plot$method!="HLLE",]
data.plot$method <- with(data.plot, reorder(method, -concord, FUN=mean))
pdf(file = "simu-mammoth-concord.pdf", width=4,height=4) 
#pdf(file = "simu-mammoth-concord2.pdf", width=6,height=6) 
ggplot(data.plot, aes(x=theta, y=concord, group=factor(method))) +
  geom_line(aes(color=method), size=rep(c(0.8,0.8, rep(0.2,8),   rep(0.2,8) ), each=length(theta.v))) + 
  geom_point(data = subset(data.plot, method == "meta-spec" | method == "meta-aver"), aes(shape=method))+ 
  ylab("concordance")
dev.off()


summary(score.concord.cor)

###### cluster simu stay good for cor and cos. use umap as reprojection
###### smiley face simu tried cos and works. use kpca for reprojection, n=300
###### mammoth data n=500
#summary(score.concord.cor)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.9867  0.9892  0.9905  0.9902  0.9911  0.9929 
0.9654  0.9781  0.9948  0.9882  0.9963  0.9991 
0.9908  0.9941  0.9966  0.9957  0.9972  0.9987 

