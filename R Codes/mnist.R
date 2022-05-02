###################################
#     Apr 28, 2022  by Rong Ma
#   A Simple Example to Get Start
###################################

library(cluster)
library(ggplot2)
library(uwot)
library(ggrepel)
library(Rtsne)
library(phateR)
###################################
######## Load mnist data 
###################################
data=read.csv("~/Download/mnist/mnist_test.csv", header =F)

N=400
info = data[,1]
set.seed(20) 
data0=data[c(which(data[,1]==1)[sample(sum(info==1),N)], which(data[,1]==2)[sample(sum(info==2),N)],
            which(data[,1]==3)[sample(sum(info==3),N)], which(data[,1]==4)[sample(sum(info==4),N)], 
             which(data[,1]==5)[sample(sum(info==5),N)]),]


table(data0[,1])
data.nolabel = data0[,-1]
info = data0[,1]
data.nolabel = data.nolabel[,which(colSums(data.nolabel)!=0)] #remove noninformative features
data = scale(data.nolabel, center=TRUE, scale = TRUE)
dim(data)
length(info)
n=dim(data)[1]

########################################
############# ensemble visualization
########################################


#generate 16 candidatet visualizations
candidate.out = candidate.visual(data, method=c("PCA", "MDS", "iMDS", "Sammon", "LLE", "HLLE", "Isomap", 
                                                "kPCA", "LEIM", "UMAP", "tSNE","PHATE"),
                                 kpca.sigma = c(0.01, 0.001), 
                                 umap.k= c(100, 200), 
                                 tsne.perplexity = c(10, 100),
                                 phate.k = c(100, 200))


#obtain eigenscores and meta-distance matrix
ensemble.out = ensemble.v.local(data.list=candidate.out[[1]], name.method = candidate.out[[2]])

setwd("meta-visualization/output/mnist")# set directory for the output
save(candidate.out, ensemble.out, file="mnist1-5.RData")

####################################
############ results
####################################

#obtain silhouette indices for candidate viz
K=dim(ensemble.out[[4]])[1]
sindex = matrix(ncol=K,nrow=n)
for(i in 1:K){
  sindex[,i]=silhouette(as.numeric(factor(info)), dmatrix = ensemble.out[[4]][i,,])[,3]
}

#boxplot for eigenscores
data.plot = data.frame(eigen.score = c(ensemble.out[[2]]), method = rep(candidate.out[[2]], each=n))
pdf(file = "cluster-score-box.pdf", width=4,height=4) 
ggplot(data.plot, aes(x=reorder(method, eigen.score, FUN=median), y=eigen.score)) +
  geom_boxplot(outlier.size = 0.5) + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  ylab("eigenscore") + xlab("method")
dev.off()

##### boxplot of sindex
set.seed(30)
ensemble.data=umap(as.dist(ensemble.out[[1]]),  n_neighbors = 200)
sindex.ensemble=silhouette(as.numeric(factor(info)),dmatrix = as.matrix(dist(ensemble.data)))[,3]
data.mat= as.matrix(dist(data))
sindex.data=silhouette(as.numeric(factor(info)),dmatrix = data.mat)[,3]
naive.data =  umap(as.dist( apply(ensemble.out[[4]], c(2,3),mean)),  n_neighbors = 200) 
sindex.naive = silhouette(as.numeric(factor(info)),dmatrix = as.matrix(dist(naive.data)))[,3]
sindex.ori = silhouette(as.numeric(factor(info)),dmatrix = as.matrix(data.mat))[,3]
sindex.ext = cbind(sindex, sindex.ensemble, sindex.naive, sindex.ori)
colnames(sindex.ext)=c(candidate.out[[2]], "meta-Spec", "meta-Aver", "ori")
data.plot=data.frame(sindex = c(sindex.ext), method = rep(c(candidate.out[[2]],"meta-Spec", "meta-Aver", "ori"), each=n))
pdf(file = "cluster-sindex-box.pdf", width=4,height=4) 
ggplot(data.plot, aes(x=reorder(method, sindex, FUN=median), y=sindex)) +
  geom_boxplot(outlier.size = 0.5) + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) + 
  ylab("silhouette index") + xlab("method")
dev.off()


#### sindex vs eigen-score
data.plot = data.frame(eigenscore = apply(ensemble.out[[2]],2, mean), sindex = colMeans(sindex), method = candidate.out[[2]])
pdf(file = "cluster-score-sindex.pdf", width=4*1.2,height=3*1.4) 
ggplot(data.plot, aes(x=eigenscore, y=sindex)) +
  geom_smooth(method=lm)+
  geom_point(size=2) + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) + 
  ylab("average silhouette index") + xlab("average eigenscore")+geom_label_repel(aes(label = method),
                                                                                 force         = 70,
                                                                                 max.overlaps = 30,
                                                                                 segment.color = 'grey50')
dev.off()
cor(colMeans(sindex), colMeans(ensemble.out[[2]]))
#0.8839462

#### why not average? why not raw?
data.plot=data.frame( sindex = c(sindex.data, sindex.naive, sindex.ensemble), 
                      method = rep(c("None", "ENSB-Aver", "ENSB-Spec"), each=n))
pdf(file = "cluster-comp3-sindex.pdf", width=3,height=3) 
ggplot(data.plot, aes(x=reorder(method, sindex), y=sindex, color=method)) +
  geom_boxplot() + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) + 
  ylab("Silhouette Index") + xlab("Method") + theme(legend.position="none")
dev.off()


### visualize candidate visualizations
title.name = candidate.out[[2]]
title.name[8]="kPCA"
title.name[9]="kPCA"
title.name[11]="UMAP"
title.name[12]="UMAP"
title.name[13]="tSNE"
title.name[14]="tSNE"
title.name[15]="PHATE"
title.name[16]="PHATE"
k=0
#repeat the following K times
k=k+1
data.plot = data.frame(dim1=candidate.out[[1]][[k]][,1], dim2=candidate.out[[1]][[k]][,2], cluster=factor(info))
pdf(file = paste0("cluster","-",candidate.out[[2]][k],".pdf"), width=5,height=4) 
ggplot(data.plot, aes(x=dim1, y=dim2, color = cluster)) + geom_point(size=0.5) + ggtitle(title.name[k])
dev.off()


############## meta-visualizations
data.plot = data.frame(dim1=ensemble.data[,1], dim2=ensemble.data[,2], cluster=factor(info))
pdf(file = paste0("cluster-META-spec.pdf"), width=5,height=4) 
ggplot(data.plot, aes(x=dim1, y=dim2)) + geom_point(size=0.5, aes(shape=cluster, color=cluster)) +scale_shape_manual(values=1:nlevels(data.plot$cluster))+ ggtitle("meta-Spec")
dev.off()

data.plot = data.frame(dim1=naive.data[,1], dim2=naive.data[,2], cluster=factor(info))
pdf(file = paste0("cluster-META-aver.pdf"), width=5,height=4) 
ggplot(data.plot, aes(x=dim1, y=dim2)) + geom_point(size=0.5, aes(shape=cluster, color=cluster)) +scale_shape_manual(values=1:nlevels(data.plot$cluster))+ ggtitle("meta-Spec")
dev.off()

############# bar plot with error
sindex.naive = c()
sindex.ensemble=c()
for(i in 1:30){
  ensemble.data=umap(as.dist(ensemble.out[[1]]),  n_neighbors = 200) 
  sindex.ensemble[i]=median(silhouette(as.numeric(factor(info)),dmatrix = as.matrix(dist(ensemble.data)))[,3])
  
  naive.data = umap(as.dist( apply(ensemble.out[[4]], c(2,3),mean)),  n_neighbors = 200)
  sindex.naive[i] = median(silhouette(as.numeric(factor(info)),dmatrix = as.matrix(dist(naive.data)))[,3])
  print(i)
}

data.plot = data.frame(sindex = c(apply(sindex, 2, median), mean(sindex.naive), mean(sindex.ensemble)), 
                       method = c(candidate.out[[2]], "meta-Aver","meta-Spec"),
                       low = c(rep(NA,16), quantile(sindex.naive,0.025),  quantile(sindex.ensemble,0.025)),
                       up = c(rep(NA,16), quantile(sindex.naive,0.975),  quantile(sindex.ensemble,0.975)))
pdf(file = "cluster-errorbar.pdf", width=4,height=3.5) 
ggplot(data.plot, aes(x=reorder(method, sindex), y=sindex)) +
  geom_col() + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  ylab("median silhouette index") + xlab("method") + geom_errorbar(aes(ymin=low, ymax=up), colour="red", width=0.5)
dev.off()

##### boxplot of sindex: t-SNE for meta-viz
set.seed(33)
ensemble.data=Rtsne(as.dist(ensemble.out[[1]]),  perplexity = 100)$Y
sindex.ensemble=silhouette(as.numeric(factor(info)),dmatrix = as.matrix(dist(ensemble.data)))[,3]
data.mat= as.matrix(dist(data))
sindex.data=silhouette(as.numeric(factor(info)),dmatrix = data.mat)[,3]
naive.data =  Rtsne(as.dist( apply(ensemble.out[[4]], c(2,3),mean)),  perplexity = 100)$Y 
sindex.naive = silhouette(as.numeric(factor(info)),dmatrix = as.matrix(dist(naive.data)))[,3]
sindex.ori = silhouette(as.numeric(factor(info)),dmatrix = as.matrix(data.mat))[,3]
sindex.ext = cbind(sindex, sindex.ensemble, sindex.naive, sindex.ori)
colnames(sindex.ext)=c(candidate.out[[2]], "meta-Spec", "meta-Aver", "ori")
data.plot=data.frame(sindex = c(sindex.ext), method = rep(c(candidate.out[[2]],"meta-Spec", "meta-Aver", "ori"), each=n))
pdf(file = "cluster-sindex-box-2.pdf", width=4,height=4) 
ggplot(data.plot, aes(x=reorder(method, sindex, FUN=median), y=sindex)) +
  geom_boxplot(outlier.size = 0.5) + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) + 
  ylab("silhouette index") + xlab("method")
dev.off()


sindex.naive = c()
sindex.ensemble=c()
for(i in 1:30){
  ensemble.data=Rtsne(as.dist(ensemble.out[[1]]),  perplexity = 100)$Y
  sindex.ensemble[i]=median(silhouette(as.numeric(factor(info)),dmatrix = as.matrix(dist(ensemble.data)))[,3])
  
  naive.data =  Rtsne(as.dist( apply(ensemble.out[[4]], c(2,3),mean)),  perplexity = 100)$Y 
  sindex.naive[i] = median(silhouette(as.numeric(factor(info)),dmatrix = as.matrix(dist(naive.data)))[,3])
  print(i)
}

data.plot = data.frame(sindex = c(apply(sindex, 2, median), mean(sindex.naive), mean(sindex.ensemble)), 
                       method = c(candidate.out[[2]], "meta-Aver(tSNE)","meta-Spec(tSNE)"),
                       low = c(rep(NA,16), quantile(sindex.naive,0.025),  quantile(sindex.ensemble,0.025)),
                       up = c(rep(NA,16), quantile(sindex.naive,0.975),  quantile(sindex.ensemble,0.975)))
pdf(file = "cluster-errorbar2.pdf", width=4,height=3.5) 
ggplot(data.plot, aes(x=reorder(method, sindex), y=sindex)) +
  geom_col() + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  ylab("median silhouette index") + xlab("method") + geom_errorbar(aes(ymin=low, ymax=up), colour="red", width=0.5)
dev.off()

##### boxplot of sindex: phate for meta-viz

set.seed(20)
ensemble.data=phate(as.dist(ensemble.out[[1]]),  knn = 100, ndim=2)$embedding
sindex.ensemble=silhouette(as.numeric(factor(info)),dmatrix = as.matrix(dist(ensemble.data)))[,3]
data.mat= as.matrix(dist(data))
sindex.data=silhouette(as.numeric(factor(info)),dmatrix = data.mat)[,3]
naive.data =  phate(as.dist( apply(ensemble.out[[4]], c(2,3),mean)),  knn = 100, ndim=2)$embedding
sindex.naive = silhouette(as.numeric(factor(info)),dmatrix = as.matrix(dist(naive.data)))[,3]
sindex.ori = silhouette(as.numeric(factor(info)),dmatrix = as.matrix(data.mat))[,3]
sindex.ext = cbind(sindex, sindex.ensemble, sindex.naive, sindex.ori)
colnames(sindex.ext)=c(candidate.out[[2]], "meta-Spec", "meta-Aver", "ori")
data.plot=data.frame(sindex = c(sindex.ext), method = rep(c(candidate.out[[2]],"meta-Spec", "meta-Aver", "ori"), each=n))
pdf(file = "cluster-sindex-box-3.pdf", width=4,height=4) 
ggplot(data.plot, aes(x=reorder(method, sindex, FUN=median), y=sindex)) +
  geom_boxplot(outlier.size = 0.5) + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) + 
  ylab("silhouette index") + xlab("method")
dev.off()


sindex.naive = c()
sindex.ensemble=c()
for(i in 1:30){
  ensemble.data=phate(as.dist(ensemble.out[[1]]),  knn = 100, ndim=2)$embedding
  sindex.ensemble[i]=median(silhouette(as.numeric(factor(info)),dmatrix = as.matrix(dist(ensemble.data)))[,3])
  
  naive.data =  phate(as.dist( apply(ensemble.out[[4]], c(2,3),mean)),  knn = 100, ndim=2)$embedding
  sindex.naive[i] = median(silhouette(as.numeric(factor(info)),dmatrix = as.matrix(dist(naive.data)))[,3])
  print(i)
}

data.plot = data.frame(sindex = c(apply(sindex, 2, median), mean(sindex.naive), mean(sindex.ensemble)), 
                       method = c(candidate.out[[2]], "meta-Aver(PHATE)","meta-Spec(PHATE)"),
                       low = c(rep(NA,16), quantile(sindex.naive,0.025),  quantile(sindex.ensemble,0.025)),
                       up = c(rep(NA,16), quantile(sindex.naive,0.975),  quantile(sindex.ensemble,0.975)))
pdf(file = "cluster-errorbar3.pdf", width=4,height=3.5) 
ggplot(data.plot, aes(x=reorder(method, sindex), y=sindex)) +
  geom_col() + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  ylab("median silhouette index") + xlab("method") + geom_errorbar(aes(ymin=low, ymax=up), colour="red", width=0.5)
dev.off()
