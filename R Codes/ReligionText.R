###################################
#     Apr 28, 2022  by Rong Ma
#   Visualizing Religious Texts
###################################

install.packages(c("cluster","ggplot2","uwot","ggrepel","Rtsne","phateR"))
library(cluster)
library(ggplot2)
library(uwot)
library(ggrepel)
library(Rtsne)
library(phateR)

###################################
######## Load data 
###################################

data = read.csv("~/Download/Data/AllBooks_baseline_DTM_Labelled.csv")

info = data[,1]
info = gsub("\\_.*", "",info)
data = data[,-1]

dim(data)
table(info)

data = data[,which(colSums(data)!=0)]
data = scale(data, center=TRUE, scale = TRUE)
dim(data)
length(info)
n=dim(data)[1]
info=factor(info)
levels(info) = c("BOE1", "BOE2", "BOP", "BOW", "BUD", "TTC", "UPA","YOG")

################### ensemble visualization


candidate.out = candidate.visual(data, method=c("PCA", "MDS", "iMDS", "Sammon", "LLE", "HLLE", "Isomap", 
                                                "kPCA", "LEIM", "UMAP", "tSNE","PHATE"),
                                 kpca.sigma = c(0.01, 0.001), 
                                 umap.k= c(30, 50), 
                                 tsne.perplexity = c(10, 50),
                                 phate.k = c(30, 50))

ensemble.out = ensemble.v.local(data.list=candidate.out[[1]], name.method = candidate.out[[2]])

setwd("~/meta-visualization/output/text") # set directory for the output
save(candidate.out, ensemble.out, file="text.RData")



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
ensemble.data=umap(as.dist(ensemble.out[[1]]),  n_neighbors = 30)
sindex.ensemble=silhouette(as.numeric(factor(info)),dmatrix = as.matrix(dist(ensemble.data)))[,3]
data.mat= as.matrix(dist(data))
sindex.data=silhouette(as.numeric(factor(info)),dmatrix = data.mat)[,3]
naive.data =  umap(as.dist( apply(ensemble.out[[4]], c(2,3),mean)),  n_neighbors = 30) 
sindex.naive = silhouette(as.numeric(factor(info)),dmatrix = as.matrix(dist(naive.data)))[,3]
sindex.ori = silhouette(as.numeric(factor(info)),dmatrix = as.matrix(data.mat))[,3]
sindex.ext = cbind(sindex, sindex.ensemble, sindex.naive, sindex.ori)
colnames(sindex.ext)=c(candidate.out[[2]], "meta-spec", "meta-aver", "ori")
data.plot=data.frame(sindex = c(sindex.ext), method = rep(c(candidate.out[[2]],"meta-spec", "meta-aver", "ori"), each=n))
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

corr=c()
for(i in 1:n){
  corr[i]=cor(sindex[i,], ensemble.out[[2]][i,])
}
summary(corr)


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
ggplot(data.plot, aes(x=dim1, y=dim2)) + geom_point(size=1.5, aes(shape=cluster, color=cluster)) +
  scale_shape_manual(values=1:nlevels(data.plot$cluster)) + ggtitle(title.name[k])
dev.off()

#visualize eigenscores
k=0
k=k+1
data.plot = data.frame(dim1=candidate.out[[1]][[k]][,1], dim2=candidate.out[[1]][[k]][,2], 
                       cluster=factor(info), eigenscore = c(ensemble.out[[2]][,k]), sindex = sindex[,k])
pdf(file = paste0("cluster","-",candidate.out[[2]][k],"-eigenscore.pdf"), width=5,height=4) 
ggplot(data.plot, aes(x=dim1, y=dim2)) + geom_point(size=1.5, aes(shape=cluster, color=eigenscore)) +
  scale_shape_manual(values=1:nlevels(data.plot$cluster)) + ggtitle(title.name[k])
dev.off()
pdf(file = paste0("cluster","-",candidate.out[[2]][k],"-eigenscore-box.pdf"), width=4,height=4) 
ggplot(data.plot, aes(x=cluster, y=eigenscore)) + geom_boxplot() + ggtitle(title.name[k])+ 
ylab("eigenscore") + xlab("cluster")+ylim(0.1,0.35)
dev.off()
pdf(file = paste0("cluster","-",candidate.out[[2]][k],"-sindex-box.pdf"), width=4,height=4) 
ggplot(data.plot, aes(x=cluster, y=sindex)) + geom_boxplot() + ggtitle(title.name[k])+ 
  ylab("silhouette index") + xlab("cluster")+ylim(-1,0.5)
dev.off()


############## meta-visualizations
data.plot = data.frame(dim1=ensemble.data[,1], dim2=ensemble.data[,2], cluster=factor(info))
pdf(file = paste0("cluster-META-spec.pdf"), width=5,height=4) 
ggplot(data.plot, aes(x=dim1, y=dim2)) + geom_point(size=1.5, aes(shape=cluster, color=cluster)) +
  scale_shape_manual(values=1:nlevels(data.plot$cluster))+ ggtitle("meta-spec")
dev.off()

data.plot = data.frame(dim1=naive.data[,1], dim2=naive.data[,2], cluster=factor(info))
pdf(file = paste0("cluster-META-aver.pdf"), width=5,height=4) 
ggplot(data.plot, aes(x=dim1, y=dim2)) + geom_point(size=1.5, aes(shape=cluster, color=cluster)) +
  scale_shape_manual(values=1:nlevels(data.plot$cluster))+ ggtitle("meta-apec")
dev.off()

############# bar plot with error
sindex.naive = c()
sindex.ensemble=c()
for(i in 1:30){
  ensemble.data=umap(as.dist(ensemble.out[[1]]),  n_neighbors = 30) 
  sindex.ensemble[i]=median(silhouette(as.numeric(factor(info)),dmatrix = as.matrix(dist(ensemble.data)))[,3])
  
  naive.data = umap(as.dist( apply(ensemble.out[[4]], c(2,3),mean)),  n_neighbors = 30)
  sindex.naive[i] = median(silhouette(as.numeric(factor(info)),dmatrix = as.matrix(dist(naive.data)))[,3])
  print(i)
}

data.plot = data.frame(sindex = c(apply(sindex, 2, median), mean(sindex.naive), mean(sindex.ensemble)), 
                       method = c(candidate.out[[2]], "meta-aver","meta-spec"),
                       low = c(rep(NA,16), quantile(sindex.naive,0.025),  quantile(sindex.ensemble,0.025)),
                       up = c(rep(NA,16), quantile(sindex.naive,0.975),  quantile(sindex.ensemble,0.975)))
pdf(file = "cluster-errorbar-small.pdf", width=4,height=3.5) 
ggplot(data.plot, aes(x=reorder(method, sindex), y=sindex)) +
  geom_col() + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  ylab("median silhouette index") + xlab("method") + geom_errorbar(aes(ymin=low, ymax=up), colour="red", width=0.5)
dev.off()

##### boxplot of sindex: t-SNE for meta-viz

sindex.naive = c()
sindex.ensemble=c()
for(i in 1:30){
  ensemble.data=Rtsne(as.dist(ensemble.out[[1]]),  perplexity = 50)$Y
  sindex.ensemble[i]=median(silhouette(as.numeric(factor(info)),dmatrix = as.matrix(dist(ensemble.data)))[,3])
  
  naive.data =  Rtsne(as.dist( apply(ensemble.out[[4]], c(2,3),mean)),  perplexity = 50)$Y 
  sindex.naive[i] = median(silhouette(as.numeric(factor(info)),dmatrix = as.matrix(dist(naive.data)))[,3])
  print(i)
}

data.plot = data.frame(sindex = c(apply(sindex, 2, median), mean(sindex.naive), mean(sindex.ensemble)), 
                       method = c(candidate.out[[2]], "meta-aver(tSNE)","meta-spec(tSNE)"),
                       low = c(rep(NA,16), quantile(sindex.naive,0.025),  quantile(sindex.ensemble,0.025)),
                       up = c(rep(NA,16), quantile(sindex.naive,0.975),  quantile(sindex.ensemble,0.975)))
pdf(file = "cluster-errorbar2.pdf", width=4,height=3.5) 
ggplot(data.plot, aes(x=reorder(method, sindex), y=sindex)) +
  geom_col() + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  ylab("median silhouette index") + xlab("method") + geom_errorbar(aes(ymin=low, ymax=up), colour="red", width=0.5)
dev.off()

##### boxplot of sindex: phate for meta-viz

sindex.naive = c()
sindex.ensemble=c()
for(i in 1:30){
  ensemble.data=phate(as.dist(ensemble.out[[1]]),  knn = 30, ndim=2)$embedding
  sindex.ensemble[i]=median(silhouette(as.numeric(factor(info)),dmatrix = as.matrix(dist(ensemble.data)))[,3])
  
  naive.data =  phate(as.dist( apply(ensemble.out[[4]], c(2,3),mean)),  knn = 30, ndim=2)$embedding
  sindex.naive[i] = median(silhouette(as.numeric(factor(info)),dmatrix = as.matrix(dist(naive.data)))[,3])
  print(i)
}

data.plot = data.frame(sindex = c(apply(sindex, 2, median), mean(sindex.naive), mean(sindex.ensemble)), 
                       method = c(candidate.out[[2]], "meta-aver(PHATE)","meta-spec(PHATE)"),
                       low = c(rep(NA,16), quantile(sindex.naive,0.025),  quantile(sindex.ensemble,0.025)),
                       up = c(rep(NA,16), quantile(sindex.naive,0.975),  quantile(sindex.ensemble,0.975)))
pdf(file = "cluster-errorbar3.pdf", width=4,height=3.5) 
ggplot(data.plot, aes(x=reorder(method, sindex), y=sindex)) +
  geom_col() + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  ylab("median silhouette index") + xlab("method") + geom_errorbar(aes(ymin=low, ymax=up), colour="red", width=0.5)
dev.off()
