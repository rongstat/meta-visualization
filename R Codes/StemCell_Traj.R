#################################################
#             Apr 29, 2022  by Rong Ma
#   Differentiation of Stem Cell (Trajectory)
################################################

library(cluster)
library(ggplot2)
library(uwot)
library(ggrepel)
library(Rtsne)
library(phateR)
library(data.table)
library(Seurat)

###################################
######## Load data 
###################################


#load RNA seq data
setwd("~/Download/Data/Muscle/mESC-differentiation_hayashi")
data=read.table("GSE98664_tpm_sailfish_mergedGTF_RamDA_mESC_differentiation_time_course.txt", header=T)
rownames(data)=data[,1]
data=data[,-1]
dim(data)

#get cell labels
info=colnames(data)
info=gsub("RamDA_mESC_","", info)
info=substr(info,1,3)
length(info)
table(info)

#normalization and QC
data <- CreateSeuratObject(counts = data, project = "TI", min.cells = 3, min.features = 200)
data <- NormalizeData(data)

#identify highly variable genes
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 500)
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)

#get the normalized data matrix
scale.data = data@assays[["RNA"]]@scale.data
scale.data.var = scale.data[match(data@assays[["RNA"]]@var.features,rownames(scale.data)),]
dim(scale.data.var)
n=dim(scale.data.var)[2]
data.names = factor(info)
levels(data.names)=c(1:5)
data=t(scale.data.var) 
n=dim(data)[1]

################### ensemble visualization


set.seed(29)
candidate.out = candidate.visual(data, kpca.sigma = c(0.002, 0.001), 
                                 umap.k= c(30, 50), 
                                 tsne.perplexity = c(10, 50),
                                 phate.k = c(30, 50),
                                 method=c("PCA", "MDS", "iMDS", "Sammon", "LLE", "HLLE", "Isomap", 
                                                "kPCA", "LEIM", "UMAP", "tSNE","PHATE"))
ensemble.out = ensemble.v.local(data.list=candidate.out[[1]], name.method = candidate.out[[2]])

setwd("~/meta-visualization/output/muscle") # set directory for the output
save(candidate.out, ensemble.out, file="traj.RData")


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
pdf(file = paste0("traj","-",candidate.out[[2]][k],".pdf"), width=5,height=4) 
ggplot(data.plot, aes(x=dim1, y=dim2)) + geom_point(size=1.5, aes(shape=cluster, color=cluster))  + ggtitle(title.name[k])
dev.off()



#obtain silhouette indices for candidate viz
K=dim(ensemble.out[[4]])[1]
sindex = matrix(ncol=K,nrow=n)
for(i in 1:K){
  sindex[,i]=silhouette(as.numeric(factor(info)),dmatrix = ensemble.out[[4]][i,,])[,3]
}
data.plot = data.frame(eigen.score = c(ensemble.out[[2]]), method = rep(candidate.out[[2]], each=n))
pdf(file = "traj-score-box.pdf", width=4,height=4) 
ggplot(data.plot, aes(x=reorder(method, eigen.score, FUN=median), y=eigen.score)) +
  geom_boxplot(outlier.size = 0.5) + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  ylab("eigenscore") + xlab("method")
dev.off()


#visualize eigenscores
k=0
k=k+1
data.plot = data.frame(dim1=candidate.out[[1]][[k]][,1], dim2=candidate.out[[1]][[k]][,2], 
                       cluster=factor(info), eigenscore = c(ensemble.out[[2]][,k]), sindex = sindex[,k])
pdf(file = paste0("traj","-",candidate.out[[2]][k],"-eigenscore.pdf"), width=5,height=4) 
ggplot(data.plot, aes(x=dim1, y=dim2)) + geom_point(size=1.5, aes(shape=cluster, color=eigenscore)) +
  scale_shape_manual(values=1:nlevels(data.plot$cluster)) + ggtitle(title.name[k])
dev.off()
pdf(file = paste0("traj","-",candidate.out[[2]][k],"-eigenscore-box.pdf"), width=4,height=4) 
ggplot(data.plot, aes(x=cluster, y=eigenscore)) + geom_boxplot() + ggtitle(title.name[k])+ 
  ylab("eigenscore") + xlab("cluster")
dev.off()
pdf(file = paste0("traj","-",candidate.out[[2]][k],"-sindex-box.pdf"), width=4,height=4) 
ggplot(data.plot, aes(x=cluster, y=sindex)) + geom_boxplot() + ggtitle(title.name[k])+ 
  ylab("silhouette index") + xlab("cluster")#+ylim(-1,0.5)
dev.off()


##### boxplot of sindex
K.mat = exp(-as.matrix(ensemble.out[[1]])^2/quantile(ensemble.out[[1]],0.2)^2)
eigen.K.mat = eigen(K.mat)
ensemble.data = eigen.K.mat$vectors[,2:3] %*% diag(eigen.K.mat$values[2:3])
sindex.ensemble=silhouette(as.numeric(factor(info)),dmatrix = as.matrix(dist(ensemble.data)))[,3]
data.mat= as.matrix(dist(data))
sindex.data=silhouette(as.numeric(factor(info)),dmatrix = data.mat)[,3]
K.mat = exp(- apply(ensemble.out[[4]], c(2,3),mean)^2/quantile( apply(ensemble.out[[4]], c(2,3),mean),0.5)^2)
eigen.K.mat = eigen(K.mat)
naive.data = eigen.K.mat$vectors[,2:3] %*% diag(eigen.K.mat$values[2:3])
sindex.naive = silhouette(as.numeric(factor(info)),dmatrix = as.matrix(dist(naive.data)))[,3]
sindex.ori = silhouette(as.numeric(factor(info)),dmatrix = as.matrix(data.mat))[,3]
sindex.ext = cbind(sindex, sindex.ensemble, sindex.naive, sindex.ori)
colnames(sindex.ext)=c(candidate.out[[2]], "meta-spec", "meta-aver", "ori")
data.plot=data.frame(sindex = c(sindex.ext), method = rep(c(candidate.out[[2]],"meta-spec", "meta-aver", "ori"), each=n))
pdf(file = "traj-sindex-box.pdf", width=4,height=4) 
ggplot(data.plot, aes(x=reorder(method, sindex, FUN=median), y=sindex)) +
  geom_boxplot(outlier.size = 0.5) + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) + 
  ylab("silhouette index") + xlab("method")
dev.off()



############## meta-visualizations
data.plot = data.frame(dim1=ensemble.data[,1], dim2=ensemble.data[,2], cluster=factor(info))
pdf(file = paste0("traj-META-spec.pdf"), width=5,height=4) 
ggplot(data.plot, aes(x=dim1, y=dim2)) + geom_point(size=1.5, aes(shape=cluster, color=cluster))  + ggtitle("meta-spec")
dev.off()

data.plot = data.frame(dim1=naive.data[,1], dim2=naive.data[,2], cluster=factor(info))
pdf(file = paste0("traj-META-aver.pdf"), width=5,height=4) 
ggplot(data.plot, aes(x=dim1, y=dim2)) + geom_point(size=1.5, aes(shape=cluster, color=cluster))  + ggtitle("meta-aver")
dev.off()


############# kendall's tau

K=dim(ensemble.out[[4]])[1]
tau=c()
for(j in 1:K){
  data.traj=candidate.out[[1]][[j]]
  u=order(svd(data.traj)$u[,1])
  u2=order(-svd(data.traj)$u[,1])
  tau[j]=max(c(cor(x=as.numeric(data.names), y=as.numeric(data.names)[u], method = "kendall"), 
               cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2], method="kendall")))
  
}
data.traj=ensemble.data
u=order(svd(data.traj)$u[,1])
u2=order(-svd(data.traj)$u[,1])
tau[K+1]=max(c(cor(x=as.numeric(data.names), y=as.numeric(data.names)[u], method = "kendall"), 
               cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2], method="kendall")))

data.traj=naive.data
u=order(svd(data.traj)$u[,1])
u2=order(-svd(data.traj)$u[,1])
tau.naive=max(c(cor(x=as.numeric(data.names), y=as.numeric(data.names)[u], method = "kendall"), 
                   cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2], method="kendall")))


data.plot = data.frame(tau = c(tau, tau.naive), method = c(candidate.out[[2]], "meta-spec","meta-aver"))
pdf(file = "traj-tau.pdf", width=4,height=3) 
ggplot(data.plot, aes(reorder(method, tau), tau)) +
  geom_col() + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  ylab("kendall's tau") + xlab("method")
dev.off()


############## tau vs sindex
data.plot = data.frame(tau = c(tau[1:K], tau[K+1], tau.naive), 
                       sindex = c(apply(sindex,2,median), median(sindex.ensemble), median(sindex.naive)), 
                       method = c(candidate.out[[2]], "meta-spec", "meta-aver")
)

pdf(file = "traj-tau-sindex.pdf", width=4,height=4) 
ggplot(data.plot, aes(x=tau, y=sindex)) +
  geom_point(size=2) + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  ylab("median silhouette index") + xlab("kendall's tau")+geom_label_repel(aes(label = method),
                                                                            force=100,
                                                                            segment.color = 'grey50') 
dev.off()


