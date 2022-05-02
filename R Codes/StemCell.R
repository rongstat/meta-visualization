#######################################
#     Apr 29, 2022  by Rong Ma
#   Mouse Embryonic Stem Cells (Cycle)
#######################################

library(cluster)
library(ggplot2)
library(uwot)
library(ggrepel)
library(Rtsne)
library(phateR)
library(Seurat)

###################################
######## Load mnist data 
###################################
setwd("~/Download/Data/E-MTAB-2805")#set directory of the example datasets.
#read data
g2m=read.table("G2M_singlecells_counts.txt", header =T)
s=read.table("S_singlecells_counts.txt", header =T)
g1=read.table("G1_singlecells_counts.txt", header =T)
dim(g2m)
dim(s)
dim(g1)

gene = g1[,1:4]
g1=g1[,-(1:4)]
s=s[,-(1:4)]
g2m=g2m[,-(1:4)]
combined = cbind(g1, s, g2m)

#remove ERCCs
combined = combined[-which(substr(gene[,1], 1,4)=="ERCC"),]
gene= gene[-which(substr(gene[,1], 1,4)=="ERCC"),]
dim(combined)
#[1] 38298   288

data=combined
colnames(data)

info = c(rep(c("G1","S","G2M"), each=96))
table(info)
#normalization and QC
data <- CreateSeuratObject(counts = data, project = "cell_cycle", min.cells = 3, min.features = 20)
data <- NormalizeData(data)

#identify highly variable genes
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)

#get the normalized data matrix
scale.data = data@assays[["RNA"]]@scale.data
scale.data.var = scale.data[match(data@assays[["RNA"]]@var.features,rownames(scale.data)),]
dim(scale.data.var)
#[1] 2000  288
data.names = factor(info)
levels(data.names)=c(1,3,2)
data=t(scale.data.var) 
n=dim(data)[1]

#####select relevant genes for cycling
cor.gene=c()
for(i in 1:2000){
  cor.gene[i]=min(c(t.test(data[which(info=="G1"),i], data[-which(info=="G1"),i])$p.value,
                    t.test(data[which(info=="G2M"),i], data[-which(info=="G2M"),i])$p.value,
                    t.test(data[which(info=="S"),i], data[-which(info=="S"),i])$p.value))
}
sum(p.adjust(cor.gene, method="BH")<0.3)
data=data[,which(p.adjust(cor.gene, method="BH")<0.3)]
dim(data)

################### ensemble visualization

#generate 16 candidatet visualizations
set.seed(30)
candidate.out = candidate.visual(data, method=c("PCA", "MDS", "iMDS", "Sammon", "LLE", "HLLE", "Isomap", 
                                          "kPCA", "LEIM", "UMAP", "tSNE","PHATE"),
                                 kpca.sigma = c(0.01, 0.001), 
                                 umap.k= c(30, 50), 
                                 tsne.perplexity = c(10, 50),
                                 phate.k = c(30,50))

#obtain eigenscores and meta-distance matrix
ensemble.out = ensemble.v.local(data.list=candidate.out[[1]], name.method = candidate.out[[2]])

setwd("~/meta-visualization/output/cycle") #set directory for the output
save(candidate.out, ensemble.out, file="cycle.RData")



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
pdf(file = paste0("cycle","-",candidate.out[[2]][k],".pdf"), width=5,height=4) 
ggplot(data.plot, aes(x=dim1, y=dim2)) + geom_point(size=1.5, aes(shape=cluster, color=cluster))  + ggtitle(title.name[k])
dev.off()



#obtain silhouette indices for candidate viz
K=dim(ensemble.out[[4]])[1]
sindex = matrix(ncol=K,nrow=n)
for(i in 1:K){
  sindex[,i]=silhouette(as.numeric(factor(info)), dmatrix = ensemble.out[[4]][i,,])[,3]
}

#boxplot for eigenscores
data.plot = data.frame(eigen.score = c(ensemble.out[[2]]), method = rep(candidate.out[[2]], each=n))
pdf(file = "cycle-score-box.pdf", width=4,height=4) 
ggplot(data.plot, aes(x=reorder(method, eigen.score, FUN=median), y=eigen.score)) +
  geom_boxplot(outlier.size = 0.5) + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  ylab("eigenscore") + xlab("method")
dev.off()


#visualize eigenscores
k=0
k=k+1
data.plot = data.frame(dim1=candidate.out[[1]][[k]][,1], dim2=candidate.out[[1]][[k]][,2], 
                       cluster=factor(info), eigenscore = c(ensemble.out[[2]][,k]), sindex = sindex[,k])
pdf(file = paste0("cycle","-",candidate.out[[2]][k],"-eigenscore.pdf"), width=5,height=4) 
ggplot(data.plot, aes(x=dim1, y=dim2)) + geom_point(size=1.5, aes(shape=cluster, color=eigenscore)) +
  scale_shape_manual(values=1:nlevels(data.plot$cluster)) + ggtitle(title.name[k])
dev.off()
pdf(file = paste0("cycle","-",candidate.out[[2]][k],"-eigenscore-box.pdf"), width=4,height=4) 
ggplot(data.plot, aes(x=cluster, y=eigenscore)) + geom_boxplot() + ggtitle(title.name[k])+ 
  ylab("eigenscore") + xlab("cluster")+ylim(0.2,0.3)
dev.off()
pdf(file = paste0("cycle","-",candidate.out[[2]][k],"-sindex-box.pdf"), width=4,height=4) 
ggplot(data.plot, aes(x=cluster, y=sindex)) + geom_boxplot() + ggtitle(title.name[k])+ 
  ylab("silhouette index") + xlab("cluster")#+ylim(-1,0.5)
dev.off()



##### boxplot of sindex
K.mat = exp(-as.matrix(ensemble.out[[1]])^2/quantile(ensemble.out[[1]],0.5)^2)
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
pdf(file = "cycle-sindex-box.pdf", width=4,height=4) 
ggplot(data.plot, aes(x=reorder(method, sindex, FUN=median), y=sindex)) +
  geom_boxplot(outlier.size = 0.5) + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) + 
  ylab("silhouette index") + xlab("method")
dev.off()



############## meta-visualizations
data.plot = data.frame(dim1=ensemble.data[,1], dim2=ensemble.data[,2], cluster=factor(info))
pdf(file = paste0("cycle-META-spec.pdf"), width=5,height=4) 
ggplot(data.plot, aes(x=dim1, y=dim2)) + geom_point(size=1.5, aes(shape=cluster, color=cluster))  + ggtitle("meta-spec")
dev.off()

data.plot = data.frame(dim1=naive.data[,1], dim2=naive.data[,2], cluster=factor(info))
pdf(file = paste0("cycle-META-aver.pdf"), width=5,height=4) 
ggplot(data.plot, aes(x=dim1, y=dim2)) + geom_point(size=1.5, aes(shape=cluster, color=cluster))  + ggtitle("meta-aver")
dev.off()


############# calculate kendall's tau
K=dim(ensemble.out[[4]])[1]

shifter <- function(x, n = 1) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}
tau = c()
for(k in 1:K){
  u1=candidate.out[[1]][[k]][,1]
  u2=candidate.out[[1]][[k]][,2]
  u=order(atan(u1/u2),decreasing =F)
  ure = order(atan(u1/u2),decreasing =T)
  u2= match(1:n,u)
  u2re = match(1:n,ure)
  tau[k]= max(c(cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2], method="kendall"),
                      cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2re], method="kendall")))
  for(i in 1:(n-1)){
    u2= match(shifter(1:n, n =i),u)
    u2re = match(shifter(1:n, n =2),ure)
    temp = max(c(cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2], method="kendall"),
                 cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2re], method="kendall")))
    if(temp >  tau[k]){
      tau[k]=temp
    }
  }
}

#kendall's tau for meta-spec
u1=ensemble.data[,1]
u2=ensemble.data[,2]
u=order(atan(u1/u2),decreasing =F)
ure = order(atan(u1/u2),decreasing =T)
u2= match(1:n,u)
u2re = match(1:n,ure)
tau[K+1]= max(c(cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2], method="kendall"),
              cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2re], method="kendall")))
for(i in 1:(n-1)){
  u2= match(shifter(1:n, n =i),u)
  u2re = match(shifter(1:n, n =2),ure)
  temp = max(c(cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2], method="kendall"),
               cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2re], method="kendall")))
  if(temp >  tau[K+1]){
    tau[K+1]=temp
  }
}
#kendall's tau for meta-aver
u1=naive.data[,1]
u2=naive.data[,2]
u=order(atan(u1/u2),decreasing =F)
ure = order(atan(u1/u2),decreasing =T)
u2= match(1:n,u)
u2re = match(1:n,ure)
tau.naive= max(c(cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2], method="kendall"),
                 cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2re], method="kendall")))
for(i in 1:(n-1)){
  u2= match(shifter(1:n, n =i),u)
  u2re = match(shifter(1:n, n =2),ure)
  temp = max(c(cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2], method="kendall"),
               cor(x=as.numeric(data.names), y=as.numeric(data.names)[u2re], method="kendall")))
  if(temp >  tau.naive){
    tau.naive=temp
  }
}

data.plot = data.frame(tau = c(tau, tau.naive), method = c(candidate.out[[2]], "meta-spec","meta-aver"))

pdf(file = "cycle-tau.pdf", width=4,height=4) 
ggplot(data.plot, aes(reorder(method, tau), tau)) +
  geom_col() + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  ylab("kendall's tau") + xlab("method")
dev.off()


##### circle plot
k=0
k=k+1
data.circ=data.frame(candidate.out[[1]][[k]][,1], candidate.out[[1]][[k]][,2])
data.circ=t(scale(t(data.circ), center = F))
col = factor(info)
levels(col)=c("red","forestgreen", "cornflowerblue")
pdf(file = paste0("cycle","-cir-",candidate.out[[2]][k],".pdf"), width=4,height=4) 
plot(data.circ, col=as.character(col),  xlab="dim 1", ylab="dim 2", font.lab=2, cex.lab = 1.2,
     main = title.name[k],cex.main =1, 
     pch=c(17, 15, 19)[as.numeric(factor(info))]) 
dev.off()


data.circ=data.frame(ensemble.data[,1], ensemble.data[,2])
data.circ=t(scale(t(data.circ), center = F))
col = factor(info)
levels(col)=c("red","forestgreen", "cornflowerblue")
pdf(file = paste0("cycle-cir-meta.pdf"), width=4,height=4) 

plot(data.circ, col=as.character(col),  xlab="dim 1", ylab="dim 2", font.lab=2, cex.lab = 1.2,
     main = "meta-spec",cex.main =1, 
    pch=c(17, 15, 19)[as.numeric(factor(info))]) 
dev.off()


############## tau vs sindex


data.plot = data.frame(tau = c(tau[1:K], tau[K+1], tau.naive), 
                       sindex = c(apply(sindex,2,median), median(sindex.ensemble), median(sindex.naive)), 
                       method = c(candidate.out[[2]], "meta-spec", "meta-aver")
)

pdf(file = "cycle-tau-sindex.pdf", width=4,height=4) 
ggplot(data.plot, aes(x=tau, y=sindex)) +
  geom_point(size=2) + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  ylab("median silhouette index") + xlab("kendall's tau")+geom_label_repel(aes(label = method),
                                                                            force=100,
                                                                            segment.color = 'grey50') 
dev.off()
 
