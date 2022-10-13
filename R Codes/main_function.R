###################################
#     Apr 28, 2022  by Rong Ma
###################################

######################################################
####### Obtain a list of candidate visualizations
######################################################

#load packages
install.packages(c("rARPACK","MASS","lle","dimRed","uwot","cluster","Rtsne","phateR"))
library(rARPACK)
library(MASS)
library(lle)
library(dimRed)
library(uwot)
library(cluster)
library(phateR)
library(Rtsne)


# candidate.visual: produces candidate visualizations based on 12 DR methods.
# kpca.sigma: bandwidth parameter for kPCA
# tsne.perplexity: perplexity for t-SNE
# phate.k: number of neighborhoods for PHATE
# umap.k: number of neighborhoods for UMAP
# cal_dist: logical, indicating whether distance matrix is needed. This can be set as FALSE if MDS, iMDS, Sammon are not used, to improve speed. 
# data: input n x p data set, with rows representing samples, and columns representing features.
# Output: a list of 2-dimensional embeddings for each candidate visualization, and their method names.

candidate.visual <- function(data, dim=2, methods= c("PCA", "MDS", "iMDS", "Sammon", "LLE", "HLLE","Isomap", 
                                                     "kPCA", "LEIM", "UMAP", "tSNE", "PHATE"), 
                             kpca.sigma = c(0.001, 0.002), 
                             umap.k= c(30, 50), 
                             tsne.perplexity = c(30, 50),
                             phate.k = c(30,50),
                             cal_dist = TRUE){
  n=dim(data)[1]
  
  dim.red.data = list()
  if(cal_dist){
    dist.data = dist(data)
  }
  
  name.method =c()
  #####################
  ####### PCA
  #####################
  
  
  i=0
  if(sum(methods == "PCA")>0){
    i=i+1
    pc.data = rARPACK::svds(as.matrix(data), k =dim)
    dim.red.data[[i]] = pc.data$u[,1:2]
    name.method =  c(name.method, "PCA")
  }
  #####################
  ####### classical MDS
  #####################
  
  if(sum(methods == "MDS")>0){
    i=i+1
    dim.red.data[[i]] = cmdscale(dist.data, k=dim)
    name.method =  c(name.method, "MDS")
  }
  #####################
  ####### isoMDS
  #####################
  
  if(sum(methods == "iMDS")>0){
    i=i+1
    imds.data = isoMDS(dist.data, k=dim)
    dim.red.data[[i]] = imds.data$points
    name.method =  c(name.method, "iMDS")
  }
  
  #####################
  ####### Sammon's nonlinear mapping
  #####################
  
  if(sum(methods == "Sammon")>0){
    i=i+1
    sam.data = sammon(dist.data, k=dim)
    dim.red.data[[i]] = sam.data$points
    name.method =  c(name.method, "Sammon")
  }
  
  #####################
  ####### LLE
  #####################
  
  if(sum(methods == "LLE")>0){
    i=i+1
    #k.sel = calc_k(data, m = 2)
    lle.data = lle(data, m=dim, k=20, reg=2)
    dim.red.data[[i]] =  lle.data$Y
    name.method =  c(name.method, "LLE")
  }
  
  #####################
  ####### HLLE
  #####################
  
  if(sum(methods == "HLLE")>0){
    i=i+1
    hlle.data <- embed(data, "HLLE", knn =20, ndim=dim)
    dim.red.data[[i]] = hlle.data@data@data
    name.method =  c(name.method, "HLLE")
  }
  
  #####################
  ####### isomap
  #####################
  
  if(sum(methods == "Isomap")>0){
    i=i+1
    imp.data <- embed(data, "Isomap", knn = 20, ndim=dim)
    dim.red.data[[i]] = imp.data@data@data
    name.method =  c(name.method, "Isomap")
  }
  
  #####################
  ####### kPCA
  #####################
  
  if(sum(methods == "kPCA")>0){
    for(j in 1:length(kpca.sigma)){
      i=i+1
      kpca.data <- embed(data,  "kPCA", kpar=list(sigma=kpca.sigma[j]), features=dim, ndim=dim)
      dim.red.data[[i]] = kpca.data@data@data
      name.method =  c(name.method, paste0("kPCA", j))
    }
  }
  
  #####################
  ####### Laplacian Eigenmap
  #####################
  
  if(sum(methods == "LEIM")>0){
    i=i+1
    lem.data <- embed(data,  "LaplacianEigenmaps", ndim=dim)
    dim.red.data[[i]] = lem.data@data@data
    name.method =  c(name.method, "LEIM")
  }
  
  #####################
  ####### UMAP
  #####################
  
  if(sum(methods == "UMAP")>0){
    for(j in 1:length(umap.k)){
      i=i+1
      umap.data <- umap(data,  n_neighbors = umap.k[j], n_components = dim)
      dim.red.data[[i]] = umap.data
      name.method =  c(name.method, paste0("UMAP", j))
    }
  }
  #####################
  ####### tSNE
  #####################
  
  if(sum(methods == "tSNE")>0){
    for(j in 1:length(tsne.perplexity)){
      i=i+1
      tsne.data <- embed(data,  "tSNE", perplexity = tsne.perplexity[j], ndim=dim)
      dim.red.data[[i]] = tsne.data@data@data
      name.method =  c(name.method, paste0("tSNE", j))
    }
  }
  
  #####################
  ####### PHATE
  #####################
  
  if(sum(methods == "PHATE")>0){
    for(j in 1:length(phate.k)){
      i=i+1
      dim.red.data[[i]] = phate(data, knn=phate.k[j], ndim=dim)$embedding
      name.method =  c(name.method, paste0("PHATE", j))
    }
  }
  
  
  
  
  
  return(list(dim.red.data, name.method))
}

#################################################
######## ensemble visualization function
#################################################

# ensemble.v.local: using spectral method for quantifying visualization quality and generating meta-visualization.
# data.list: a list of 2-dimensionoal embeddings, which is created by candidate.visual.
# name.method: the names of the candidate visualizations.
# original.data: an option to use the original data for the quality assessment, instead of using eigenscores.
# randomize: logical, indicating if a randomized procedure should be use to speed up the algorithm for dealing with large data.
# prop: numeric between 0 and 1, indicating proportion of samples used for the randomized procedure. Only used if randomize = TRUE.
# Output: a list containing (1) a meta-distance for meta-visualization, (2) eigenscores for candidate visualizations,
#         (3) names of the method ordered by averaged eigenscores, (4) distance matrices for candidate visualizations.

ensemble.v.local <- function(data.list, name.method=NA, original.data=NA){

  
  n=dim(data.list[[1]])[1]
  K=length(data.list)
  
  
  dist.list=array(dim=c(K,n,n))
  for(i in 1:K){
    dist.list[i,,]=as.matrix(dist(data.list[[i]]))
  }
  
  if(is.na(original.data)){
    ########## obtain weights
    
    ensemble.mat = matrix(ncol=n,nrow=n)
    weight = matrix(ncol=K,nrow=n)
    embed.mat.norm = list()
    for(j in 1:n){
      comp.mat = matrix(ncol=K, nrow=K)
      for(i in 1:K){
        embed.mat.norm[[i]] = dist.list[i,,j]/sqrt(sum(dist.list[i,,j]^2))
        for(k in 1:K){
          comp.mat[i,k] = sum(dist.list[k,,j]*dist.list[i,,j])/sqrt(sum(dist.list[k,,j]^2))/sqrt(sum(dist.list[i,,j]^2))
        }
      }
      weight[j,] = abs(eigen(comp.mat)$vectors[,1])
      
      ensemble.mat[,j] = apply(do.call(cbind, embed.mat.norm), 1, weighted.mean, w = weight[j,])*sum(weight[j,])
      if(j/1000==floor(j/1000)){
        print(j)
      }
    }
    
    
    return(list(ensemble.dist.mat=(ensemble.mat+t(ensemble.mat))/2, eigenscore=weight, 
                name.method=name.method[order(colMeans(weight), decreasing = T)], dist.list=dist.list))
    
  }else{
    
    ensemble.mat = matrix(ncol=n,nrow=n)
    weight = matrix(ncol=K,nrow=n)
    data.mat= as.matrix(dist(original.data))
    for(j in 1:n){
      
      eigen.score = c()
      embed.mat.norm = list()
      for(i in 1:K){
        eigen.score[i]=sum(data.mat[,j]*dist.list[i,,j])/sqrt(sum(data.mat[,j]^2))/sqrt(sum(dist.list[i,,j]^2))
        embed.mat.norm[[i]] = dist.list[i,,j]/sqrt(sum(dist.list[i,,j]^2)) 
      }
      
      eigen.score[which(eigen.score<0)]=0
      weight[j,] = eigen.score^2/sum(eigen.score^2)
      
      ensemble.mat[,j] = apply(do.call(cbind, embed.mat.norm), 1, weighted.mean, w = weight[j,])*sum(weight[j,])
      
    }
    
    return(list(ensemble.dist.mat=(ensemble.mat+t(ensemble.mat))/2, eigenscore=weight, 
                name.method=name.method[order(colMeans(weight), decreasing = T)], dist.list=dist.list))
    
  }

}

