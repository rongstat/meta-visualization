'''
Starter code for Python implementation of meta-visualization

Eric Sun (2022)
'''

from sklearn.metrics import pairwise_distances
from sklearn.neighbors import NearestNeighbors
from scipy.linalg import eigh
from sklearn.manifold import TSNE, MDS
import numpy as np

    
def meta_viz(visualizations, projection_method=None):
    '''
    Inputs:
       visualizations [list of numpy arrays]
        - list of numpy arrays of size n x 2 (columns being the two-dimensional coordinates of visualization)
        - each row should be matched to the same observation/sample across all elements in the list
       projection_method [str or None]
        - options are "tSNE" or "MDS" for reprojecting the meta-distance matrix into a meta-visualization
    
    Approach:
        1. Computes euclidean distance matrices from 2D embeddings
        2. Ensembles the euclidean distance according to eigenvector with largest eigenvalue
        3. Uses projection_method ["MDS", "tSNE"] to transform euclidean distance matrix to 2D embeddings
        
    Returns:
        meta_visualization [numpy array] - size nx2 two-dimensional visualization genereated using the meta-visualization approach
            - returned only if projection_method is not None
        meta_distances [numpy array] - size nxn pairwise euclidean distance matrix generated using meta-visualization approach
    '''
    # define number of samples
    n = visualizations[0].shape[0]
    K = len(visualizations)
    
    # Iterate and record distance matrices
    X_distance_matrix_list = []
    
    # compute pairwise distance matrix
    for X_embedded in visualizations:
        assert X_embedded.shape[0] == n, "All visualization arrays need to have the same number of rows"
        assert X_embedded.shape[1] == 2, "All visualization arrays need to have two columns"
        X_distance = pairwise_distances(X_embedded)
        X_distance_matrix_list.append(X_distance)
             
    # Compute weights for meta-visualization
    meta_distances = np.zeros((n,n))
    weights = np.zeros((n,K))
    for j in range(n):
        # fill in comparison matrix
        comp_mat = np.zeros((K,K))
        for i in range(K):
            for k in range(K):
                comp_mat[i,k] = np.sum(X_distance_matrix_list[k][:,j]*X_distance_matrix_list[i][:,j])/np.sqrt(np.sum(X_distance_matrix_list[k][:,j]**2))/np.sqrt(np.sum(X_distance_matrix_list[i][:,j]**2))
        # Eigenscore
        w, v = np.linalg.eig(comp_mat)
        weights[j,:] = np.abs(v[:,0])
    
        # Ensembles distance matrices
        matrix_norms = []
        for i in range(K):
            matrix_norms.append(X_distance_matrix_list[i][:,j]/np.sqrt(np.sum(X_distance_matrix_list[i][:,j]**2)))
        
        temp = np.zeros(matrix_norms[0].shape)
        for i in range(K):
            temp += matrix_norms[i]*weights[j,i]
            
        meta_distances[:,j] = temp
    
    meta_distances = np.nan_to_num((meta_distances+meta_distances.T)/2)
    
    # Re-project on 2D embedding space
    if projection_method is None:
        return(meta_distances)
    else:
        if projection_method == "tSNE":
            tsne = TSNE(n_components=2, metric="precomputed", verbose=0).fit(meta_distances)
            meta_visualization = tsne.embedding_
        elif projection_method == "MDS":
            mds = skm.MDS(n_components=2, verbose=0, dissimilarity="precomputed").fit(meta_distances)
            meta_visualization = mds.embedding_
    
        return(meta_visualization, meta_distances)