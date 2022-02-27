## Distance metrics using Numba lib in python to speed up computation because R is slow
## and I couldnt figure out how to use DoParrallel libarary


import pandas as pd
import numpy as np
import numpy.ma as ma
import numpy.linalg as lin
import matplotlib.pyplot as plt
import numba as nb
from numba import prange
import scipy.spatial.distance as sc


  



## Weighted Cosine Sim
@nb.njit
def cosine_similarity_weighted(vec1, vec2, mask1, mask2, weight1, weight2):
    mask = mask1 & mask2
    if len(mask) > 0:
        vec1_ma = vec1[mask]
        vec2_ma = vec2[mask]
        vec1_w = weight1[mask]/np.sum(weight1[mask])
        vec2_w = weight2[mask]/np.sum(weight2[mask])
        vec1_ma = np.multiply(vec1_ma, vec1_w)
        vec2_ma = np.multiply(vec2_ma, vec2_w)
        numerator = np.dot(vec1_ma, vec2_ma)
        vec1_norm = np.sqrt(np.dot(vec1_ma, vec1_ma))
        vec2_norm = np.sqrt(np.dot(vec2_ma, vec2_ma))
        dist = numerator / (vec1_norm * vec2_norm + 10e-10)
        return dist
    return 0
@nb.njit(parallel = True)
def sparse_cosine_similarity_weight(matrix, mask, weight):
    assert matrix.shape == mask.shape
    #assert matrix.shape == weight.shape
    dim1, dim2 = matrix.shape
    result = np.zeros((dim1,dim1))
    for i in prange(dim1):
        for j in prange(dim1):
            result[i,j] = cosine_similarity_weighted(matrix[i],matrix[j],mask[i],mask[j],weight[i],weight[j])
    return result



## Cosine Sim

@nb.njit
def cosine_similarity(vec1, vec2, mask1, mask2):
    mask = mask1 & mask2
    if len(mask) > 0:
        vec1_ma = vec1[mask]
        vec2_ma = vec2[mask]
        numerator = np.dot(vec1_ma, vec2_ma)
        vec1_norm = np.sqrt(np.dot(vec1_ma, vec1_ma))
        vec2_norm = np.sqrt(np.dot(vec2_ma, vec2_ma))
        dist = numerator / (vec1_norm * vec2_norm + 10e-10)
        return dist
    return 0

@nb.njit(parallel = True)
def sparse_cosine_similarity(matrix, mask):
    assert matrix.shape == mask.shape
    dim1, dim2 = matrix.shape
    result = np.zeros((dim1,dim1))
    for i in prange(dim1):
        for j in prange(dim1):
            result[i,j] = cosine_similarity(matrix[i],matrix[j],mask[i],mask[j])
    return result


## Dot Product

@nb.njit
def dot_similarity(vec1, vec2, mask1, mask2):
    mask = mask1 & mask2
    vec1_ma = vec1[mask]
    if len(vec1_ma) > 0:
        vec2_ma = vec2[mask]
        norm_dot = np.dot(vec1_ma, vec2_ma)/len(vec1_ma)
        #vec1_norm = np.sqrt(np.dot(vec1_ma, vec1_ma))
        #vec2_norm = np.sqrt(np.dot(vec2_ma, vec2_ma))
        #dist = numerator / (vec1_norm * vec2_norm + 10e-10)
        return norm_dot
    return 0
  
@nb.njit(parallel = True)
def sparse_dot(matrix, mask):
    assert matrix.shape == mask.shape
    dim1, dim2 = matrix.shape
    result = np.zeros((dim1,dim1))
    for i in prange(dim1):
        for j in prange(dim1):
            result[i,j] = dot_similarity(matrix[i],matrix[j],mask[i],mask[j])
    return result

## Euclidian Distance divided by vector length 

#@nb.njit
def euclid_similarity(vec1, vec2, mask1, mask2):
    mask = mask1 & mask2
    vec1_ma = vec1[mask]
    if len(vec1_ma) > 0:
        vec2_ma = vec2[mask]
        vec1_ma = np.array(vec1_ma)
        vec2_ma = np.array(vec2_ma)
        spat_dist = np.linalg.norm(vec1_ma-vec2_ma)/len(vec1_ma)
        return spat_dist
    return 0


#@nb.njit(parallel = True)
def sparse_euclid(matrix, mask):
    assert matrix.shape == mask.shape
    dim1, dim2 = matrix.shape
    result = np.zeros((dim1,dim1))
    for i in range(dim1):
        for j in prange(dim1):
            result[i,j] = euclid_similarity(matrix[i],matrix[j],mask[i],mask[j])
    return result







