# -*- coding: utf-8 -*-
"""
Created on Sun May 13 16:27:26 2018

@author: Patrick
"""

# Mahalanobis Distance Tutorial
# http://people.revoledu.com/kardi/tutorial/Similarity/MahalanobisDistance.html

import numpy as np
 
def generateMahalanobis_testData():
    """
    Test data is same as in link above.
    """
    g1 = np.array([[2,2,6,7,4,6,5,4,2,1],
                  [2,5,5,3,7,4,3,6,5,3]])
    g1 = g1.T
    g2 = np.array([[6,7,8,5,5],
                  [5,4,7,6,4]])
    g2 = g2.T
    return g1, g2

def mahalanobisDistance(g1, g2, axis = 0):
    """
    Computes the Mahalanobis distance between two vectors.
    
    input format should be:
        g1 - vector with features in column and samples in row
        g2 == g1
        
    """
    if axis == 1:
        g1 = g1.T
        g2 = g2.T
    g1_mu = np.mean(g1, axis = 0)
    g2_mu = np.mean(g2, axis = 0)
    
    g1_centered = g1-g1_mu
    g2_centered = g2-g2_mu
    
    g1_cov = np.cov(g1_centered.T, bias = True)
    g2_cov = np.cov(g2_centered.T, bias = True)
    
    g1_w = (float(g1.shape[0])/(float(g1.shape[0])+float(g2.shape[0])))
    g2_w = (float(g2.shape[0])/(float(g1.shape[0])+float(g2.shape[0])))
    
    pooled_cov = (g1_w*g1_cov)+(g2_w*g2_cov)
    try:
        inv_cov = np.linalg.inv(pooled_cov)
    except:
        inv_cov = 1/pooled_cov
    mean_diff = np.array([g1_mu-g2_mu])
    
    mahalanobis = np.sqrt(np.dot(np.dot(mean_diff, inv_cov), mean_diff.T))
    return mahalanobis

if __name__ == "__main__":
    g1, g2 = generateMahalanobis_testData()
    D = mahalanobisDistance(g1, g2)
    print(D)
    #should be 1.41
    D = mahalanobisDistance(g1.T, g2.T, axis = 1)
    print(D)
    #should be 1.41 also.
    
"Brick In R for MCMC"
"PyMC3 for MCMC"
