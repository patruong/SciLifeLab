# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 15:15:29 2018

@author: Patrick
"""

from __future__ import print_function


import numpy as np
from matplotlib.pyplot import scatter
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import scipy as sp

from sklearn.preprocessing import scale
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from scipy import stats
from IPython.display import display, HTML


import statsmodels.api as sm

from os import listdir
from os.path import isfile, join
from os import chdir
from os import getcwd
from IPython.display import display, HTML


# To check distribution QQ-plot
from scipy.stats import probplot
from scipy import stats
import pylab

path = "C:\\Users\\Patrick\\Desktop\\Thesis In\\LukasK"

chdir(path)
print(getcwd())



df = pd.read_csv("P3128_LUKAS_PROTEIN_April2018.csv", sep = ",", index_col = 0, na_values = ["ND"])
df = pd.read_csv("P1902_LUKAS_PROTEIN_April2018.csv", sep = ",", index_col = 0, na_values = ["ND"])

# remove all bulk cells
df = df[df["Cell Number"] < 2] # remove all cell numbers larger than 1
df = df[df["Cell Number"] > 0] # remove all cell number smaller than 1

"Exclude certain columns"
# check 0:12 for both data set
# [0:12] is from cell tag : All Events FSC-H Mean
exclusion_list = [i for i in df.columns[0:11]]
exclusion_list.extend(["FACS_CHANNEL_CCR7",
                      "FACS_CHANNEL_DEXTRAMER",
                      "FACS_CHANNEL_LIN_NEG"])

# exclude RNA transcripts for Heydar stuff
for i in df.columns:
    if i[:4] == "RNA_":
        exclusion_list.append(i)

# Put back clonality indicator
exclusion_list.remove("clonality")

boolean_list_to_exclude_columns = ~df.columns.isin(exclusion_list)


df = df[df.columns[boolean_list_to_exclude_columns]]

# Rename columns for convenience
orig_columns = df.columns
new_columns = ["V"+str(i) for i in range(1, len(df.columns)+1)]
columns_map = pd.DataFrame(data = zip(new_columns, orig_columns), columns = ["New", "Original"])

df.columns = new_columns  # rename column names to be similar to R naming convention
df.V1 = df.V1.astype(str)

# Specify dependent and independent variables
X = df.loc[:, df.columns[1]:]  # independent variables data
y = df.V1 # dependent variable data


#############
# ANALYSIS ##
#############


## MATRIX SCATTERPLOT

pd.tools.plotting.scatter_matrix(df.loc[:, df.columns[1]:], diagonal="kde")
plt.tight_layout()
plt.show()

pd.tools.plotting.scatter_matrix(df.loc[:, df.columns[1]:], diagonal="hist")
plt.tight_layout()
plt.show()

sns.lmplot("V4", "V5", df, hue="V1", fit_reg=False);

ax = df[["V2","V3","V4","V5","V6"]].plot()
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5));

#expressions of the different genes for different protein markers
df.loc[:, df.columns[1]:].transpose().plot.bar()



df.mean()
df.std()

def printMeanAndSdByGroup(variables, groupvariable):
    data_groupby = variables.groupby(groupvariable)
    print("## Means:")
    display(data_groupby.apply(np.mean))
    print("\n## Standard deviations:")
    display(data_groupby.apply(np.std))
    print("\n## Sample sizes:")
    display(pd.DataFrame(data_groupby.apply(len)))

printMeanAndSdByGroup(X, y)

# Clonality mean, std and sizes
X.groupby(y).mean()
X.groupby(y).std()
X.groupby(y).apply(len)

def calcWithinGroupsVariance(variable, groupvariable): #own formula for understanding
    # find out how many values the group variable can take
    unique_groups = groupvariable.unique()
    num_unique_groups = len(unique_groups)
    # get the mean and standard deviation for each group
    num_total = 0
    denom_total = 0
    for i in unique_groups:
        group_i = variable[groupvariable==i]
        len_group_i = len(group_i)
        # get the standard deviation for group i
        sd_i = np.std(group_i)
        # Within-group variance formula
        num_i = (len_group_i)*sd_i**2
        denom_i = (len_group_i)
        # Summation procedure in within-group variance formula
        num_total = num_total + num_i
        denom_total = denom_total + denom_i
    V_w = num_total / (denom_total - num_unique_groups)
    return V_w

def calcBetweenGroupsVariance(variable, groupvariable): #own formula for understanding
    
    # find out how many values the group variable can take
    unique_groups = groupvariable.unique()
    num_unique_groups = len(unique_groups)
    
    #calculate the overall grand mean
    grand_mean = np.mean(variable)
    
    # the the mean and standard deviation for each group
    num_total = 0
    denom_total = 0
    for i in unique_groups:
        group_i = variable[groupvariable == i]
        len_group_i = len(group_i)
        # get the mean and standard deviation for group i
        mean_i = np.mean(group_i)
        std_i = np.std(group_i)
        # Between-group variance formula
        num_i = (len_group_i) * ((mean_i - grand_mean)**2)
        denom_i = (len_group_i)
        # Summation procedure in between-group variance formula
        num_total = num_total + num_i
        denom_total = denom_total + denom_i
    # valvulate the between-groups variance
    V_b = num_total / (num_unique_groups - 1) 
    return(V_b)

# Calculate Seperation acieved by all variables in a multivariate dataset
def calcSeparations(variables, groupvariable):
    # Calculate the separation for each variable
    for i in variables:
        variable_i = variables[i]
        V_w = calcWithinGroupsVariance(variable_i, groupvariable) 
        V_b = calcBetweenGroupsVariance(variable_i, groupvariable)
        sep = V_b/V_w
        print( "variable", i, "V_w=", V_w, "V_b=", V_b, "separation=", sep)



def calcWithinGroupsCovariance(variable1, variable2, groupvariable):
    # find out how many values the group variables can take
    unique_groups = groupvariable.unique()
    num_unique_groups = len(unique_groups)
    cov_w = 0.0
    
    # get the covarance of variable 1 and variable 2 for each groups
    for i in unique_groups:
        group_i_var1 = variable1[groupvariable == i]
        group_i_var2 = variable2[groupvariable == i]
        mean_var1 = np.mean(group_i_var1) #for each group in var1
        mean_var2 = np.mean(group_i_var2) #for each group in var2
        len_group_i = len(group_i_var1)
        # get the covariance for this group
        cov_j = 0.0
        for q,k in zip(group_i_var1, group_i_var2):
            cov_j += (q - mean_var1)*(k - mean_var2)
        cov_group_i = cov_j 
        cov_w += cov_group_i
    totallength = len(variable1)
    cov_w = cov_w / (totallength - num_unique_groups)
    return cov_w

def calcBetweenGroupsCovariance(variable1, variable2, groupvariable):
    # find out how many values the goup variable can take
    unique_groups = groupvariable.unique()
    num_unique_groups = len(unique_groups)
    # calculate the grand means
    var1_Gmean = np.mean(variable1)
    var2_Gmean = np.mean(variable2)
    # calculate the between-groups covariance
    cov_b = 0.0
    for i in unique_groups:
        group_i_var1 = variable1[groupvariable == i]
        group_i_var2 = variable2[groupvariable == i]
        mean_var1 = np.mean(group_i_var1)
        mean_var2 = np.mean(group_i_var2)
        len_group_i = len(group_i_var1)
        cov_i = (mean_var1 - var1_Gmean) * (mean_var2 - var2_Gmean) * len_group_i
        cov_b += cov_i
    cov_b = cov_b / (num_unique_groups - 1)
    return cov_b

"""
var_wg = calcWithinGroupsVariance(X,y)
var_bg = calcBetweenGroupsVariance(X,y)
cov_wg = calcWithinGroupsCovariance(X,y)
cov_bg = calcBetweenGroupsCovariance(X,y)
"""


#############
## MANOVA ##
############

IV = X.columns 

var_wg = pd.DataFrame(columns = IV, index = ["var_wg"])
for i in IV:
    var_wg_i = calcWithinGroupsVariance(X[i], y)
    var_wg[i] = var_wg_i

var_bg = pd.DataFrame(columns = IV, index = ["var_bg"])
for i in IV:
    var_bg_i = calcBetweenGroupsVariance(X[i], y)
    var_bg[i] = var_bg_i

cov_wg = pd.DataFrame(columns = IV, index = IV) # eye is var_wg
for i in IV:
    for j in IV:
        cov_wg_ij = calcWithinGroupsCovariance(X[i],X[j],y)
        cov_wg[i][j] = cov_wg_ij

cov_bg = pd.DataFrame(columns = IV, index = IV) # eye is var_bg
for i in IV:
    for j in IV:
        cov_bg_ij = calcBetweenGroupsCovariance(X[i],X[j],y)
        cov_bg[i][j] = cov_bg_ij
        
        
var_wg["V2"]["var_wg"]+var_wg["V3"]["var_wg"]+var_bg["V2"]["var_bg"]+var_bg["V3"]["var_bg"]+cov_bg["V2"]["V3"]+cov_wg["V2"]["V3"]

"check for multinomial goodness of fit"
"https://www.biostat.wisc.edu/~kbroman/teaching/labstat/fourth/notes02.pdf

#######################################
## MULTIVARIATE GAUSSIAN DISTRIBUTION #
#######################################

df_classes = df.V1
df_unique = df_classes.value_counts()

# Permute clonality and log FACS values
def log_df(df, celltag_col = "V1"):
    multiIndex = zip(df.index, df[celltag_col])
    index = pd.MultiIndex.from_tuples(multiIndex)
    index.names = ["cellTag", celltag_col]
    df_m = pd.DataFrame(df.drop([celltag_col], axis = 1).values, index = index, columns = df.columns[1:])
    
    df_adjuster = df_m.min()
    # adjust so that minimal value for each columns i 1 by adding the min + 1
    df_newValues = df_m - df_adjuster + 1 
    df_log = np.log(df_newValues)
    return df_log
    
# Permute clonality and log FACS values
def permute_log_df(df, celltag_col = "V1"):
    multiIndex = zip(df.index, np.random.permutation(df[celltag_col]))
    index = pd.MultiIndex.from_tuples(multiIndex)
    index.names = ["cellTag", celltag_col]
    df_m = pd.DataFrame(df.drop([celltag_col], axis = 1).values, index = index, columns = df.columns[1:])
    
    df_adjuster = df_m.min()
    # adjust so that minimal value for each columns i 1 by adding the min + 1
    df_newValues = df_m - df_adjuster + 1 
    df_log = np.log(df_newValues)
    return df_log





###################################################################
# DO I REALLY NEED TO SIMULATE THE VALUES TO COMPARE DIFFERNECE? ##
###################################################################
    
# Function for return Multivariate normal distribution dataframe
def MVN_pdf(df, clonality, clone_col = "V1", sample_size = 1000):
    """
    Input should be multiIndex df from log_df or permute_log_df
    
    clonality => clonality strain e.g. "E4"
    sample_size = samples to draw from MVN
    """
    
    #df_log = log_df(df, clone_col)
    df = df[df.index.get_level_values(clone_col) == clonality]
    #means = pd.DataFrame(np.mean(df_log), columns = None)
    #cov = pd.DataFrame(np.cov(df_log.transpose()), index = df_log.columns, columns = df_log.columns)
    means = df.mean()
    cov = df.cov()
    samples = sample_size
    multiVar_norm = np.random.multivariate_normal(means, cov, samples)
    df_MVN = pd.DataFrame(multiVar_norm)
    return df_MVN

clonality = "E4"
clone_col = "V1"
sample_size = 1000

"""
Input should be multiIndex df form log_df or permute_log_df

clonality ==> clonality strain e.g. "E4"
sample_size = samples to draw from MVN

Output

"""
df_log = log_df(df, celltag_col = "V1")
    df_mvn = df_log[df_log.index.get_level_values(clone_col) == clonality]
    df_mvn_complement = df_log[df_log.index.get_level_values(clone_col) != clonality]

    #means = pd.DataFrame(np.mean(df_log), columns = None)
    #cov = pd.DataFrame(np.cov(df_log.transpose()), index = df_log.columns, columns = df_log.columns)
    means = df_mvn.mean()
    cov = df_mvn.cov()
    
    #calculate same mena + cov for complement
    samples = sample_size
    multiVar_norm = np.random.multivariate_normal(means, cov, samples)
    df_MVN = pd.DataFrame(multiVar_norm)

"""
Challenge is to find a test for multivariate distribution equality
"""
    
df_log = log_df(df)
df_mvn = df_log[df_log.index.get_level_values(clone_col) == clonality]
df_mvn_complement = df_log[df_log.index.get_level_values(clone_col) != clonality]

X = df_mvn.transpose()
X = X.T
Y = df_mvn_complement.transpose()
Y = Y.T

E_orig,T_orig = energyStatistic(X,Y)


e_stat_list = []
iterations = 10
 
for i in range(iterations):
    df_log = permute_log_df(df)
    df_mvn = df_log[df_log.index.get_level_values(clone_col) == clonality]
    df_mvn_complement = df_log[df_log.index.get_level_values(clone_col) != clonality]
    X = df_mvn#.transpose()
    Y = df_mvn_complement#.transpose()
    E,T = energyStatistic(X,Y)
    e_stat_list.append(E)
    
#E_df = pd.DataFrame(E_orig)
E_df = pd.DataFrame([E_orig for i in range(len(e_stat_list))])
E_perm_df = pd.DataFrame(e_stat_list)

boolean_table = (E_df-E_perm_df).transpose()>0
1-boolean_table.transpose().sum()/boolean_table.transpose().count()
# Similarity between two matrices
#https://math.stackexchange.com/questions/507742/distance-similarity-between-two-matrices

def absolute_diff(X,Y):
    absolute_diff = abs(X-Y)
    res = absolute_diff.sum().sum()
    return res

def ss_diff(X,Y):
    ss = np.power((X-Y),2)
    ss_sqrt = np.sqrt(ss)
    res = ss_sqrt.sum().sum()
    return res

def maximalValue_diff(X,Y):
    abs_diff = abs(X-Y)
    res = abs_diff.max().max()
    return res


#### TAKE A SUBSET
df_log = log_df(df, "V1")

#random vs random case
def boot_random(clone, iterations = 100, sample_size = 1000, clone_col = "V1", standardize = True):
    """
    clone => clonality e.g. "E4"
    """
    diffs = pd.Series(index = range(iterations))
    for i in range(iterations):
        df_log = permute_log_df(df, clone_col)
        df_perm = permute_log_df(df, clone_col)
        
        cGroup_MVN = MVN_pdf(df_log, clone, clone_col = clone_col, sample_size = sample_size)
        pGroup_MVN = MVN_pdf(df_perm, clone, clone_col = clone_col, sample_size = sample_size)
        
        diff = ss_diff(cGroup_MVN, pGroup_MVN)
        diffs[i] = diff
    if standardize == True:
        diffs = (diffs - diffs.mean())/diffs.std()
    return diffs


#random vs random case
def boot_clone(clone, iterations = 100, sample_size = 1000, clone_col = "V1", standardize = True):
    """
    clone => clonality e.g. "E4"
    """
    df_log = log_df(df, clone_col)
    diffs = pd.Series(index = range(iterations))
    for i in range(iterations):

        df_perm = permute_log_df(df, clone_col)
        
        cGroup_MVN = MVN_pdf(df_log, clone, clone_col = clone_col, sample_size = sample_size)
        pGroup_MVN = MVN_pdf(df_perm, clone, clone_col = clone_col, sample_size = sample_size)
        
        diff = ss_diff(cGroup_MVN, pGroup_MVN)
        diffs[i] = diff
    if standardize == True:
        diffs = (diffs - diffs.mean())/diffs.std()
    return diffs


baseline = boot_random("E4")
permutation_test = boot_clone("E4")


clone_col = "V1"
clone = "E4"
sample_size = 100

df_log = permute_log_df(df, clone_col)
df_perm = permute_log_df(df, clone_col)
        
cGroup_MVN = MVN_pdf(df_log, clone, clone_col = clone_col, sample_size = sample_size)
pGroup_MVN = MVN_pdf(df_perm, clone, clone_col = clone_col, sample_size = sample_size)
        
cGroup_MVN = cGroup_MVN.transpose()
pGroup_MVN = pGroup_MVN.transpose()

# Norm functions
def absoluteNorm(X):
    """
    Takes a vector and calculates the absolute norm
    
    if X = U - V, then we have the absolute distance between U and V.
    """"
    absNorm = abs(X)
    return absNorm

def euclideanNorm(X):
    """
    Takes a vector and calculates the Euclidean norm
    
    if X = U - V, then euclideanNorm is the Euclidean distance between U and V
    """
    
    squaredVec = np.square(X)
    sumVec = np.sum(squaredVec)
    eucNorm = np.sqrt(sumVec)
    return eucNorm

def mahalanobisNorm(X, Y):
    """
    Takes one vector and calculates the Mahalanobis distance of it OR
    takes two random vectors of the same distribution and calculates the 
    Mahalanobis dissimilarity between them.
    
    Note vectors should be in matrix notation np.array[[]], because 
    we need to transpose multiplicate.
    
    X,Y - should be verticle vector / verticle multi-dim vector.
    
    Notice: numpy transpose operations are reversed
    """
    if X.ndim != 2:
        error_dim = X.ndim
        raise ValueError("Vector ndim must be 2. Vector ndim is "+str(error_dim)+
                             ". Check with np.array().ndim. X vector Error.")
    elif Y.ndim != 2:
        error_dim = Y.dim
        raise ValueError("Vector ndim must be 2. Vector ndim is "+str(error_dim)+
                             ". Check with np.array().ndim. Y vector Error.")
        
    e = (X-Y)
    v_stacked = np.vstack([X,Y])
    cov = np.cov(v_stacked.transpose())
    inv_cov = np.linalg.inv(cov)
    D = np.sqrt(np.sum(np.dot(e,inv_cov)*e, axis = 1))
    
    e = (Y-X)
    v_stacked = np.vstack([Y,X])
    cov = np.cov(v_stacked.transpose())
    inv_cov = np.linalg.inv(cov)
    D = np.sqrt(np.sum(np.dot(e,inv_cov)*e, axis = 1))
    D= np.dot(e,inv_cov)*e
    
    #http://mccormickml.com/2014/07/21/mahalanobis-distance/
    
    e = (X-Y)
    mu = np.mean(e, axis = 1)
    cov = np.dot((e.T - mu).T,(e.T - mu))/(e.T.shape[0]-1)
    inv_cov = np.linalg.inv(cov)
    D = np.sum(np.dot(e.T, inv_cov)*e.T, axis = 1)
    
    #Mahalanobis tutorial
    #http://people.revoledu.com/kardi/tutorial/Similarity/MahalanobisDistance.html    
    
    return D

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



# energy statistic - testing for equal distributions wikipedia
'''
def energyStatistic(X,Y, distance_function = euclideanNorm):
    """    
    Input
        X,Y - two vectors or matrices
        
            columns - x,y values
            rows - samples
        
        distance_function - norm function to use
        
        
    n = len(cGroup_MVN.columns)
    m = len(pGroup_MVN.columns)
    
    A = 0
    for i in cGroup_MVN:
        for j in pGroup_MVN:
            temp = abs(cGroup_MVN[i] - pGroup_MVN[j])
            A += temp
    A = A / (n*m)
            
    
    B = 0 
    for i in cGroup_MVN:
        for j in cGroup_MVN:
            temp = abs(cGroup_MVN[i] - cGroup_MVN[j])
            B += temp
    B = B / (n*n)
    
    C = 0
    for i in pGroup_MVN:
        for j in pGroup_MVN:
            temp = abs(pGroup_MVN[i] - pGroup_MVN[j])
            C += temp
    C = C / (m*m)
    
    E = 2*A - B - C #if F is zero ==> same distribution of p and c
    T = ((n*m)/(n+m))*F
    """
    X = X.transpose()
    Y = Y.transpose()
    
    n = len(X.columns)
    m = len(Y.columns)
    
    #n = len(X.index)
    #m = len(Y.index)
    A = 0
    
    for i in X:
        for j in Y:
            diff = X[i] - Y[j]
            dist = distance_function(diff)
            A += dist
    A = A / (n*m)
            
    
    B = 0 
    for i in X:
        for j in X:
            diff = X[i] - X[j]
            dist = distance_function(diff)
            B += dist
    B = B / (n*n)
    
    C = 0
    for i in Y:
        for j in Y:
            diff = Y[i] - Y[j]
            dist = distance_function(diff)
            C += dist
    C = C / (m*m)
    
    E = 2*A - B - C #if F is zero ==> same distribution of p and c
    T = ((n*m)/(n+m))*E
    return E, T
'''

def energyStatistic(X,Y, distance_function = mahalanobisDistance):
    """    
    Input
        X,Y - two vectors or matrices
        
            columns - x,y values
            rows - samples
        
        distance_function - norm function to use
        
    """
    X = X.transpose()
    Y = Y.transpose()
    
    print("We want index to be features V2...V11")
    print("Printing index")
    print(X.index)
    
    n = len(X.columns)
    m = len(Y.columns)
    
    #n = len(X.index)
    #m = len(Y.index)
    A = 0
    
    for i in X:
        for j in Y:
            #diff = X[i] - Y[j]
            X_i = np.array([X[i].values]).T
            Y_j = np.array([Y[j].values]).T
            dist = distance_function(X_i, Y_j)
            A += dist
    A = A / (n*m)
            
    
    B = 0 
    for i in X:
        for j in X:
            #diff = X[i] - X[j]
            X_i = np.array([X[i].values]).T
            X_j = np.array([X[j].values]).T
            dist = distance_function(X_i, X_j)
            B += dist
    B = B / (n*n)
    
    C = 0
    for i in Y:
        for j in Y:
            #diff = Y[i] - Y[j]
            X_i = np.array([Y[i].values]).T
            Y_j = np.array([Y[j].values]).T
            dist = distance_function(Y[i], Y[j])
            C += dist
    C = C / (m*m)
    
    E = 2*A - B - C #if F is zero ==> same distribution of p and c
    T = ((n*m)/(n+m))*E
    return E, T



clone_col = "V1"
clones = df_unique.index
sample_size = 100

iterations = 3
E_stats = pd.DataFrame(columns = clones, index = df.columns[1:])
T_stats = pd.DataFrame(columns = clones, index = df.columns[1:])
for j in range(iterations):
    for i in clones:
    
        df_log = permute_log_df(df, clone_col)
        df_perm = permute_log_df(df, clone_col)
                
        cGroup_MVN = MVN_pdf(df_log, clone, clone_col = clone_col, sample_size = sample_size)
        pGroup_MVN = MVN_pdf(df_perm, clone, clone_col = clone_col, sample_size = sample_size)
                
        cGroup_MVN = cGroup_MVN.transpose()
        pGroup_MVN = pGroup_MVN.transpose()
        
        E, T = energyStatistic(cGroup_MVN, pGroup_MVN)
        E.index = df.columns[1:]
        T.index = df.columns[1:]
        
        E_stats[i] = E
        T_stats[i] = T
    E_stats = E_stats + E_stats
    T_stats = T_stats + T_stats
E_stats = E_stats / iterations
T_stats = T_stats / iterations    




"""
G "E4"
not G

draw many samples
G "Permute"
not G

Compare how many

G "E4" not G 
is lower than 
G "Permute" not G

Look up more at Energy Distance?

Look up P-value for permutation test
Look up Permutation Statistics


ToDO:
    
    Implement
    - Bhattacharyya Distance
    - Hellinger Distance
    - Kullback-Leibler Divergence
"""


############################
## BHATTACHARYYA DISTANCE ##
############################

# https://stats.stackexchange.com/questions/296361/intuition-of-the-bhattacharya-coefficient-and-the-bhattacharya-distance
# https://en.wikipedia.org/wiki/Bhattacharyya_distance
# Bhattacharyya distance generalizes the Mahalanobis distance
# Mahalanobis ditsance ==> same covariance matrix

X = np.random.multivariate_normal([0,0], [[1,0],[0,1]], 1000)
X2 = np.random.multivariate_normal([0,0], [[1,0],[0,1]], 100)
Y = np.random.multivariate_normal([2,5], [[5,2],[2,7]], 15)

#P-value for Bhattacharyya and Hellinger distance = 1/max(sample_size_1, sample_size_2)

def BhattacharyyaDistance(X,Y, axis = 0, weighted_cov = False):
    """
    input:
            columns - features
            rows - samples
    """
    if axis == 1:
        X = X.T
        Y = Y.T
        
    X_mu = np.mean(X, axis = 0)
    Y_mu = np.mean(Y, axis = 0)
    
    X_cov = np.cov(X.T, bias = True)
    Y_cov = np.cov(Y.T, bias = True)
    
    #weighted_cov = False
    if weighted_cov == True:
        X_w = (float(X.shape[0])/(float(X.shape[0])+float(Y.shape[0])))
        Y_w = (float(Y.shape[0])/(float(X.shape[0])+float(Y.shape[0])))
        pooled_cov = (X_w*X_cov)+(Y_w*Y_cov)
    else:
        pooled_cov = (X_cov + Y_cov)/2
    
    inv_cov = np.linalg.inv(pooled_cov)
    
    pooled_det = np.linalg.det(pooled_cov)
    X_det = np.linalg.det(X_cov)
    Y_det = np.linalg.det(Y_cov)
    
    comp1 = (float(1)/float(8))*np.dot(np.dot(np.array([(X_mu - Y_mu)]), inv_cov), 
             np.array([(X_mu - Y_mu)]).T)
    comp2 = (float(1)/float(2))*np.log(pooled_det / np.sqrt(np.dot(X_det, Y_det)))
    
    bhattacharyya = comp1 + comp2
    
    return bhattacharyya
##### END BHATTACHARYYA #####

########################
## HELLINGER DISTANCE ##
########################
    
X = np.random.multivariate_normal([0,0], [[1,0],[0,1]], 10)
Y = np.random.multivariate_normal([2,5], [[5,2],[2,7]], 15)

def HellingerDistance(X,Y, axis = 0, weighted_cov = False):
    """
    input:
            columns - features
            rows - samples
    """
    if axis == 1:
        X = X.T
        Y = Y.T
    
    X_mu = np.mean(X, axis = 0)
    Y_mu = np.mean(Y, axis = 0)
    
    X_cov = np.cov(X.T, bias = True)
    Y_cov = np.cov(Y.T, bias = True)
    
    #weighted_cov = False
    if weighted_cov == True:
        X_w = (float(X.shape[0])/(float(X.shape[0])+float(Y.shape[0])))
        Y_w = (float(Y.shape[0])/(float(X.shape[0])+float(Y.shape[0])))
        pooled_cov = (X_w*X_cov)+(Y_w*Y_cov)
    else:
        pooled_cov = (X_cov + Y_cov)/2
    
    inv_cov = np.linalg.inv(pooled_cov)
    
    comp1 = np.divide(np.dot(np.power(np.linalg.det(X_cov),float(1)/float(4)), 
                             np.power(np.linalg.det(Y_cov),float(1)/float(4))),
            np.power(np.linalg.det(pooled_cov), float(1)/float(2)))
            
    comp2 = np.exp(-(float(1)/float(8))*np.dot(np.dot(np.array([(X_mu - Y_mu)]), inv_cov), 
                     np.array([(X_mu - Y_mu)]).T))  
    
    hellinger = 1 - comp1*comp2 

    return hellinger

#### END HELLINGER #####



#######################

diff = ss_diff(cGroup_MVN, pGroup_MVN)
diffs[i] = diff



"""TEST CASES
d1 = boot("E4")
d2 = boot("F3")
d3 = boot("A7")
    
d1.plot()
d2.plot()
d3.plot()  
"""

"""
df_log = log_df(df, "V1")
df_1 = df_log[df_log.index.get_level_values("V1") == "E4"]
means = df_1.mean()
cov = df_1.cov()
samples = 1000
multiVar_norm = np.random.multivariate_normal(means, cov, samples)
df_MVN_1 = pd.DataFrame(multiVar_norm)

df_perm = permute_log_df(df, "V1")
df_2 = df_perm[df_perm.index.get_level_values("V1") == "E4"]
means = df_2.mean()
cov = df_2.cov()
samples = 1000
multiVar_norm = np.random.multivariate_normal(means, cov, samples)
df_MVN_2 = pd.DataFrame(multiVar_norm)
"""

ss_diff(df_MVN_1, df_MVN_2)
"""
df_log[df_log.index.get_level_values("V1") == "F3"]
#means = pd.DataFrame(np.mean(df_log), columns = None)
#cov = pd.DataFrame(np.cov(df_log.transpose()), index = df_log.columns, columns = df_log.columns)
means = df_log.mean()
cov = df_log.cov()
samples = 1000
multiVar_norm = np.random.multivariate_normal(means, cov, samples)
df_MVN_2 = pd.DataFrame(multiVar_norm)
"""




#fit multivariate normal distirbution to these parameters
# perhaps enought to check difference between these parameters???


#df_m[df_m.index.get_level_values("V1") == "F3"]



#############################
## MULTINOMIAL DISTRIBUTION #
#############################

df_classes = df.V1
df_unique = df_classes.value_counts()
df_probs = df_unique/df_unique.sum()

#sampling procedure
draws = np.random.multinomial(n=100, pvals=df_probs) #draw form each group
df_draws = pd.Series(data = draws, index = df_unique.index)

df_samples = pd.DataFrame(columns = df.columns)
for i in range(len(df_draws)):
    df_temp = df[df.V1 == df_draws.index[i]].sample(df_draws[i])
    df_samples = df_samples.append(df_temp)
    
#check the distribution of gene expression

#check distribution of V2 - RNA
df.V2.mean()
df.V2.std()

#log the data

df.drop(["V1"], axis = 1).min()
#find all min smaller than 0

df_adjuster = ((df.drop(["V1"], axis = 1).min()<0)*1)*df.drop(["V1"], axis = 1).min()

# adjust so that minimal value for each columns i 1 by adding the min + 1
df_newValues = (df.drop(["V1"], axis = 1) - df_adjuster +1 )

df_logValues = np.log(df_newValues) #negative rexpressions?

df_log = df_logValues.join(df.V1)
#df_log = df_logValues[new_columns] # rearrange columns to new_columns


#NOW DATA IS NORMALLY DISTRIBUTED!!



### 
#df[df.columns != "V2"]

QQ = probplot(df_log.V3, dist=stats.norm, plot=pylab)
QQ = probplot(df.V3, dist=stats.norm, plot=pylab)

def empirical_dist(series, bins = 100, title = "Empirical and Fitted probability distribution functions",
                   xlabel = "FACS values"):
    mean = series.mean()
    std = np.sqrt(series.var())

    plt.figure(figsize = (8,6)) 
    h = sorted(series)
    h = pd.DataFrame(h)
    fit = stats.norm.pdf(h, mean, std)
    
    plt_label = "N(%.2f, %.2f)" % (mean, std)
    plt.plot(h,fit, '-', label = plt_label)
    plt.legend()
    plt.hist(h, normed = True, bins = bins, label = "Empirical dist")#, plt.legend()
    plt.title(title)   
    plt.ylabel("Probability")
    plt.xlabel(xlabel)


