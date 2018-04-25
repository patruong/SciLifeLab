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
    return Covb

