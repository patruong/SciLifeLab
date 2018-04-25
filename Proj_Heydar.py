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

# Specify dependent and independent variables
X = df.loc[:, df.columns[1]:]  # independent variables data
y = df.clonality  # dependent variable data