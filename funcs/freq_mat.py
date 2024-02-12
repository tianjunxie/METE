#!/usr/bin/env python
from collections import Counter
from sklearn.ensemble import RandomForestClassifier
import numpy as np
import re
import sys
import pprint
np.set_printoptions(linewidth=150)

def findndel(myList, verbose):
    ### Find all empty elements
    empty_elements = [i for i, x in enumerate(myList) if not x]

    # ### turned off for deployment
    # ### Print their indices
    # if verbose:
    #     print(empty_elements)

    ### Delete them
    myList = [x for x in myList if x != '']
    
    return(myList)

def add2list(a,myList):
    # a = a.replace(' ','')
    if a not in myList:
        myList.append(a)

def search_col(arr,fet):
    ### get rank of the matrix
    rank0 = np.linalg.matrix_rank(arr)
    ## Loop over the columns
    for i in range(arr.shape[1]):
        # Delete the current column
        new_arr = np.delete(arr, i, axis=1)
        # Check the rank of the new array
        rank1 = np.linalg.matrix_rank(new_arr)
        ### potential uneccessary feature
        if rank1 == rank0:
        ### and rare appearance:
            if (arr[:,i].sum()) <=2:
                print(f'Rank unaffected after deleting feature {fet[i]}')

def find_linearly_dependent_columns(arr):
    # Compute the rank of the matrix formed by the columns
    rank = np.linalg.matrix_rank(arr)

    # Check if the rank is less than the number of columns
    if rank < arr.shape[1]:
        # Find the linearly dependent columns
        dependent_columns = []
        for i in range(arr.shape[1]):
            subset = np.delete(arr, i, axis=1)
            subset_rank = np.linalg.matrix_rank(subset)
            if subset_rank == rank:
                dependent_columns.append((i, i+1))
        # print("Linearly dependent column pairs:")
        # for pair in dependent_columns:
        #     print(pair)
        return dependent_columns
    else:
        return None


def rank_most_occurring_columns(pairs, num_cols, top_n):
    counter = Counter()

    for pair in pairs:
        counter.update(pair)

    top_columns = [column for column, count in counter.most_common(top_n) if column < num_cols]
    return top_columns


def rank_features(matrix):
    num_samples, num_features = matrix.shape
    
    # Create a dummy target variable
    target = np.zeros(num_samples)
    
    # Create a Random Forest Classifier model
    model = RandomForestClassifier()
    
    # Fit the model
    model.fit(matrix, target)
    
    # Get feature importances
    feature_importances = model.feature_importances_
    
    # Sort features based on importance scores
    sorted_indices = np.argsort(feature_importances)[::-1]
    
    # Return the sorted indices
    return sorted_indices

def get_freq_mat(text, verbose):
    ### pre-processing
    ### rid of any empty entries- no 2nd NN molecs
    text = findndel(text, verbose)
    # for i in text:
    #     i = i.replace('\'','').replace('\'','').replace(' ','')
    #     i = i.splitlines()
    #     line = i.split('\,')
    #     print(line) 
    
    ### elem construct
    col = []
    for i in text:
        x = i.replace('\'','').replace('\'','').replace(' ','')
        line = re.split('\,',x)
        # line = re.split('\,',i)
        for j in line:
            ele = re.split('\:',j)
    ## ele[0] will be used for building the columns for GROUPS ARRAY
            add2list(ele[0],col)

    ## initialize the matrix with 0s
    M = len(col)
    N = len(text)
    res = np.zeros((N,M))
    linn = 0
    ## repeat again for ele[1] this time for frequency
    ## and ele[1] will be used to fill the frequency rows
    for i in text:
        # line = re.split('\,',i)
        x = i.replace('\'','').replace('\'','').replace(' ','')
        line = re.split('\,',x)
        for j in line:
            ele = re.split('\:',j)
            for k in range(len(col)):
                if ele[0] == col[k]:
                    res[linn,k] = ele[1]
        linn = linn + 1
    num_rows, num_cols = res.shape    

    # # ### find linearly dependent columns
    # search_col(res, col)
    
    # # dependent_columns = find_linearly_dependent_columns(res)
    # # print("Linearly dependent columns found:", dependent_columns)
    # # top_n = 10
    # # top_columns = rank_most_occurring_columns(dependent_columns, res.shape[1], top_n)
    # # print("Top", top_n, "most occurring columns:")
    # # print(top_columns)

    # sorted_indices = rank_features(res)
    
    # print("Ranked features (from most important to least):")
    # for i in sorted_indices:
    #     print(col[i])
    # print(sorted_indices)
    
    # ### turned off for deployment
    # if verbose:
    #     print("shape: ", num_rows, num_cols)
    #     print("rank: ", np.linalg.matsrix_rank(res))
    #     print(col)
    #     # print(res)
    #     pprint.pprint(res)
    return(col,res)




