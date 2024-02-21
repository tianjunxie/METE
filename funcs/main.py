#!/usr/bin/env python

### Tianjun Xie
### NCAT 2023-2024

### Main function that serves dual training and predicting purposes


### Usage
### Supports mix and match
### >> import main
### >> main.predict(['CC(=O)C[Pt]', 'CCC[Ir]']) or main.predict(['CC(C(C(C)([Pt])O[Pt])=O)=O','C([Ru])C(C([Ru])([Ru]))([Ru])CC'])
###




### Dependencies and their functions
### freq_mat: feature matrix generator
### ga1: molecule fragementator and linear evaluator @ ga lvl1
### ga2: molecule fragementator and nonlinear evaluator @ ga lvl2
### GPR: Guassian Process Regression class

### Main data loadable
### json files, discretized group data @ ga lvl1
### pkl files, GPR models @ lvl2

from rdkit import Chem
import numpy as np
import pandas as pd
import re
import sys

import freq_mat
import ga1
import ga2
# import ui
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw

import warnings
warnings.filterwarnings('ignore')
### no line wrapping
np.set_printoptions(threshold=sys.maxsize)


def que(self, groups, property_set_name):
    """
    group query (for training)
    
    Verify groups present.
    """
    missing_groups = [group for group in groups \
                          if property_set_name not in self[group]]
    if missing_groups:
        # f_report = open('./report.txt', 'w')
        # f_report.write(missing_groups+"\n")
        # f_report.flush()
        print(missing_groups)
        return(missing_groups)


### GA defragmentation
### choice of GA1 or GA2
### and properties calculator
def calc(name, lvl):
    """
    Given a smiles string and ga lvl number, analyze the smiles and calculate properties
    ***To improve

    Parameters
    ----------
    name : str
        RDKit compatible SMILES string
    lvl : int
        level of ga fragmentation, 1 or 2

    Returns
    -------
    str(cleaner[0]) : str
        sanitized SMILES string that is RDKit compatible
    cp_new : list
        list of thermochemisty
    dict(descriptors) : dict
        dictionary of group names, and occurances

    """    
    name = name.replace('[{Ru}]','[Ru]')
    name = name.replace('[{Pt}]','[Pt]')
    if name == '[Pt]' or name == '[Ru]':
        return(0)
    elif (name == 'H([Pt])') or (name == 'H[Pt]'):
        return(-10.587)
    elif (name == 'H([Ru])') or (name == 'H[Ru]'):
        return(-18.450008104)
    elif (name == 'C(=O)(=O)') or (name == 'O=C=O'):
        return(1)
    else:
        ### new GA routines
        ### 1st or 2nd order GA, based on input 'lvl'
        cps = np.zeros(17)
        descriptors = ga2.frag(name, lvl)
        ### smiles clean-up
        cleaner = re.split('{|}',str(descriptors))
        cleaner.pop()     ### pop the last element
        cleaner.pop(0)    ### addtional text cleanup
        # print(name, dict(descriptors))
        return( str(cleaner[0]), cps,  dict(descriptors))    ### testing for GA verbal extraction for now

def composition(sml):
    """
    Given a smiles string, analyze the elemental composition

    Parameters
    ----------
    sml : str
        RDKit compatible SMILES string

    Returns
    -------
    atomic_count.items : list
        composition numbers by order

    """

    from collections import defaultdict
    mol =Chem.MolFromSmiles(sml)
    mh = Chem.AddHs(mol)
    atomic_count = defaultdict(lambda: 0)
    for atom in mh.GetAtoms():
        atomic_count[atom.GetAtomicNum()] += 1
    s = sorted(atomic_count.items())
    print(s)


def train(molecs, verbose=False):
    ### main driver ##
    ### modes
    ### 1: train, takes a spreadsheet and yields training matricies
    ### 2: predict, provide a ML model and predict
    ### verbose: controls verbosity of the fragmentation information and others

    ### training mode
    text = []
    # df = pd.read_csv('trainingPt.csv')
    df = pd.read_csv('trainingIr.csv')
    # df = pd.read_csv('trainingRu.csv')
    smiles_strings = df['SMILES'].to_numpy()
    for i in smiles_strings:
        # i = i.replace('{M}','Pt')
        if (i == 'H([Pt])') or (i == 'H[Pt]'):
            continue
        elif (i == 'H([Ru])') or (i == 'H[Ru]'):
            continue
        ### sanitize the SMILES so it complies with pgradd
        m = Chem.MolFromSmiles(i.split()[0])
        a = Chem.MolToSmiles(m)
        # print(a)

        # # ### get compostions
        # composition(a)

    #     ### GA lvl1
    #     smiles,ener, groups = calc(a, 1)
    #     text.append(smiles)
    # print('GA1 Fragmentation')

        ### GA lvl2
        smiles,ener, groups= calc(a, 2)
        text.append(smiles)
    print('------GA2 Fragmentation------')

    ### verbose switch
    verbose = True
    ### get the frequency mat
    freq_mat.get_freq_mat(text,verbose)
        
        
def predict(molecs, verbose=False):
### predicting mode, default
    print("OUTPUT: Hf298, S298, Cp[100-1500K]\n")
    # ### Prioritize keyboard input
    # if len(sys.argv) > 1:
    #     molecs = sys.argv[1:]
    # else:
    #     # ## Taking the last molecs from the training set
    #     molecs = [
    #         # 'C([Pt])(O([Pt]))C(=O)C(O)C([Pt])(=O)',
    #         # '[Pt]OC([Pt])C(O)C(C([Pt])([Pt])O)=O',
    #         # 'CC(C(C(C)([Pt])O[Pt])=O)=O',
    #         'CC([Pt])CO',
    #     ]
    #     # molecs = [
    #     #     'C([Ru])C([Ru])(C([Ru])([Ru])([Ru]))CC',
    #     #     'C([Ru])C([Ru])(C([Ru]))C([Ru])([Ru])C',
    #     #     'C([Ru])C(C([Ru])([Ru])([Ru]))C([Ru])C',
    #     # ]

    ### evaluate at ga lvl1
    nasa9, raw = ga1.eval(molecs, verbose)
    if verbose:
        for arr in raw:
            print(np.round(arr, 2))
    # print(np.array(nasa9))

    ### evaluate at ga lvl2
    groups_2,y_out = ga2.ml_eval(molecs, verbose)
    if verbose:
        for arr in y_out:
            print(np.round(arr, 2))
    # [print(x) for x in groups_2]
    # ### need to concatenate the ouput
    # if not np.all(y_out == 0):
    #     y_out = np.concatenate(y_out,axis=0)
    #     print(y_out)

    ### Adding final calibration onto lvl1
    # y_correct = np.array(raw) +np.array(y_out)
    # Perform element-wise addition between corresponding arrays
    y_correct = [arr1 + arr2 for arr1, arr2 in zip(raw, y_out)]
    # print(y_correct)
    for ind, ele in enumerate(molecs):
        print(ele, *y_correct[ind].round(2))
    return(y_correct)

    # ms = Chem.MolFromSmiles(molecs)
    # oin = Draw.MolToImage(ms) ### this returns a PIL object
    # oin = Draw.ShowMol(ms)
    ### oin = Draw.MolToMPL(ms)
    ### Draw.MolToImageFile(ms, "_smiles.png")

    # oin.show()





