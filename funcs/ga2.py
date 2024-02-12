# -*- coding: utf-8 -*-
### script for group additivity
### fragmentation up to 2nd nearest neighbors


from collections import defaultdict
from rdkit import Chem
import numpy as np
import os
import re
import freq_mat
import GPR

def Get_metal(text):
    if 'Pt' in text:
        return 'Pt'
    if 'Ru' in text:
        return 'Ru'
    if 'Ir' in text:
        return 'Ir'

def atom_sani(listin):
    """
    Sanitizes the atom list and output

    Parameters
    ----------
    listin : str
        input string of words

    Returns
    -------
        sanitized_list : str
        sanitzed string for RDKit

    """
    common_names = ['C','H','O','Pt','Ru','Ir','=C','=O']
    counts = []
    sanitized_list = []
    for i in common_names:
        counts.append(listin.count(i))
    for idx, ele in enumerate(counts):
        if ele == 1:
            sanitized_list.append('('+common_names[idx]+')')
        if ele > 1:
            sanitized_list.append('('+common_names[idx]+')'+str(counts[idx]))
    # for i in sanitized_list:
    #     print(i,sep='',end='')
    # print('\n')
    return(''.join(sanitized_list))


def get_neighbor(i,mol):
    """
    Neighbor atoms finder, 1st nearest

    Parameters
    ----------
    i : int
        atom object in a RDKit molecule
    mol : obj
        RDKit molecule obj

    Returns
    -------
        [atom_neigh, idx_neigh] : [str, int]
        neigboring atom symbol and the index of the neighbors

    """
    idx_neigh = []
    atom_list = []
    BO={Chem.BondType.SINGLE:1,
        Chem.BondType.DOUBLE:2,
        Chem.BondType.TRIPLE:3,
        Chem.BondType.AROMATIC:1.5,
        Chem.BondType.UNSPECIFIED:0
        }
    # iatom = i.GetSymbol()
    natom = i.GetAtomicNum()
    # if (iatom =='C') or \
    #     (iatom =='O'):
    ### read in center heteratom
    ### and only for light atoms
    if natom <= 16:
        # print(mol.GetAtomWithIdx(i).GetSymbol(),end='')
        neigh = i.GetNeighbors()
        ### loop over neighbors
        for j in neigh:
            ### Getting bond order GetBondBetweenAtoms
            # print(i.GetIdx(),j.GetIdx())
            bo = BO[mol.GetBondBetweenAtoms(i.GetIdx(),j.GetIdx()).GetBondType()]
            if bo == 1:
                jatom = j.GetSymbol()
            elif bo == 2:
                jatom = '=' + j.GetSymbol()
            elif bo == 3:
                jatom = '#' + j.GetSymbol()
            atom_list.append(jatom)
            atom_neigh = atom_sani(atom_list)
            # print('(',jatom,')', sep='')
            if ('C' in jatom) or \
               ('O' in jatom):
                idx_neigh.append(j.GetIdx())
                # print([atom_neigh, idx_neigh])
        return([atom_neigh, idx_neigh])


### func to defrag the molec based on
### the SMILES string (str) and level (int)
def frag(smiles, lvl):
    """
    Fragment a molecule based on the choice of ga lvls

    Parameters
    ----------
    smiles : str
        sanitzied, RDKit allowed SMILES string
    lvl : int
        GA levels, choose between 1 or 2.

    Returns
    -------
    dict(counts) : dict
        dictionary of group names, and occurances

    """
    mol=Chem.MolFromSmiles(smiles)
    ### add H back
    mol = Chem.AddHs(mol)
    metals = ['Pt','Ru','Ir']
    # subgrph = Chem.rdmolops.FindAllSubgraphsOfLengthN(mol, 3)
    # for row in subgrph:
    #   for item in row:
    #       print(item, end=' ')
    #   print()
    
    # Get all bonds in the molecule
    # bonds = [(x.GetBeginAtomIdx(), x.GetEndAtomIdx()) for x in mol.GetBonds()]
    
    # l_2 = []
    # r_2 = []
    # for i in range(len(mol.GetAtoms())):
    # ### mark as central atom
    # ### Going Left
    #     for idx1 in bonds:
    #         if idx1[1] == i:
    #             l_1_atom = idx1[0]
    #     for idx2 in bonds:
    #         if idx2[1] == l_1_atom:
    #             l_2.append([idx2[0],l_1_atom,i])
    # ### Going Right
    #     for idx3 in bonds:
    #         if idx3[0] == i:
    #             r_1_atom = idx3[1]
    #     for idx4 in bonds:
    #         if idx4[0] == r_1_atom:
    #             r_2.append([i,r_1_atom,idx4[1]])
    
    # mol.GetAtomWithIdx(0).GetSymbol()
    # [0].GetSymbol()
    
    atom_neigh_list = []
    atom_neigh2_list = []
    idx_neigh_list= []
    idx_neigh2_list = []

    ### new loop
    for i in mol.GetAtoms():
        ### loop for C hearatoms
        hearatom = i.GetSymbol()
        if (hearatom == 'C'):
            ### first nearest neighbor atoms and indices
            atom_neigh, idx_neigh = get_neighbor(i, mol)
            ### following convention: hearatom + neighbor atoms
            atom_neigh_list.append(hearatom + atom_neigh)
            idx_neigh_list.append(idx_neigh)
            # print(i.GetIdx(), idx_neigh,atom_neigh_list)
            # for j in idx_neigh:
            #     atom_neigh2, idx_neigh2=get_neighbor(mol.GetAtomWithIdx(j), mol)
            #     # print(f"{hearatom},{atom_neigh} >> {mol.GetAtomWithIdx(j).GetSymbol()},{atom_neigh2}")
            #     atom_neigh2_list.append(mol.GetAtomWithIdx(j).GetSymbol() + atom_neigh2)
            #     idx_neigh2_list.append(idx_neigh2)
                
        ### loop for O hearatoms
        elif (hearatom == 'O'):
            ### first neighbor atom and indices
            atom_neigh, idx_neigh=get_neighbor(i, mol)
            ### following pgradd convention: hearatom + neighbor atoms
            if ('C' not in atom_neigh):
                atom_neigh_list.append(hearatom + atom_neigh)
                idx_neigh_list.append(idx_neigh)
   #################################################################################
    ### second nearest neighbor atom and indices
    for nnlist2 in idx_neigh_list:
        if len(nnlist2) >= 2: ### that has a 2nd NN
            # print(nnlist2)
            # Outer loop to iterate over each element in the array j
            for i in range(len(nnlist2)):
                # Inner loop to iterate over the remaining elements
                for j in range(len(nnlist2)):
                    # Skip the element being 'removed'
                    if j == i:
                        continue
                    # get 2NN neighbors
                    # print("Remaining element:", nnlist2[j])
                    atom_neigh2, idx_neigh2=get_neighbor(mol.GetAtomWithIdx(nnlist2[j]), mol)
                    # print(f"{hearatom},{atom_neigh} >> {mol.GetAtomWithIdx(j).GetSymbol()},{atom_neigh2}")
                    nn2out = mol.GetAtomWithIdx(nnlist2[j]).GetSymbol() + atom_neigh2
                    if not nn2out.startswith("O"):
                        atom_neigh2_list.append(nn2out)
                        idx_neigh2_list.append(idx_neigh2)
                    ### else if it is O based groups
                    else:
                        if not any(item in nn2out for item in metals):
                            atom_neigh2_list.append(nn2out)
                            idx_neigh2_list.append(idx_neigh2)    
    
    # use dict comprehension to count the number of times each element occurs in the list
    counts = {x: atom_neigh_list.count(x) for x in atom_neigh_list}
    ### now converts to defaultdict objects
    # d_lvl_1 = defaultdict(int, **counts)
    
    d_lvl_1 = dict(counts)
    counts = {x: atom_neigh2_list.count(x) for x in atom_neigh2_list}
    ### now converts to defaultdict objects
    # d_lvl_2 = defaultdict(int, **counts)
    d_lvl_2 = dict(counts)
    if lvl == 1:
        return(d_lvl_1)
    elif lvl == 2:
        return(d_lvl_2)
    # print(d_lvl_2)
    # print('\n')
    

def ml_eval(molecs, verbose):
    """
    ML(GPR) evaluator of the R matrix, at ga lvl2

    Parameters
    ----------
    molecs : list str
        list of smiles strings of molecules of interest.
    verbose : boolean
        switch for verbose ouput.

    Returns
    -------
    group_out, y_out : list str, list float
        groups names at ga lvl2, R matrix regression results 

    """
    
    ### 18 permittable features for the Pt database
    ga2_pt_features=np.array(
    ['C(C)(Pt)3', 'O(C)(H)', 'O(=C)', 'C(C)(H)2(Pt)', 'C(C)(H)3', 'C(C)(H)(Pt)2', 'C(C)2(H)(O)', 'C(C)(H)2(O)', 'C(C)(H)(O)(Pt)', 'C(C)2(O)(Pt)', 'C(C)(H)(=O)', 'C(C)2(=O)', 'C(C)(Pt)(=O)', 'C(C)2(H)2', 'C(C)(O)(Pt)2', 'C(C)2(Pt)2', 'C(C)2(H)(Pt)']
    )
    ### 9 permittable features for the Ru database
    ga2_ru_features=np.array(
    ['C(C)(H)2(Ru)', 'C(C)(H)3', 'C(C)2(H)2', 'C(C)(H)(Ru)2', 'C(C)2(H)(Ru)', 'C(C)2(Ru)2', 'C(C)(Ru)3', 'C(C)3(H)', 'C(C)3(Ru)']
    )
    ### 13 permittable features for the Ir database
    ga2_ir_features=np.array(
    ['C(C)(Ir)3', 'C(C)(H)(Ir)2', 'C(C)2(Ir)2', 'C(C)2(H)2', 'O(C)2', 'C(C)(H)2(O)', 'C(C)(O)(Ir)2', 'C(C)2(H)(Ir)', 'C(C)(H)2(Ir)', 'C(C)(H)3', 'O(C)(H)', 'O(=C)', 'O(H)(O)']
        )
    y_out = []
    group_out = []
    for i in molecs:
        molec = i
        M = Get_metal(molec)
        ### loading pre-trained models
        if M == 'Pt' or \
           M == None:
            ga2_features = ga2_pt_features
            local = './Data/gpr_pt_model.pkl'
            Gdrive = '/drive/My Drive/Colab Notebooks/METE/Data/gpr_pt_model.pkl'
            if os.path.exists(local):
                gprfilename = local
            elif os.path.exists(Gdrive):
                gprfilename = Gdrive

        elif M == 'Ru':
            ga2_features = ga2_ru_features
            local = './Data/gpr_ru_model.pkl'
            Gdrive = '/drive/My Drive/Colab Notebooks/METE/Data/gpr_ru_model.pkl'
            if os.path.exists(local):
                gprfilename = local
            elif os.path.exists(Gdrive):
                gprfilename = Gdrive
                
        elif M == 'Ir':
            ga2_features = ga2_ir_features
            local = './Data/gpr_ir_model.pkl'
            Gdrive = '/drive/My Drive/Colab Notebooks/METE/Data/gpr_ir_model.pkl'
            if os.path.exists(local):
                gprfilename = local
            elif os.path.exists(Gdrive):
                gprfilename = Gdrive

        ### convert the input to smiles, if not already
        descriptors = frag(molec, 2)
        cleaner = re.split('{|}', str(descriptors))
        cleaner.pop()     ### pop the last element
        cleaner.pop(0)    ### addtional text cleanup
        smiles = str(cleaner[0])
        text = []
        text.append(smiles)

        molec = Chem.MolFromSmiles(molec.split()[0])
        molec = Chem.MolToSmiles(molec)
        col,res = freq_mat.get_freq_mat(text,verbose)
        freq = np.zeros(len(ga2_features))
        for i, elei in enumerate(ga2_features):
            for j, elej in enumerate(col):
                if elei==elej:
                    freq[i]+=res[0][j]
                    
        all_zero = np.all(freq == 0)
        
        ### only if the molecule has 2NN
        if not all_zero:
            ### Get predictions
            ## test with X_test and y_test
            X_test = np.array([freq])
            # print(X_test)
            ### Create an instance of the GPRModel class
            ### Plus WhiteKernel(0.1)
            ### If pickle files exisits, will load trained ML model
            gpr = GPR.GPRModel(autosave=False, filename=gprfilename)
            y_pred, std = gpr.predict(X_test)
            if verbose:
                print('Groups at GA2: {}'.format(text))
            # print("Correction from GPR-GA2: {:.2f} kcal/mol,  STDEV {:.5f} \n".format(float(-y_pred), float(std)))
            ### output name, the full calibration on the thermochesmitry data, Hf0, S0, and Cp@[100,1500] K
            group_out.append(dict(descriptors)) ### groups name and occurance
            y_out.append(y_pred)
        else:
            if verbose:
                print('No GA2 contributions')
            group_out, y_out = 0 , 0
    return(group_out, y_out)
        # ### Save the model
        # ### default is False
        # if gpr.autosave:
        #     gpr.save(gprfilename)