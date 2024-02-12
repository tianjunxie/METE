import ga2
import json
import numpy as np
import os

def eval(mols, verbose):
    """
    Defragment molecules and evaluate using ga lvl1 data
    
    Parameters
    ----------
    mols : list str
        list of molecules in SMILES string format

    Returns
    -------
    output : list float
        list of stand enthalpy, entropy, and NASA9 coefficient for heat capacity
             list float
        raw list of tand enthalpy, entropy, heat capacities at [100, 1500] K

    """
    output= []
    raw = []
    for mol in mols:
        ### fragmentation at lvl1
        groups_1 = ga2.frag(mol, 1)
        if verbose:
            print('Groups at GA1: {}'.format(groups_1))
        ### getting json file
        if "Pt" in mol:
            local = './Data/Pt.json'
            Gdrive = '/drive/My Drive/Colab Notebooks/METE/Data/Pt.json'
            if os.path.exists(local):
                json_file = local
            elif os.path.exists(Gdrive):
                json_file = Gdrive
            
        elif "Ru" in mol:
            local = './Data/Ru.json'
            Gdrive = '/drive/My Drive/Colab Notebooks/METE/Data/Ru.json'
            if os.path.exists(local):
                json_file = local
            elif os.path.exists(Gdrive):
                json_file = Gdrive

        elif "Ir" in mol:
            local = './Data/Ir.json'
            Gdrive = '/drive/My Drive/Colab Notebooks/METE/Data/Ir.json'
            if os.path.exists(local):
                json_file = local
            elif os.path.exists(Gdrive):
                json_file = Gdrive
        ### initialization
        Ho = 0
        So = 0
        a0 = 0
        a1 = 0
        a2 = 0
        a3 = 0
        a4 = 0
        a5 = 0
        a6 = 0
        thermo = np.zeros(17)
        ### load according JSON entries
        with open(json_file, 'r') as f:
            ga1_data = json.load(f)
            # print(mol)
            ### loop over the groups within molecule
            for group_key, occur_num in groups_1.items():
                if group_key not in ga1_data:
                    print(f'Warning, no {group_key} found in current database!')
                    break
                else:
                    ### loop over the groups within ga lvl1 data
                    for key, value in ga1_data.items():
                        ### accounting group contributions
                        if group_key == key:
                            Ho += value['Ho'] * occur_num
                            So += value['So'] * occur_num
                            a0 += value['a0'] * occur_num
                            a1 += value['a1'] * occur_num
                            a2 += value['a2'] * occur_num
                            a3 += value['a3'] * occur_num
                            a4 += value['a4'] * occur_num
                            a5 += value['a5'] * occur_num
                            a6 += value['a6'] * occur_num
                            thermo += np.array(value['raw']) * occur_num
                            # print(key, value['a0'], occur_num)
            output.append((Ho,So,a0,a1,a2,a3,a4,a5,a6))
            raw.append((thermo))
    return(output, raw)
