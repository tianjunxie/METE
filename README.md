# METE
Machine-learning Enabled Thermochemistry Estimator (METE)

# Dependency

RDKit

# Use the Jupytor Notebook(Recommended)
open the 'METE_local' file

# Use Python IDE (e.g., Spyder)
cd ./funcs
import main
smiles = ['OC([Pt])([Pt])C[Pt]', 'CCC[Ir]CO','CC(C)C[Ru]']
main.predict(smiles,verbose=False)


# Todos
1. Visualization of H, Cp, and S curve 
2. train module fix-up
3. convert serveral functions to classes
