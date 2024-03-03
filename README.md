# METE
Machine-learning Enabled Thermochemistry Estimator (METE)

METE is a python package that takes in a sequence of SMILES strings and ouput thermodynamical energetics.
Including Standard Enthalpy, Standard Entropy, Heat Capacities from [100-1500] K.
METE can also plot interactive plots of the quantaties above as a function of teperature.

# Dependency

RDKit

# Use the Jupytor Notebook(Recommended)
open the 'METE_local' file

# Use Python IDE (e.g., Spyder)
## Getting numerical results
cd ./funcs
import main
smiles = ['OC([Pt])([Pt])C[Pt]', 'CCC[Ir]CO','CC(C)C[Ru]']
thermo = main.predict(smiles,verbose=False) (Default verbose is off)
## Plotting
import plotting
nasa = plotting.nasa9fit(thermo[0][0])
nasa.C(800)
nasa.H(600)
nasa.S(500)
nasa.plotC()


# Citation
Tianjun Xie, Gerhard R. Wittreich, Matthew T. Curnan, Geun Ho Gub, Kayla N. Seals, Justin S. Tolbert. Journal of Chemical Information and Modeling 2024. Submitted.
