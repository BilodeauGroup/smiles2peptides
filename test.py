#%%
from smiles2peptides.interface import Smiles2Peptide

smiles2pep = Smiles2Peptide()

# Build a peptide molecule from a sequence
mol = smiles2pep.get_peptide("ACDWY{am}")
smiles2pep.get_plot_aminoacids(mol, highlight_bonds=True)


# %%
