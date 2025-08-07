# Smiles2Peptides



## Installation Guide

To create the environment with all the required packages, simply run the following script:

```sh
bash setup_smiles2peptides.shp_env.sh 
```

### 1. Create a Conda Environment - Manually
Alternatively, you can create the environment step by step by running the following commands manually in the terminal:

```sh
conda create --name smiles2peptides python=3.7.16
```

### 2. Activate the Environment

```sh
conda activate smiles2peptides
```

### 3. Install Dependencies

```sh
conda install conda-forge/label/cf202003::rdkit
```
```sh
conda install pytorch::pytorch
```
```sh
conda install anaconda::openpyxl
```
```sh
conda install conda-forge::scikit-learn
```
```sh
conda install anaconda::ipykernel
```
```sh
pip install https://github.com/BilodeauGroup/smiles2peptides.git
```

## Usage

Once installed, you can import and use the package in your Python scripts:

```python
from smiles2peptides.interface import Smiles2Peptide
```

## About the Smiles2Peptide class
Smiles2Peptide is the central interface for building peptide molecules and analyzing amino acids using RDKit. It combines the functionalities of peptide construction and amino acid processing through internal modules.

## Instantiating the class
Before using it, you need to create an instance of the class:

```python
smiles2peptide = Smiles2Peptide()
```
Optionally, you can provide a custom amino acid dictionary by passing the path to an Excel file:
```python
smiles2pep = Smiles2Peptide(custom_dict_path="path/to/custom_amino_acids.xlsx")
```

## About the amino acid dictionary
The dictionary contains the amino acids definitions used for building peptides. It can be:

Default: Loads the standard amino acid dictionary included with the package.

Custom: You can provide your own Excel file as a custom dictionary.

Important: The custom dictionary must follow the expected structure, with amino acids defined in the CHUCKLES format, including Map Numbers

This structure ensures the peptide builder can interpret and construct molecules correctly.

## Example usage:
```python
# Instantiate (load default amino acid dictionary)
smiles2pep = Smiles2Peptide()

# Build a peptide molecule from a sequence
mol = smiles2pep.get_peptide("ACDEFGHIKLMNPQRSTVWY")

# Show atomic features tensor shape
smiles2pep.describe_peptide_atomic_features(mol)

# Visualize amino acid fragments in the peptide
smiles2pep.get_plot_aminoacids(mol, highlight_bonds=True)

# Get adjacency matrix at atomic level
adj_matrix = smiles2pep.get_peptide_atomic_adjacency_matrix(mol)

# Get amino acid features tensor
aa_features = smiles2pep.get_amino_acid_features("ACDEFGHIKLMNPQRSTVWY")
```


## Peptide Notation



- **Non-natural amino acids** are enclosed in `{Xyz}`.

- **Modifications** such as acetylation and amidation are also enclosed in `{}`, e.g., `{ac}` for acetylation and `{am}` for amidation.

For a full list of supported amino acids, refer to **`amino_acid_library.xlsx`**.

## Author

[Daniel Garzón Otero](https://github.com/danielgarzonotero)
