#!/bin/bash

echo "[INFO] Creating conda environment 'smiles2peptides' with Python 3.7.16..."
conda create --yes --name smiles2peptides python=3.7.16
echo "[INFO] Environment created."

echo "[INFO] Activating the environment..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate smiles2peptides
echo "[INFO] Environment activated."

echo "[INFO] Installing RDKit from conda-forge/label/cf202003..."
conda install --yes conda-forge/label/cf202003::rdkit

echo "[INFO] Installing PyTorch from pytorch channel..."
conda install --yes pytorch::pytorch

echo "[INFO] Installing Openpyxl from anaconda..."
conda install --yes anaconda::openpyxl

echo "[INFO] Installing scikit-learn from conda-forge..."
conda install --yes conda-forge::scikit-learn

echo "[INFO] Installing ipykernel..."
pip install ipykernel

echo "[INFO] Installing smiles2peptides package from GitHub..."
pip install git+https://github.com/danielgarzonotero/smiles2peptides.git

echo "[INFO] Setup complete. You have now the 'smiles2peptides' environment ready to use."