#!/bin/bash

echo "[INFO] Creating conda environment 'smiles2peptides' with Python 3.7.16..."
conda create --yes --name smiles2peptides python=3.7.16
echo "[INFO] Environment created."

echo "[INFO] Activating the environment..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate smiles2peptides
echo "[INFO] Environment activated."

echo "[INFO] Installing RDKit from conda-forge/label/cf202003..."
pip install rdkit

echo "[INFO] Installing PyTorch from:..."
pip3 install torch torchvision

echo "[INFO] Installing Openpyxl from anaconda..."
pip install openpyxl

echo "[INFO] Installing scikit-learn..."
pip install scikit-learn

echo "[INFO] Installing ipykernel..."
pip install ipykernel

echo "[INFO] Installing pandas...
pip install pandas

echo "[INFO] Installing smiles2peptides package from GitHub..."
pip install git+https://github.com/BilodeauGroup/smiles2peptides.git

echo "[INFO] Setup complete. You have now the 'smiles2peptides' environment ready to use."