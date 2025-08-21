#!/bin/bash
set -e  # Exit immediately if a command exits with a non-zero status
set -o pipefail

ENV_NAME="smiles2peptides"
PYTHON_VERSION="3.7.16"

echo "[INFO] Creating conda environment '$ENV_NAME' with Python $PYTHON_VERSION..."
conda create --yes --name $ENV_NAME python=$PYTHON_VERSION
echo "[INFO] Environment created."

echo "[INFO] Activating the environment..."
# Initialize conda for the current shell
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate $ENV_NAME
echo "[INFO] Environment activated."

echo "[INFO] Installing RDKit..."
pip install rdkit || { echo "[ERROR] RDKit installation failed"; exit 1; }

echo "[INFO] Installing PyTorch..."
pip install torch torchvision || { echo "[ERROR] PyTorch installation failed"; exit 1; }

echo "[INFO] Installing Openpyxl..."
pip install openpyxl || { echo "[ERROR] Openpyxl installation failed"; exit 1; }

echo "[INFO] Installing scikit-learn..."
pip install scikit-learn || { echo "[ERROR] scikit-learn installation failed"; exit 1; }

echo "[INFO] Installing ipykernel..."
pip install ipykernel || { echo "[ERROR] ipykernel installation failed"; exit 1; }

echo "[INFO] Installing pandas..."
pip install pandas || { echo "[ERROR] pandas installation failed"; exit 1; }

echo "[INFO] Installing smiles2peptides package from GitHub..."
pip install git+https://github.com/BilodeauGroup/smiles2peptides.git || { echo "[ERROR] smiles2peptides installation failed"; exit 1; }

echo "[INFO] Setup complete. The '$ENV_NAME' environment is ready to use."
echo "[INFO] To activate the environment, run: conda activate $ENV_NAME"
