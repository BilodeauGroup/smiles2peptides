# smiles2peptides

## Installation Guide

Follow the steps below to set up the environment and install the required dependencies for using `smiles2peptides`.

### 1. Create a Conda Environment

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
conda install anaconda::pandas
conda install anaconda::ipykernel
conda install anaconda::openpyxl
```

### 4. Install `smiles2peptides`

```sh
pip install git+https://github.com/danielgarzonotero/smiles2peptides.git
```

## Usage

Once installed, you can import and use the library in your Python scripts:

```python
import smiles2peptides
# Add usage examples here
```

## License

This project is licensed under the MIT License.

## Author

[Daniel Garzón Otero](https://github.com/danielgarzonotero)
