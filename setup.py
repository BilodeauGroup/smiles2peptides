from setuptools import setup, find_packages

setup(
    name="smiles2peptides", 
    version="0.1",  
    packages=find_packages(),  
    install_requires=[  
        "pandas",  # Especificando la versión para evitar incompatibilidades
        "rdkit",  # Igualmente, especifica la versión exacta de RDKit
    ],
    entry_points={
        "console_scripts": [
            "smiles2peptides=smiles2peptides.main:main",  
        ],
    },
    author="Daniel Garzon Otero",
    author_email="vvd9fd@virginia.edu",
    description="Un módulo para convertir SMILES a péptidos y análisis de estructuras moleculares",
    long_description=open('README.md').read(),  
    long_description_content_type="text/markdown",  
    url="https://github.com/danielgarzonotero/smiles2peptides.git",  
    python_requires='>=3.7.16',  # Asegurarse de que funcione con Python 3.7 o superior
)

