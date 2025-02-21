from setuptools import setup, find_packages

setup(
    name="smiles2peptides", 
    version="0.1",  
    packages=find_packages(),  
    install_requires=[  
        "pandas",  
        "rdkit",  
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
    url="https://github.com/tu_usuario/smiles2peptides",  
    python_requires='>=3.6',  # La versión mínima de Python que soporta tu módulo
)