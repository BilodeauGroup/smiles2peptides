from setuptools import setup, find_packages

setup(
    name="smiles2peptides", 
    version="0.1",  
    packages=find_packages(include=["smiles2peptides", "smiles2peptides.*", "data", "functions"]),  # Asegúrate de incluir todos los subdirectorios relevantes
    install_requires=[  
        "pandas>=1.3.5",  # Especifica una versión para evitar incompatibilidades
        "rdkit>=2023.3.2",  # Especifica una versión para evitar incompatibilidades
    ],
    author="Daniel Garzon Otero",
    author_email="vvd9fd@virginia.edu",
    description="Un módulo para convertir SMILES a péptidos y análisis de estructuras moleculares",
    long_description=open('README.md').read() if 'README.md' in locals() else 'Descripción del paquete',  
    long_description_content_type="text/markdown",  
    url="https://github.com/danielgarzonotero/smiles2peptides.git",  
    python_requires='>=3.7.16',  # Asegúrate de que funcione con Python 3.7 o superior
)
