from setuptools import setup, find_packages

# Leer el contenido de README.md si está disponible
try:
    with open("README.md", "r", encoding="utf-8") as f:
        long_description = f.read()
except FileNotFoundError:
    long_description = "Un módulo para convertir SMILES a péptidos y análisis de estructuras moleculares"

setup(
    name="smiles2peptides", 
    version="0.1",  
    packages=find_packages(),  
    install_requires=[  
        "pandas>=1.3.5",
        "rdkit>=2023.3.2",
    ],
    include_package_data=True,  # Asegura que se incluyan archivos no Python
    package_data={
        "smiles2peptides": ["amino_acid_library.xlsx"],  # Especifica el archivo a incluir
    },
    author="Daniel Garzon Otero",
    author_email="vvd9fd@virginia.edu",
    description="Un módulo para convertir SMILES a péptidos y análisis de estructuras moleculares",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/danielgarzonotero/smiles2peptides.git",
    python_requires='>=3.7.16',
)


