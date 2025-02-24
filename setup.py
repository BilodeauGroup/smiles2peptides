from setuptools import setup, find_packages

try:
    with open("README.md", "r", encoding="utf-8") as f:
        long_description = f.read()
except FileNotFoundError:
    long_description = "This library is to obtain the RDKit molecule from the amino acid sequence including the non-natural ones present in the library"

setup(
    name="smiles2peptides", 
    version="0.1",  
    packages=find_packages(),  
    install_requires=[  
        "pandas>=1.3.5",
        "rdkit>=2023.3.2",
    ],
    include_package_data=True, 
    package_data={
        "smiles2peptides": ["amino_acid_library.xlsx"], 
    },
    author="Daniel Garzon Otero",
    author_email="vvd9fd@virginia.edu",
    description="A module to obtain the RDKit molecule from an amino acid sequence",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/danielgarzonotero/smiles2peptides.git",
    python_requires='>=3.7.16',
)


