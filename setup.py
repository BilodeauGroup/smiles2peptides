from setuptools import setup, find_packages

setup(
    name="smiles2peptides",
    version="0.1",
    packages=find_packages(include=["smiles2peptides"]),
    install_requires=[
        "pandas>=1.3.5",
        "rdkit>=2023.3.2",
    ],
    include_package_data=True,  # Importante para incluir archivos no Python
    package_data={"": ["data/*.xlsx"]},  # Asegura que los archivos de data/ sean incluidos
    author="Daniel Garzon Otero",
    author_email="vvd9fd@virginia.edu",
    description="Un módulo para convertir SMILES a péptidos y análisis de estructuras moleculares",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/danielgarzonotero/smiles2peptides",
    python_requires=">=3.7.16",
)
