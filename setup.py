from setuptools import setup, find_packages

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
        "smiles2peptides": ["data.xlsx"],  # Especifica el archivo a incluir
    },
    author="Daniel Garzon Otero",
    author_email="vvd9fd@virginia.edu",
    description="Un módulo para convertir SMILES a péptidos y análisis de estructuras moleculares",
    long_description=open('README.md').read() if 'README.md' in locals() else 'Descripción del paquete',
    long_description_content_type="text/markdown",
    url="https://github.com/danielgarzonotero/smiles2peptides.git",
    python_requires='>=3.7.16',
)
