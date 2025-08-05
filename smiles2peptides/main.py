#%%
import os
import pandas as pd
from rdkit import Chem  
from rdkit.Chem import Draw  
from IPython.display import display 
import re

import sys
import numpy as np
import torch
from IPython.display import Image, display


import os
import pandas as pd

class AminoAcidDictionary:
    """
    Loads and provides access to amino acid data stored in an Excel file.
    
    Each amino acid entry contains:
        - Chirality: The chirality configuration of the amino acid.
        - SMILES: The SMILES representation of the amino acid.
        - Exception: Special flags or notes regarding the amino acid.
    
    Usage:
        - Initialize the class optionally specifying the path to the Excel file.
        - Query the dictionary for amino acid data by notation.
        - Check existence of an amino acid notation.
        - Iterate over all entries or get a full copy of the dictionary.
    """
    
    def __init__(self, dictionary_path_file=None):
        """
        Initializes the dictionary by loading data from the specified Excel file.
        If no file path is provided, loads from a default file named 'amino_acids_v4.xlsx' 
        located in the same directory as this script.
        
        Args:
            excel_path (str, optional): Path to the Excel file containing amino acid data.
        """
        
        if dictionary_path_file is None:
            base_dir = os.path.dirname(os.path.abspath(__file__))
            dictionary_path_file = os.path.join(base_dir, "amino_acids_v4.xlsx")
        
        self.excel_path = dictionary_path_file
        self._dictionary = self._load_dictionary()
    
    def _load_dictionary(self):
        """
        Reads the Excel file and builds an internal dictionary.
        
        Returns:
            dict: Mapping from amino acid notation (str) to a tuple containing (Chirality (str), SMILES (str), Exception (str)).
        """
        df = pd.read_excel(self.excel_path)
        
        required_columns = {'Notation', 'Chirality', 'SMILES', 'Exception'}
        if not required_columns.issubset(df.columns):
            raise ValueError(f"Excel file must contain columns: {required_columns}")
        
        dictionary = {}
        for _, row in df.iterrows():
            notation = str(row['Notation']).strip()
            chirality = str(row['Chirality']).strip() if pd.notnull(row['Chirality']) else ""
            smiles = str(row['SMILES']).strip() if pd.notnull(row['SMILES']) else ""
            exception = str(row['Exception']).strip() if pd.notnull(row['Exception']) else ""
            
            if notation:
                dictionary[notation] = (chirality, smiles, exception)
        
        return dictionary
    
    def get(self, notation):
        """
        Retrieves the data tuple for a given amino acid notation.
        
        Args:
            notation (str): The amino acid notation to look up.
        
        Returns:
            tuple: (Chirality, SMILES, Exception) if found, else None.
        """
        return self._dictionary.get(notation)
    
    def has(self, notation):
        """
        Checks if a notation exists in the dictionary.
        
        Args:
            notation (str): The amino acid notation to check.
        
        Returns:
            bool: True if exists, else False.
        """
        return notation in self._dictionary
    
    def items(self):
        """
        Returns an iterable of all items in the dictionary.
        
        Returns:
            iterable: Yields tuples of (notation, (Chirality, SMILES, Exception))
        """
        return self._dictionary.items()
    
    def as_dict(self):
        """
        Returns a shallow copy of the internal amino acid dictionary.
        
        Returns:
            dict: Copy of the internal dictionary mapping notation to data tuples.
        """
        return self._dictionary.copy()
    
    def __contains__(self, key):
        """
        Implements the 'in' operator for checking if a notation exists in the dictionary.
        
        Args:
            key (str): Amino acid notation.
        
        Returns:
            bool: True if the notation exists, False otherwise.
        """
        return key in self.dictionary


class PeptideUtils:
    """
    Utility class for processing peptide sequences and SMILES strings.
    Includes methods for reading amino acid data, parsing sequences, and modifying SMILES.
    """
    
    @staticmethod
    def extract_characters(sequence):
        """
        Extracts characters or residue blocks from a sequence.
        Groups inside `{}` are treated as single units.
        
        Args:
            sequence (str): Input sequence string.
        
        Returns:
            list: List of characters and grouped residues.
        """
        pattern = r'\{[^}]*\}|[a-zA-Z]'
        return re.findall(pattern, sequence)
    
    @staticmethod
    def removing_O_and_H(smile, index, characters, character=None):
        """
        Adjusts the SMILES string by removing terminal oxygen and modifying nitrogen groups.
        
        Args:
            smile (str): Input SMILES string.
            index (int): Position in the sequence.
            characters (list): List of characters in the sequence.
            character (str, optional): Current residue (used for warnings).
        
        Returns:
            str: Modified SMILES string.
        """
        if index < len(characters) and smile.endswith("O"):
            smile = smile[:-1]
        
        if index > 1:
            if smile.startswith("[N:1]") and character is not None:
                print(f"Warning: The amino acid '{character}' cannot form a peptide bond because the nitrogen's bonds are already saturated.")
            if smile.startswith("[NH2:1]"):
                smile = smile.replace("[NH2:1]", "[NH1:1]", 1)
            elif smile.startswith("[NH2]"):
                smile = smile.replace("[NH2]", "[NH1]", 1)
            elif smile.startswith("[NH1:1]"):
                smile = smile.replace("[NH1:1]", "[N:1]", 1)
            elif smile.startswith("[NH3:1]"):
                smile = smile.replace("[NH3:1]", "[NH2:1]", 1)
            

        return smile
    
    @staticmethod
    def handle_special_case(character, index, characters, special_smile, dictionary):
        """
        Handles special SMILES modifications (e.g., PEG, acetyl groups) at N-terminus or mid-chain.
        
        Args:
            character (str): Current character in sequence.
            index (int): Position in sequence.
            characters (list): Entire character list.
            special_smile (str): SMILES fragment to insert.
            dictionary (dict): Amino acid dictionary.
        
        Returns:
            tuple: (Modified SMILES, updated index)
        """
        if character not in dictionary:
            raise ValueError(f"The character '{character}' is not found in the dictionary.")
        
        n_pattern = r"^\[(NH2|NH1|N|NH3):1\]|\[NH2\]"
        
        if index == 0:
            if index + 1 < len(characters) and characters[index + 1] in dictionary:
                next_smile = dictionary[characters[index + 1]][1]
                match = re.match(n_pattern, next_smile)
                if match:
                    prefix = match.group(0)
                    rest = next_smile[len(prefix):]
                    next_smile = prefix + special_smile + rest
                    next_smile = PeptideUtils.removing_O_and_H(next_smile, index + 2, characters)
                    return next_smile, index + 1
                else:
                    raise ValueError(
                        f"The SMILE of the character '{characters[index + 1]}' does not start with a valid N-terminus group.")
            else:
                raise ValueError(f"There is no valid character after '{character}'.")
        else:
            prev_smile = dictionary[characters[index - 1]][1]
            if re.match(n_pattern, prev_smile):
                prev_smile = prev_smile[:1] + special_smile + prev_smile[1:]
                prev_smile = PeptideUtils.removing_O_and_H(prev_smile, index, characters)
                return prev_smile, index
            else:
                raise ValueError(
                    f"The SMILE of the character '{characters[index - 1]}' does not start with a valid N-terminus group.")
                
    @staticmethod
    def show_molecule(mol):
        drawer = Draw.MolDraw2DCairo(800, 800)
        options = drawer.drawOptions()
        options.addAtomIndices = False
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        img = drawer.GetDrawingText()
        display(Image(data=img))
    
    @staticmethod
    def validate_exception_position(character, exception_value, position, last_position):
        if exception_value == 1 and position != 0:
            raise ValueError(f"Error: {character} can only be at the beginning.")
        if exception_value == 2 and position != last_position:
            raise ValueError(f"Error: {character} can only be at the end.")
        if exception_value == 3 and position != 0:
            raise ValueError(f"Error: The modification {character} can only be placed at the beginning.")

class AminoAcidUtils:
    """
    Utility class for handling peptide bonds in RDKit molecules.
    Provides methods to identify, highlight, and display peptide bonds.
    
    """
    @staticmethod
    def get_peptide_bonds(mol):
        """
        Identifies peptide bonds in a molecule based on atomMap numbers.
        Peptide bonds are defined as bonds between atoms with atomMap numbers 1 and 2.
        
        Args:
            mol (Chem.Mol): RDKit molecule object.
        
        Returns:
            list: List of bond indices that are identified as peptide bonds.
        
        Raises:
            ValueError: If the molecule does not contain atoms with atomMap numbers 1 or 2.
        """
        # Get atom indices with atomMap numbers 1 and 2 (used to identify peptide bonds)
        atoms_map1 = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomMapNum() == 1]
        atoms_map2 = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomMapNum() == 2]
        
        # Raise error if either list is empty
        if not atoms_map1:
            raise ValueError("No atoms found with atomMapNum 1 in the molecule.")
        if not atoms_map2:
            raise ValueError("No atoms found with atomMapNum 2 in the molecule.")
        
        # Identify bonds connecting atoms from map1 to map2 — i.e., peptide bonds
        peptidic_bonds = []
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if (a1 in atoms_map1 and a2 in atoms_map2) or (a1 in atoms_map2 and a2 in atoms_map1):
                peptidic_bonds.append(bond.GetIdx())
        
        return peptidic_bonds
    
    @staticmethod
    def highlight_peptide_bonds(mol, peptidic_bonds):
        """ 
        Highlights peptide bonds in a molecule by displaying them in a different color.
        
        Args:
            mol (Chem.Mol): RDKit molecule object.
            peptidic_bonds (list): List of bond indices that are peptide bonds.
        Returns:
            None: Displays the molecule with highlighted peptide bonds.
        Raises:
            ValueError: If the list of peptide bonds is empty.
        
        """
        # Raise error if either list is empty
        if not peptidic_bonds:
            raise ValueError("No peptide bonds found to highlight.")
        
        # Optional: highlight identified peptide bonds (cyan color)
        highlight_bond_colors = {idx: (0.0, 1.0, 1.0) for idx in peptidic_bonds}
        
        # Display the molecule with highlighted peptide bonds
        img = Draw.MolToImage(
                            mol,
                            size=(600, 600),
                            highlightBonds=peptidic_bonds,
                            highlightBondColors=highlight_bond_colors
                            )
        display(img)
    
    @staticmethod
    def display_amino_acids(mol, peptidic_bonds):
        """
        Displays the amino acids in a molecule by highlighting them.
        
        Args:
            mol (Chem.Mol): RDKit molecule object.
            peptidic_bonds (list): List of bond indices that are peptide bonds.
        
        Returns:
            None: Displays the molecule with highlighted amino acids.
        """
        # Fragment the molecule at the peptide bonds without adding dummy atoms
        amino_acids = Chem.FragmentOnBonds(mol, peptidic_bonds, addDummies=False)
        
        # Display the resulting fragments
        img = Draw.MolToImage(amino_acids, size=(600, 600))
        display(img)
        
        return amino_acids
    
    @staticmethod
    # Fragment the molecule at the peptide bonds
    def get_amino_acid_fragments_vector(mol, peptidic_bonds):
        """ 
        Creates a tensor that maps each atom in the molecule to its corresponding amino acid fragment index.
        Ideal for scatter processing in PyTorch.
        This method fragments the molecule at the peptide bonds and assigns each atom to its corresponding fragment index.
        Args:
            mol (Chem.Mol): RDKit molecule object.
            peptidic_bonds (list): List of bond indices that are peptide bonds. 
        Returns:
            torch.Tensor: A tensor of shape (num_atoms,) mapping each atom to its amino acid index.
        """
        amino_acids = Chem.FragmentOnBonds(mol, peptidic_bonds, addDummies=False)
        
        # Extract individual fragments (amino acids)
        fragments = list(Chem.GetMolFrags(amino_acids, asMols=True))
        
        # Create a vector assigning each atom to its corresponding fragment index
        peptide_idx = []
        for i, frag in enumerate(fragments):
            num_atoms = frag.GetNumAtoms()
            peptide_idx.extend([i] * num_atoms)

        amino_acid_fragments_vector = torch.tensor(peptide_idx, dtype=torch.long)
        
        return amino_acid_fragments_vector


class PeptideBuilder:
    """
    Class responsible for constructing a peptide molecule from a given sequence.
    Optionally displays a visual representation of the resulting molecule.
    """   
    @staticmethod
    def peptide_builder(sequence, show_display=False, amino_acid_library=None):
        sequence = sequence.replace('\u200b', '').replace(' ', '')
        characters = PeptideUtils.extract_characters(sequence)
        
        concatenated_smile = ""
        
        last_index = len(characters) - 1
        
        for i, character in enumerate(characters):
            if character not in amino_acid_library:
                raise ValueError(f"The amino acid '{character}' is not found in the dictionary.")
            
            exception_value = amino_acid_library[character][2]  # asume que ya es int o None
            
            if exception_value is not '':
                exception_value = int(float(exception_value))
                PeptideUtils.validate_exception_position(character, exception_value, i, last_index)
            
            if exception_value == 3:
                special_smile = amino_acid_library[character][1]
                next_smile, new_i = PeptideUtils.handle_special_case(character, i, characters, special_smile, amino_acid_library)
                concatenated_smile += next_smile
                i = new_i  
            elif exception_value == 4:
                # lógica pendiente
                pass
            else:
                smile = amino_acid_library[character][1]
                smile = PeptideUtils.removing_O_and_H(smile, i + 1, characters, character)
                concatenated_smile += smile.strip()
        
        concatenated_smile = concatenated_smile.replace('\u200b', '').replace(' ', '')
        
        if sequence.startswith("{FITC-Ahx}"):
            editable = Chem.EditableMol(Chem.MolFromSmiles(concatenated_smile))
            editable.RemoveBond(3, 4)
            editable.AddBond(3, 5, order=Chem.rdchem.BondType.SINGLE)
            mol = editable.GetMol()
            
        else:
            mol = Chem.MolFromSmiles(concatenated_smile)
        
        if show_display and mol:
            PeptideUtils.show_molecule(mol)
        
        return mol



class AminoAcidBuilder:
    """
    Class for handling peptide bond processing in RDKit molecules.
    Provides methods to identify peptide bonds, highlight them, and optionally fragment the molecule.
    """
    @staticmethod
    def amino_acid_plotting(mol,
                            highlight_bonds=False,
                            obtain_amino_acids=False,
                            ): 
        """
        Handles peptide bond processing in a molecule by identifying peptide bonds (atomMap 1 and 2) 
        and optionally performing visualization or fragment extraction.
        
        Parameters:
            mol (Chem.Mol): RDKit molecule object.
            highlight_bonds (bool): If True, highlights peptidic bonds in the molecule. defaults = False.
            obtain_amino_acids (bool): If True, returns the molecule fragmented at peptidic bonds. defaults = False.
        Returns:
            Chem.Mol (optional): Fragmented molecule if `obtain_amino_acids` is True.
            torch.Tensor (optional): Vector assigning each atom to a fragment index if `get_fragment_amino_acid_vector` is True.
        """
        
        peptidic_bonds = AminoAcidUtils.get_peptide_bonds(mol)
        
        if highlight_bonds:
            AminoAcidUtils.highlight_peptide_bonds(mol, peptidic_bonds)
            
        if obtain_amino_acids:
            AminoAcidUtils.display_amino_acids(mol, peptidic_bonds)
    
    @staticmethod
    def get_amino_acid_mapping_vector(mol, peptidic_bonds):
        """ Returns a tensor mapping each atom to its amino acid fragment."""
        
        peptidic_bonds = AminoAcidUtils.get_peptide_bonds(mol)
        
        return AminoAcidUtils.get_amino_acid_fragments_vector(mol, peptidic_bonds)




peptides = [
    "{pra}{FITC-Ahx}{Pra}nT",
    "{PEG2}pr{GlcNAc-T}l",
    "v{photo-L}{phospho-Y}l{4-tert-butyl-p}iK",
    "{acm-C}wM{4&5-hydrox-L}{(N-me)-a}{iso-Q}",
    "{4-hydrox-p}{seleno-C}R{trime-L}{3&4-dihydrox-F}G",
    "K{p-carboxyl-F}{photo-L}c{phospho-Y}{me-f}",
    "{nor-r}{photo-M}l{GlcNAc-S}H{me-f}",
    "i{h-s}{(N-me)-A}{Pra}w{acm-C}l",
    "d{p-carboxyl-f}{photo-L}{me-Y}K{C-me}",
    "{photo-l}{p-carboxyl-f}m{3-me-f}w{4-hydrox-p}",
    "a{Pra}{trime-L}y{photo-M}{iso-Q}",
    "{(N-me)-a}H{3-me-f}{photo-M}l{p-carboxyl-F}",
    "M{phospho-Y}{4-hydrox-P}d{acm-c}{Pra}",
    "n{4-hydrox-p}{GlcNAc-S}K{(N-me)-A}{4-tert-butyl-p}",
    "{succinyl-K}{photo-M}{4-hydrox-p}R{me-f}t",
    "G{me-F}{Dip-a}{photo-M}{p-carboxyl-f}{4-hydrox-p}",
    "{seleno-C}v{Pra}{p-carboxyl-F}{GlcNAc-T}w",
    "p{photo-L}{me-f}{trime-L}K{3&4-dihydrox-F}",
    "{4-hydrox-P}y{photo-M}{h-S}F",
    "W{p-carboxyl-f}M{photo-M}q{C-me}",
    "{phospho-S}{4-hydrox-p}l{succinyl-K}{3-me-f}V",
    "{photo-l}{seleno-C}{(N-me)-A}{4-hydrox-P}dF",
    "c{me-y}{photo-M}{3&4-dihydrox-F}{acm-C}t",
    "{GlcNAc-asn}{photo-M}n{h-S}{4-tert-butyl-p}l",
    "{iso-D}p{Pra}{me-f}H{p-carboxyl-F}",
    "{acm-c}{4-hydrox-P}G{(N-me)-A}R{photo-M}",
    "{4-hydrox-p}{photo-L}{succinyl-K}y{3-me-f}t",
    "L{me-f}{photo-M}{4&5-hydrox-L}m{iso-Q}",
    "{trime-L}H{4-hydrox-p}{photo-l}d{C-me}",
    "{Pra}{seleno-C}{4-hydrox-p}n{phospho-Y}L",
    "t{acm-C}{photo-M}{me-f}{4-tert-butyl-p}k",
    "{4-hydrox-p}V{3-me-f}{photo-l}y",
    "{4-hydrox-P}{photo-M}H{acm-c}m",
    "{succinyl-K}{photo-M}{p-carboxyl-f}{h-s}T",
    "F{3-me-f}{4-hydrox-p}{me-f}{Pra}w",
    "{photo-l}n{4&5-hydrox-L}M{trime-L}y",
    "{photo-M}{me-f}{phospho-Y}{4-hydrox-p}G{C-me}",
    "L{p-carboxyl-F}{photo-L}{Pra}h",
    "{(N-me)-a}{seleno-C}{photo-M}R{4-hydrox-p}w",
    "{phospho-S}{photo-l}{4-hydrox-P}{3-me-f}yM",
    "{photo-L}F{h-s}{GlcNAc-T}{Pra}c",
    "{trime-L}{4-hydrox-p}y{photo-M}{me-f}k",
    "R{(N-me)-A}{4-hydrox-p}{Pra}{4-tert-butyl-p}l",
    "a{photo-M}{me-F}{4-hydrox-p}{succinyl-K}n",
    "K{4-hydrox-p}{photo-l}q{(N-me)-a}",
    "{4-hydrox-P}{Pra}{acm-C}{photo-M}Lw",
    "{photo-l}G{p-carboxyl-f}{me-f}d",
    "{photo-M}{4-hydrox-p}{GlcNAc-asn}R{3-me-f}v",
    "{me-y}{photo-L}{4-hydrox-P}{succinyl-K}F",
    "{p-carboxyl-f}t{photo-M}K{h-s}{(N-me)-A}"
]



class Smiles2Peptide:
    """
    Central interface for peptide construction and amino acid analysis using RDKit.
    Combines functionality from PeptideBuilder and AminoAcidBuilder.
    """
    
    def __init__(self, custom_dict_path=None):
        
        """
        Initializes the Smiles2Peptide class with optional dictionary path.
        If no path is provided, it defaults to loading the standard amino acid dictionary.
        Args:
            dict_path (str, optional): Path to the amino acid dictionary Excel file.
        Initializes the PeptideBuilder and AminoAcidBuilder classes.
        Initializes the amino acid dictionary from the specified path or defaults to the standard dictionary.
        """
        
        self.peptide_builder = PeptideBuilder
        self.amino_acid_builder = AminoAcidBuilder
        
        if custom_dict_path:
            self.dictionary = AminoAcidDictionary(custom_dict_path).as_dict()
        else:
            self.dictionary = AminoAcidDictionary().as_dict()
        
        
    def build_peptide(self, sequence, plot_peptide=False, amino_acid_library=None):
        """
        Build a peptide molecule from a sequence.
        
        Args:
            sequence (str): Peptide sequence.
            show_display (bool): Whether to visualize the molecule.
        
        Returns:
            Chem.Mol: RDKit molecule.
        """
        return self.peptide_builder.peptide_builder(sequence, plot_peptide, amino_acid_library=self.dictionary)
    
    def plot_aminoacids(self, mol, highlight_bonds=False, obtain_amino_acids=False):
        """
        Visualize or fragment the peptide molecule.\
        
        Args:
            mol (Chem.Mol): RDKit molecule.
            highlight_bonds (bool): Highlight peptide bonds.
            obtain_amino_acids (bool): Display amino acid fragments.
        """
        return self.amino_acid_builder.amino_acid_plotting(mol, highlight_bonds, obtain_amino_acids)
    
    def get_fragment_mapping(self, mol):
        """
        Get a tensor mapping each atom to its amino acid fragment.
        
        Args:
            mol (Chem.Mol): RDKit molecule.
        
        Returns:
            torch.Tensor: Mapping vector.
        """
        return self.amino_acid_builder.get_amino_acid_mapping_vector(mol, None)




for pep in peptides:
    print('----------------------')
    print(pep)
    smiles2peptides = Smiles2Peptide() #aqui se instancia la clase para manejar peptidos, se puede usar ademas un diccionario personalizado
    
    mol = smiles2peptides.build_peptide(pep, plot_peptide=True)
    
    smiles2peptides.plot_aminoacids(mol, highlight_bonds=True, obtain_amino_acids=True)
    
    print("Amino Acid mapping Vector:", smiles2peptides.get_fragment_mapping(mol))
    



# %%
