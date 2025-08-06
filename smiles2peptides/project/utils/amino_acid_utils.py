from rdkit import Chem  
from rdkit.Chem import Draw  
from IPython.display import display 
import torch
import numpy as np

from project.utils.peptide_utils import PeptideUtils


class AminoAcidUtils:
    """
    Utility class for handling peptide bonds in RDKit molecules.
    Provides methods to identify, highlight, and display peptide bonds.
    
    """
    @staticmethod
    def util_peptide_bonds(mol):
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
    def util_highlight_peptide_bonds(mol, peptidic_bonds):
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
    def util_display_amino_acids(mol, peptidic_bonds):
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
    def util_amino_acid_mapping_vector(mol, peptidic_bonds,device):
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
        
        amino_acid_fragments_vector = torch.tensor(peptide_idx, dtype=torch.long, device = device)
        
        return amino_acid_fragments_vector
    
    @staticmethod
    def util_amino_acid_adjacency_matrix(sequence, device, architecture):
        """
        Constructs an adjacency matrix for a peptide sequence, representing linear connections between amino acids.
        Args:
            sequence (str): Peptide sequence string.
            device (torch.device): Device to store the resulting tensor (e.g., 'cpu' or 'cuda').
            linear (bool): If True, constructs a linear adjacency matrix. Defaults to True.
        Returns:
            torch.Tensor: Adjacency matrix representing connections between amino acids.
        Raises:
            ValueError: If the sequence is empty or contains invalid characters.
        Raises:
            TypeError: If the sequence is not a string.
        
        """
        if not isinstance(sequence, str):
            raise TypeError("Input sequence must be a string.")
        
        sequence = sequence.replace('\u200b', '').replace(' ', '')
        characters = PeptideUtils.util_extract_characters(sequence)
        num_amino_acids= len(characters)
        
        if num_amino_acids == 0:
            raise ValueError("Input sequence is empty or contains no valid amino acids.")
        
        if architecture == 'linear':
            edges = []
            for i in range(num_amino_acids - 1):
                edges.append((i, i + 1))
            graph_edges = [[x[0] for x in edges], [x[1] for x in edges]]
        
        elif architecture == 'cyclic':
            edges = []
            for i in range(num_amino_acids - 1):
                edges.append((i, i + 1))
            # Connect the last amino acid to the first to form a cycle
            if num_amino_acids > 1:
                edges.append((num_amino_acids - 1, 0))
            graph_edges = [[x[0] for x in edges], [x[1] for x in edges]]
        
        #TODO: Staples?
        
        else:
            raise ValueError(f"Unknown architecture: {architecture}")
        
        return torch.tensor(graph_edges, dtype=torch.long, device = device) 
    
    @staticmethod
    def util_amino_acid_features_tensor(sequence, amino_acid_library, device):
        """
        Constructs a tensor representing the chirality of each amino acid in a peptide sequence.
        Args:
            sequence (str): Peptide sequence string.
            amino_acid_library (AminoAcidDictionary, optional): Dictionary containing amino acid data.
            device (str or torch.device): Device to store the resulting tensor (e.g., 'cpu' or 'cuda').
        Returns:
            torch.Tensor: Tensor of shape (num_amino_acids, 2) where each row represents the chirality of an amino acid.
        Raises:
            ValueError: If the sequence is empty or contains invalid characters.
            TypeError: If the sequence is not a string.
        """
        sequence = sequence.replace('\u200b', '').replace(' ', '')
        characters = PeptideUtils.util_extract_characters(sequence)
        
        chirality = []
        
        for aminoacid in characters:
            if aminoacid in amino_acid_library:
                config = amino_acid_library[aminoacid][0]
                if config == 'L':
                    chirality.append([1, 0])  # L -> (1, 0)
                elif config == 'D':
                    chirality.append([0, 1])  # D -> (0, 1)
            else:
                print(f"Amino acid {aminoacid} not found in dictionary.")
                chirality.append([0, 0])  # Default value if not found
                
        #TODO Anadir mas features?
        
        return torch.tensor(np.array(chirality, dtype=np.float32), device=device)