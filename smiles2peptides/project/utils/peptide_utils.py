import os
import pandas as pd
from rdkit import Chem  
from rdkit.Chem import Draw  
from IPython.display import display 
import re

import numpy as np
import torch
from IPython.display import Image, display
from collections import defaultdict
from sklearn.preprocessing import OneHotEncoder

class PeptideUtils:
    """
    Utility class for processing peptide sequences and SMILES strings.
    Includes methods for reading amino acid data, parsing sequences, and modifying SMILES.
    """
    
    @staticmethod
    def util_extract_characters(sequence):
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
    def util_removing_O_and_H(smile, index, characters, character=None):
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
    def util_handle_modifications(character, index, characters, special_smile, dictionary):
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
                    next_smile = PeptideUtils.util_removing_O_and_H(next_smile, index + 2, characters)
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
                prev_smile = PeptideUtils.util_removing_O_and_H(prev_smile, index, characters)
                return prev_smile, index
            else:
                raise ValueError(
                    f"The SMILE of the character '{characters[index - 1]}' does not start with a valid N-terminus group.")
    
    @staticmethod
    def util_forming_cycle(smile, characters):
        """
        Adjusts the SMILES string for cyclic peptides by ensuring proper nitrogen and oxygen handling.
        Args:
            smile (str): Input SMILES string.
            characters (list): List of characters in the sequence.
        Returns:
            str: Modified SMILES string for cyclic peptides.    
        
        """
        if smile.startswith("[N:1]") and characters[0] is not None:
            print(f"Warning: The amino acidS '{characters[0]}' and '{characters[-1]}' cannot form a peptide bond please check the Hidrogens saturation.")
            
        if smile.startswith("[NH2:1]"):
            smile = smile.replace("[NH2:1]", "[NH1:1]9", 1)
        elif smile.startswith("[NH2]"):
            smile = smile.replace("[NH2]", "[NH1]9", 1)
        elif smile.startswith("[NH1:1]"):
            smile = smile.replace("[NH1:1]", "[N:1]9", 1)
        elif smile.startswith("[NH3:1]"):
            smile = smile.replace("[NH3:1]", "[NH2:1]9", 1)
        
        if smile.endswith("O"):
            smile = smile.rsplit("O", 1)[0] + "9"
        
        else:
            print(f"Warning: The amino acid '{characters[-1]}'does not end with an oxygen atom, please check the structure.")
        
        return smile
    
    @staticmethod
    def util_show_molecule(mol):
        drawer = Draw.MolDraw2DCairo(800, 800)
        options = drawer.drawOptions()
        options.addAtomIndices = False
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        img = drawer.GetDrawingText()
        display(Image(data=img))
    
    @staticmethod
    def util_validate_exception_position(character, exception_value, position, last_position):
        if exception_value == 1 and position != 0:
            raise ValueError(f"Error: {character} can only be at the beginning.")
        if exception_value == 2 and position != last_position:
            raise ValueError(f"Error: {character} can only be at the end.")
        if exception_value == 3 and position != 0:
            raise ValueError(f"Error: The modification {character} can only be placed at the beginning.")
        if exception_value == 4 and position != last_position:
            raise ValueError(f"Error: The modification {character} must be the last element.")
    
    @staticmethod
    def util_atomic_features(amino_acids_mol):
        # Extracts atomic and bond features from a list of RDKit molecules.
        atom_features = {
            "atomic_number": [],
            "aromaticity": [],
            "num_bonds": [],
            "bonded_hydrogens": [],
            "hybridization": [],
            "implicit_valence": [],
        }
        bond_features = {
            "bond_type": [],
            "in_ring": [],
            "conjugated": [],
            "bond_aromatic": [],
            "valence_contribution_i": [],
            "valence_contribution_f": [],
        }
        # Accumulate features from all amino acid molecules
        for aa_mol in amino_acids_mol:
            atom_features["atomic_number"].extend([atom.GetAtomicNum() for atom in aa_mol.GetAtoms()])
            atom_features["aromaticity"].extend([int(atom.GetIsAromatic()) for atom in aa_mol.GetAtoms()])
            atom_features["num_bonds"].extend([atom.GetDegree() for atom in aa_mol.GetAtoms()])
            atom_features["bonded_hydrogens"].extend([atom.GetTotalNumHs() for atom in aa_mol.GetAtoms()])
            atom_features["hybridization"].extend([str(atom.GetHybridization()) for atom in aa_mol.GetAtoms()])
            atom_features["implicit_valence"].extend([atom.GetImplicitValence() for atom in aa_mol.GetAtoms()])
            # Bond features
            for bond in aa_mol.GetBonds():
                bond_features["bond_type"].append(bond.GetBondTypeAsDouble())
                bond_features["in_ring"].append(int(bond.IsInRing()))
                bond_features["conjugated"].append(int(bond.GetIsConjugated()))
                bond_features["bond_aromatic"].append(int(bond.GetIsAromatic()))
                bond_features["valence_contribution_i"].append(int(bond.GetValenceContrib(bond.GetBeginAtom())))
                bond_features["valence_contribution_f"].append(int(bond.GetValenceContrib(bond.GetEndAtom())))
        
        # Fit OneHotEncoders for each feature type
        def fit_encoder(values):
            encoder = OneHotEncoder()
            encoder.fit(np.array(list(set(values))).reshape(-1, 1))
            return encoder
        
        encoders = {
            "atomic_number": fit_encoder(atom_features["atomic_number"]),
            "aromaticity": fit_encoder(atom_features["aromaticity"]),
            "num_bonds": fit_encoder(atom_features["num_bonds"]),
            "bonded_hydrogens": fit_encoder(atom_features["bonded_hydrogens"]),
            "hybridization": fit_encoder(atom_features["hybridization"]),
            "implicit_valence": fit_encoder(atom_features["implicit_valence"]),
            "bond_type": fit_encoder(bond_features["bond_type"]),
            "in_ring": fit_encoder(bond_features["in_ring"]),
            "conjugated": fit_encoder(bond_features["conjugated"]),
            "bond_aromatic": fit_encoder(bond_features["bond_aromatic"]),
            "valence_contribution_i": fit_encoder(bond_features["valence_contribution_i"]),
            "valence_contribution_f": fit_encoder(bond_features["valence_contribution_f"]),
        }
        
        # Create node features dictionary
        node_features_dict = defaultdict(list)
        for atom, aromatic, bonds, hydrogen, hybrid, impli_vale in zip(
                                                                        atom_features["atomic_number"],
                                                                        atom_features["aromaticity"],
                                                                        atom_features["num_bonds"],
                                                                        atom_features["bonded_hydrogens"],
                                                                        atom_features["hybridization"],
                                                                        atom_features["implicit_valence"]
                                                                    ):
            
            node_key = f"{atom}_{aromatic}_{bonds}_{hydrogen}_{hybrid}_{impli_vale}"
            
            feature_node = np.concatenate([
                                            encoders["atomic_number"].transform([[atom]]).toarray()[0],
                                            encoders["aromaticity"].transform([[aromatic]]).toarray()[0],
                                            encoders["num_bonds"].transform([[bonds]]).toarray()[0],
                                            encoders["bonded_hydrogens"].transform([[hydrogen]]).toarray()[0],
                                            encoders["hybridization"].transform([[hybrid]]).toarray()[0],
                                            encoders["implicit_valence"].transform([[impli_vale]]).toarray()[0],
                                        ])
            
            # Store the feature vector in the dictionary
            node_features_dict[node_key] = feature_node
            
        # Create edge features dictionary
        edge_features_dict = defaultdict(list)
        for bond, ring, conjugat, aroma, valence_i, valence_f in zip(
                                                                    bond_features["bond_type"],
                                                                    bond_features["in_ring"],
                                                                    bond_features["conjugated"],
                                                                    bond_features["bond_aromatic"],
                                                                    bond_features["valence_contribution_i"],
                                                                    bond_features["valence_contribution_f"]
                                                                ):
            edge_key = f"{bond:.1f}_{ring:.1f}_{conjugat:.1f}_{aroma:.1f}_{valence_i:.1f}_{valence_f:.1f}"
            
            feature_edge = np.concatenate([
                                            encoders["bond_type"].transform([[bond]]).toarray()[0],
                                            encoders["in_ring"].transform([[ring]]).toarray()[0],
                                            encoders["conjugated"].transform([[conjugat]]).toarray()[0],
                                            encoders["bond_aromatic"].transform([[aroma]]).toarray()[0],
                                            encoders["valence_contribution_i"].transform([[valence_i]]).toarray()[0],
                                            encoders["valence_contribution_f"].transform([[valence_f]]).toarray()[0],
                                        ])
            
            # Store the feature vector in the dictionary
            edge_features_dict[edge_key] = feature_edge
        
        return node_features_dict, edge_features_dict
    
    @staticmethod
    def util_extract_node_and_edge_keys(mol):
        """
        Extracts node and edge key features from an RDKit molecule.
        
        Parameters:
            mol (rdkit.Chem.Mol): An RDKit molecule object.\
        
        Returns:
            node_keys_features (list of str): Encoded string keys for atom-level features.
            edge_key_features (list of str): Encoded string keys for bond-level features.
        """
        # Atom-level (node) features
        atomic_number = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
        aromaticity = [int(atom.GetIsAromatic()) for atom in mol.GetAtoms()]
        num_bonds = [atom.GetDegree() for atom in mol.GetAtoms()]
        bonded_hydrogens = [atom.GetTotalNumHs() for atom in mol.GetAtoms()]
        hybridization = [str(atom.GetHybridization()) for atom in mol.GetAtoms()]
        implicit_valence = [atom.GetImplicitValence() for atom in mol.GetAtoms()]
        
        node_keys_features = [
                            f"{atomic}_{aromatic}_{bonds}_{hydrogen}_{hybrid}_{impli_vale}"
                            for atomic, aromatic, bonds, hydrogen, hybrid, impli_vale in zip(
                                                                                            atomic_number,
                                                                                            aromaticity,
                                                                                            num_bonds,
                                                                                            bonded_hydrogens,
                                                                                            hybridization,
                                                                                            implicit_valence
                                                                                        )
                            ]
        
        # Bond-level (edge) features
        edge_keys_features = []
        for bond in mol.GetBonds():
            bond_type = bond.GetBondTypeAsDouble()
            in_ring = int(bond.IsInRing())
            conjugated = int(bond.GetIsConjugated())
            bond_aromatic = int(bond.GetIsAromatic())
            valence_contribution_i = int(bond.GetValenceContrib(bond.GetBeginAtom()))
            valence_contribution_f = int(bond.GetValenceContrib(bond.GetEndAtom()))
            
            edge_key = f"{bond_type:.1f}_{in_ring:.1f}_{conjugated:.1f}_{bond_aromatic:.1f}_{valence_contribution_i:.1f}_{valence_contribution_f:.1f}"
            edge_keys_features.append(edge_key)
        
        return node_keys_features, edge_keys_features
    
    @staticmethod
    def util_atomic_features_tensors(node_keys_features, edge_key_features, node_ft_dict, edge_ft_dict, device):
        """
        Builds PyTorch tensors for node and edge features using provided feature dictionaries and keys.
        
        Parameters:
            node_keys_features (list of str): Keys for node features.
            edge_key_features (list of str): Keys for edge features.
            node_ft_dict (dict): Dictionary mapping node keys to feature arrays.
            edge_ft_dict (dict): Dictionary mapping edge keys to feature arrays.
            device (str or torch.device): Device to place the tensors on ('cpu' or 'cuda').
        
        Returns:
            nodes_features (torch.Tensor): Tensor of shape [num_nodes, node_feature_dim].
            edges_features (torch.Tensor): Tensor of shape [num_edges, edge_feature_dim].
        """
        
        
        missing_node_keys = [key for key in node_keys_features if key not in node_ft_dict]
        missing_edge_keys = [key for key in edge_key_features if key not in edge_ft_dict]
        
        if missing_node_keys:
            raise KeyError(
                f"Missing node keys: {missing_node_keys}. "
                "Node features not found in the library follow this format: "
                "{atomic_number}_{aromatic_atom_flag}_{number_of_bonds}_{number_of_hydrogens}_{hybridization}_{implicit_valence}. "
                "Please add examples for the missing keys."
            )
        
        if missing_edge_keys:
            raise KeyError(
                f"Missing edge keys: {missing_edge_keys}. "
                "Edge features not found in the library follow this format: "
                "{bond_type}_{in_ring_flag}_{conjugated}_{aromatic_flag}_{valence_contribution_to_atom_i}_{valence_contribution_to_atom_f}. "
                "Please add examples for the missing keys."
            )
        
        nodes_features = torch.tensor(
                                        np.array([node_ft_dict[key] for key in node_keys_features]),
                                        dtype=torch.float32,
                                        device=device
                                    )
        edges_features = torch.tensor(
                                        np.array([edge_ft_dict[key] for key in edge_key_features]),
                                        dtype=torch.float32,
                                        device=device
                                    )
        
        return nodes_features, edges_features
    
    @staticmethod
    def util_atomic_adjacency_matrix(mol, device):
        """
        Constructs an adjacency matrix for the atoms in an RDKit molecule.
        Parameters:
            mol (rdkit.Chem.Mol): An RDKit molecule object.
            device (str or torch.device): Device to place the tensor on ('cpu' or 'cuda').
        Returns:
            torch.Tensor: Adjacency matrix of shape [num_atoms, num_atoms]. 
        """
        edges=[]
        for bond in mol.GetBonds():
            edges.append((bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()))
        
        graph_edges = [[x[0] for x in edges],[x[1] for x in edges]]
        
        return torch.tensor(graph_edges, dtype=torch.long, device=device)
