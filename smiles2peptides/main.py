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
from collections import defaultdict
from sklearn.preprocessing import OneHotEncoder


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
            dict: Mapping from amino acid notation (str) to a tuple containing (Chirality (str), SMILES (str), Exception (str), Type (str)).
        """
        df = pd.read_excel(self.excel_path)
        
        required_columns = {'Notation', 'Chirality', 'SMILES', 'Exception', 'Type'}
        if not required_columns.issubset(df.columns):
            raise ValueError(f"Excel file must contain columns: {required_columns}")
        
        dictionary = {}
        for _, row in df.iterrows():
            notation = str(row['Notation']).strip()
            chirality = str(row['Chirality']).strip() if pd.notnull(row['Chirality']) else ""
            smiles = str(row['SMILES']).strip() if pd.notnull(row['SMILES']) else ""
            exception = str(row['Exception']).strip() if pd.notnull(row['Exception']) else ""
            type_ = str(row['Type']).strip() if pd.notnull(row['Type']) else ""
            
            if notation:
                dictionary[notation] = (chirality, smiles, exception, type_)
        
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
    
    def get_all_amino_acids_notations(self):
        """
        Returns a list of all amino acid notations in the dictionary where Type is 'Amino Acid'.
        Returns:
            list: List of notation strings.
        """
        return [
                k for k, v in self._dictionary.items()
                if len(v) > 3 and v[3] == "Amino Acid"
                ]


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
    def util_handle_special_case(character, index, characters, special_smile, dictionary):
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
            atom_features["hybridization"].extend([atom.GetHybridization().real for atom in aa_mol.GetAtoms()])
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
        hybridization = [atom.GetHybridization().real for atom in mol.GetAtoms()]
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
    def util_atomic_features_tensors(node_keys_features, edge_key_features, node_ft_dict, edge_ft_dict, device="cpu"):
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
            raise KeyError(f"Missing node keys in node_ft_dict: {missing_node_keys}")
        if missing_edge_keys:
            raise KeyError(f"Missing edge keys in edge_ft_dict: {missing_edge_keys}")
        
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
    def util_amino_acid_adjacency_matrix(sequence, device, architecture='linear'):
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
    



class PeptideBuilder:
    """
    Class responsible for constructing a peptide molecule from a given sequence.
    Optionally displays a visual representation of the resulting molecule.
    """   
    @staticmethod
    def builder_peptide(sequence, show_display=False, amino_acid_library=None):
        sequence = sequence.replace('\u200b', '').replace(' ', '')
        characters = PeptideUtils.util_extract_characters(sequence)
        
        concatenated_smile = ""
        
        last_index = len(characters) - 1
        
        for i, character in enumerate(characters):
            if character not in amino_acid_library:
                raise ValueError(f"The amino acid '{character}' is not found in the dictionary.")
            
            exception_value = amino_acid_library[character][2]  # asume que ya es int o None
            
            if exception_value is not '':
                exception_value = int(float(exception_value))
                PeptideUtils.util_validate_exception_position(character, exception_value, i, last_index)
            
            if exception_value == 3:
                special_smile = amino_acid_library[character][1]
                next_smile, new_i = PeptideUtils.util_handle_special_case(character, i, characters, special_smile, amino_acid_library)
                concatenated_smile += next_smile
                i = new_i  
            elif exception_value == 4:
                # lógica pendiente
                pass
            else:
                smile = amino_acid_library[character][1]
                smile = PeptideUtils.util_removing_O_and_H(smile, i + 1, characters, character)
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
            PeptideUtils.util_show_molecule(mol)
        
        return mol
    
    @staticmethod
    def builder_atomic_features(amino_acids_mol):
        return  PeptideUtils.util_atomic_features(amino_acids_mol)
    
    @staticmethod
    def builder_peptide_atomic_features(mol):
        atoms_keys, edges_keys = PeptideUtils.util_extract_node_and_edge_keys(mol)
        return atoms_keys, edges_keys
    
    @staticmethod
    def builder_atomic_features_tensors(node_keys_features, edge_key_features, node_ft_dict, edge_ft_dict, device):
        atom_features_tensors, edge_features_tensors = PeptideUtils.util_atomic_features_tensors(
                                                                                                node_keys_features,
                                                                                                edge_key_features,
                                                                                                node_ft_dict,
                                                                                                edge_ft_dict,
                                                                                                device
                                                                                                )
        return atom_features_tensors, edge_features_tensors
    



class AminoAcidBuilder:
    """
    Class for handling peptide bond processing in RDKit molecules.
    Provides methods to identify peptide bonds, highlight them, and optionally fragment the molecule.
    """
    @staticmethod
    def builder_amino_acid_plotting(mol,
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
        
        peptidic_bonds = AminoAcidUtils.util_peptide_bonds(mol)
        
        if highlight_bonds:
            AminoAcidUtils.util_highlight_peptide_bonds(mol, peptidic_bonds)
            
        if obtain_amino_acids:
            AminoAcidUtils.util_display_amino_acids(mol, peptidic_bonds)
    
    
    @staticmethod
    def builder_amino_acid_mapping_vector(mol, peptidic_bonds, device='cpu'):
        """ Returns a tensor mapping each atom to its amino acid fragment."""
        
        peptidic_bonds = AminoAcidUtils.util_peptide_bonds(mol)
        
        return AminoAcidUtils.util_amino_acid_mapping_vector(mol, peptidic_bonds, device)
    
    
    @staticmethod
    def builder_amino_acid_adjacency_matrix(sequence, device='cpu', architecture='linear'):
        """
        Constructs an adjacency matrix for a peptide sequence, representing linear connections between amino acids.
        
        Args:
            sequence (str): Peptide sequence string.
            device (torch.device): Device to store the resulting tensor (e.g., 'cpu' or 'cuda').
            architecture (str): Type of adjacency matrix ('linear' or 'cyclic'). Defaults to 'linear'.
        
        Returns:
            torch.Tensor: Adjacency matrix representing connections between amino acids.
        """
        return AminoAcidUtils.util_amino_acid_adjacency_matrix(sequence, device, architecture)



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
            custom_dict_path (str, optional): Path to the amino acid dictionary Excel file.
        """
        self.peptide_builder = PeptideBuilder
        self.amino_acid_builder = AminoAcidBuilder
        
        # Initialize the amino acid dictionary
        self._init_amino_acid_dictionary(custom_dict_path)
        
        # Build the amino acids molecules from the dictionary
        self.amino_acids_mol = self._build_amino_acids_mol()
        # Calcula y guarda los features como atributos
        atomic_nodes_features, atomic_edge_features = self.peptide_builder.builder_atomic_features(self.amino_acids_mol)
        self.library_atomic_nodes_features = atomic_nodes_features
        self.library_atomic_edge_features = atomic_edge_features
    
    def _init_amino_acid_dictionary(self, custom_dict_path):
        """ Initializes the amino acid dictionary from a specified path or defaults to a standard file.
        Args:
            custom_dict_path (str, optional): Path to the amino acid dictionary Excel file.
        """
        aa_dict_instance = AminoAcidDictionary(custom_dict_path) if custom_dict_path else AminoAcidDictionary()
        self.dictionary = aa_dict_instance.as_dict()
        self.amino_acids = aa_dict_instance.get_all_amino_acids_notations()
        
    def _build_amino_acids_mol(self):
        """ Constructs RDKit molecule objects for each amino acid in the dictionary.
        Returns:
            list: List of RDKit molecule objects for each amino acid.
        """
        amino_acids_mol = []
        for aa in self.amino_acids:
            mol = self.peptide_builder.builder_peptide(
                                                        sequence=aa,
                                                        show_display=False,
                                                        amino_acid_library=self.dictionary
                                                    )
            amino_acids_mol.append(mol)
        return amino_acids_mol
    
    def get_library_atomic_nodes_features(self):
        return self.library_atomic_nodes_features
    
    def get_library_atomic_edge_features(self):
        return self.library_atomic_edge_features
    
    def get_description_library_atomic_features(self):
        """
        Prints the length (dimension) of the atomic node and edge feature vectors.
        """
        node_features = self.get_library_atomic_nodes_features()
        edge_features = self.get_library_atomic_edge_features()
        node_dim = len(next(iter(node_features.values())))
        edge_dim = len(next(iter(edge_features.values())))
        print(f"Node feature vector length: {node_dim}")
        print(f"Edge feature vector length: {edge_dim}")
    
    
    def get_peptide_tensor_atomic_features(self, mol,device ='cpu'):
        """
        Constructs atomic feature tensors for a given peptide molecule.
        This method extracts atomic features from the molecule and converts them into PyTorch tensors.
        It uses the peptide builder to obtain the atomic features and then converts them into tensors.
        This is useful for preparing molecular data for machine learning models or graph neural networks.

        Args:
            mol (rdkit.Chem.Mol): RDKit molecule object.
            device (str, optional): Device to place the tensors on. Default is 'cpu'.

        Returns:
            tuple: A tuple containing two tensors:
                - atom_features_tensors: Tensor of shape [num_atoms, feature_dim] for atomic features.
                - edge_features_tensors: Tensor of shape [num_edges, feature_dim] for bond features.
        """
        node_ft_dict = self.library_atomic_nodes_features
        edge_ft_dict= self.library_atomic_edge_features
        node_keys_features, edge_key_features = self.peptide_builder.builder_peptide_atomic_features(mol)
        
        atom_features_tensors, edge_features_tensors = self.peptide_builder.builder_atomic_features_tensors(node_keys_features, edge_key_features, node_ft_dict, edge_ft_dict, device)
        
        return atom_features_tensors, edge_features_tensors
    
    def describe_peptide_atomic_features(self, mol, device='cpu'):
        """
        Prints the shape and a preview of the atomic and bond feature tensors for a given molecule.
        
        Args:
            mol (rdkit.Chem.Mol): RDKit molecule object.
            device (str, optional): Device to place the tensors on. Default is 'cpu'.
            preview_rows (int, optional): Number of rows to preview from each tensor.
        """
        atom_features, bond_features = self.get_peptide_tensor_atomic_features(mol, device=device)
        
        print("Atom Features:")
        print(f"  - Shape: {atom_features.shape} (atoms × feature dimensions)")
        print(f"  - Number of atoms: {atom_features.shape[0]}")
        print(f"  - Feature dimension: {atom_features.shape[1]}")
        print("\nBond Features:")
        print(f"  - Shape: {bond_features.shape} (bonds × feature dimensions)")
        print(f"  - Number of bonds: {bond_features.shape[0]}")
        print(f"  - Feature dimension: {bond_features.shape[1]}")
    
    
    def get_peptide(self, sequence, plot_peptide=False):
        """
        Build a peptide molecule from a sequence.
        
        Args:
            sequence (str): Peptide sequence.
            show_display (bool): Whether to visualize the molecule.
        
        Returns:
            Chem.Mol: RDKit molecule.
        """
        return self.peptide_builder.builder_peptide(sequence, plot_peptide, amino_acid_library=self.dictionary)
    
    def get_plot_aminoacids(self, mol, highlight_bonds=False, obtain_amino_acids=False):
        """
        Visualize or fragment the peptide molecule.\
        
        Args:
            mol (Chem.Mol): RDKit molecule.
            highlight_bonds (bool): Highlight peptide bonds.
            obtain_amino_acids (bool): Display amino acid fragments.
        """
        return self.amino_acid_builder.builder_amino_acid_plotting(mol, highlight_bonds, obtain_amino_acids)
    
    def get_amino_acid_atom_mapping(self, mol):
        """
        Get a tensor mapping each atom to its amino acid fragment.
        
        Args:
            mol (Chem.Mol): RDKit molecule.
        
        Returns:
            torch.Tensor: Mapping vector.
        """
        return self.amino_acid_builder.builder_amino_acid_mapping_vector(mol, None)
    
    def get_amino_acid_adjacency_matrix(self, sequence, device='cpu', architecture='linear'):
        """
        Get an adjacency matrix for the peptide sequence.
        
        Args:
            sequence (str): Peptide sequence.
            device (torch.device): Device for the resulting tensor.
            architecture (str): Type of adjacency matrix ('linear' or 'cyclic').
        
        Returns:
            torch.Tensor: Adjacency matrix.
        """
        return self.amino_acid_builder.builder_amino_acid_adjacency_matrix(sequence, device, architecture)



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

smiles2peptides = Smiles2Peptide()

for pep in peptides:
    print(pep)
    #aqui se instancia la clase para manejar peptidos, se puede usar ademas un diccionario personalizado
    
    mol = smiles2peptides.get_peptide(pep, plot_peptide=True)
    nodes_features, edges_features = smiles2peptides.get_peptide_tensor_atomic_features(mol, device="cuda:0")
    smiles2peptides.get_plot_aminoacids(mol, highlight_bonds=True, obtain_amino_acids=True)
    
    #print("Amino Acid mapping Vector:", smiles2peptides.get_amino_acid_atom_mapping(mol))
    #print("\nAmino Acid adjacecy Vector:", smiles2peptides.get_amino_acid_adjacency_matrix(pep, device="cuda:0", architecture='linear'))
    #print('Nodes: ',nodes_features)
    #print('Edges:',edges_features)
    smiles2peptides.describe_peptide_atomic_features(mol, device="cuda:0")
    print('----------------------\n')



smiles2peptides.get_description_library_atomic_features()

# %%
