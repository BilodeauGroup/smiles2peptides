import pandas as pd
from rdkit import Chem  

from smiles2peptides.project.utils.peptide_utils import PeptideUtils


class PeptideBuilder:
    """
    Class responsible for constructing a peptide molecule from a given sequence.
    Optionally displays a visual representation of the resulting molecule.
    """   
    @staticmethod
    def builder_peptide(sequence, show_display, amino_acid_library):
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
                if character == "{am}":
                    smile = amino_acid_library[character][1]
                    concatenated_smile += smile
                #TODO lógica pendiente para otros casos
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
    
    @staticmethod
    def builder_atomic_adjacency_matrix(mol, device):
        return PeptideUtils.util_atomic_adjacency_matrix(mol, device)

