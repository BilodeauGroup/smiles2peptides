#%%
import os
import pandas as pd
from rdkit import Chem  
from rdkit.Chem import Draw  
from IPython.display import display 
import re

# Get the path of the amino acid library Excel file
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
EXCEL_PATH = os.path.join(BASE_DIR, "amino_acid_library.xlsx")


def create_dictionary_from_excel():
    """
    Reads an Excel file and creates a dictionary where the keys are the values
    from the 'Amino Acid' column and the values are tuples with 'Chirality' and 'SMILE'.
    
    Returns:
        dict: Dictionary with the structure {Amino Acid: (Chirality, SMILE)}.
    """
    df = pd.read_excel(EXCEL_PATH)  # Usa la ruta dentro del paquete
    dictionary = {
        str(row['Amino Acid']).strip():
            (str(row['Chirality']).strip(), str(row['SMILES']).strip()) 
        for _, row in df.iterrows()
    }
    return dictionary


#Proccesing Sequence Section

def extract_characters(sequence):
    """
    Extracts individual characters from a sequence, considering anything inside
    curly braces `{}` as a single character.
    
    Args:
        sequence (str): The input sequence.
        
    Returns:
        list: A list of extracted characters, including groups inside curly braces as units.
    """
    # Regular expression to find groups inside `{}` or individual characters
    pattern = r'\{[^}]*\}|[a-zA-Z]'
    # Find all matches
    characters = re.findall(pattern, sequence)
    return characters


def process_sequence(sequence, dictionary):
    
    """
    Processes a sequence to obtain its concatenated SMILE based on the dictionary.
    
    Args:
        sequence (str): Input sequence.
        dictionary (dict): Dictionary with {Amino Acid: (Chirality, SMILE)}.
        
    Returns:
        str: Concatenated SMILES for the sequence.
    """
    sequence = sequence.replace('\u200b', '').replace(' ', '')
    
    def handle_special_case(character, index, characters, special_smile):
        """ Helper function to handle {PEG} or {ac} cases. """
        if character not in dictionary:
            raise ValueError(f"The character '{character}' is not found in the dictionary.")
        
        if index == 0:
            # Check next character
            if index + 1 < len(characters) and characters[index + 1] in dictionary:
                next_smile = dictionary[characters[index + 1]][1]
                if next_smile.startswith("N"):
                    next_smile = next_smile[:1] + special_smile + next_smile[1:]
                    next_smile = remove_o_if_needed(next_smile, index + 2, characters)
                    return next_smile, index + 1  # Skip next character
                else:
                    raise ValueError(f"The SMILE of the character '{characters[index + 1]}' does not start with 'N'.")
            else:
                raise ValueError(f"There is no valid character after '{character}'.")
        else:
            # Handle non-first character
            prev_smile = dictionary[characters[index - 1]][1]
            if prev_smile.startswith("N"):
                prev_smile = prev_smile[:1] + special_smile + prev_smile[1:]
                prev_smile = remove_o_if_needed(prev_smile, index, characters)
                return prev_smile, index
            else:
                raise ValueError(f"The SMILE of the character '{characters[index - 1]}' does not start with 'N'.")
    
    def remove_o_if_needed(smile, index, characters):
        """ Helper function to remove 'O' if the next character is not the last one. """
        if index < len(characters) and smile.endswith("O"):
            smile = smile[:-1]  # Remove the "O"
        return smile
    
    # Extract characters from the sequence
    characters = extract_characters(sequence)
    
    # Initialize the concatenated SMILE
    concatenated_smile = ""  
    i = 0  # Index to iterate over the characters
    
    while i < len(characters):
        character = characters[i]
        
        if character.startswith("{PEG"):
            peg_smile = dictionary[character][1]  # Get the associated SMILE for PEG
            next_smile, i = handle_special_case(character, i, characters, peg_smile)
            concatenated_smile += next_smile
        
        elif character.startswith("{ac}"):
            ac_smile = dictionary[character][1] # Get the associated SMILE for acetylation
            next_smile, i = handle_special_case(character, i, characters, ac_smile)
            concatenated_smile += next_smile
        else:
            # Normal handling for characters that are not {PEG} or {ac}
            if character in dictionary:
                smile = dictionary[character][1] # Get the associated SMILE
                # Remove "O" if the next character is not the last one
                smile = remove_o_if_needed(smile, i + 1, characters)
                concatenated_smile += str(smile).strip()
            else:
                raise ValueError(f"The character '{character}' is not found in the dictionary.")
        
        i += 1  # Increment index after processing
    
    concatenated_smile.replace('\u200b', '').replace(' ', '')
    
    if sequence.startswith("{FITC-Ahx}"):
        editable_peptide = Chem.MolFromSmiles(concatenated_smile)
        editable_mol = Chem.EditableMol(editable_peptide) 
        editable_mol.RemoveBond(3, 4)
        editable_mol.AddBond(3, 5, order=Chem.rdchem.BondType.SINGLE)
        peptide_molecule  = editable_mol.GetMol()
    else: 
        peptide_molecule = Chem.MolFromSmiles(concatenated_smile)
    
    return  peptide_molecule

#Main Function

# Define the main function that generates the RDKit molecule from a sequence
def generating_rdkit_mol(sequence, show_display=False):
    """
    Generates an RDKit molecule from a sequence using the amino acid dictionary.
    
    Parameters:
        sequence (str): The amino acid sequence or components.
        show_display (bool): If True, displays the generated molecule image.
    
    Returns:
        Chem.Mol: An RDKit molecule representing the full sequence.
    """
    # Create an amino acid dictionary from the Excel file
    dictionary = create_dictionary_from_excel()
    
    # Process the sequence to convert it into an RDKit molecule
    rdkit_mol = process_sequence(sequence, dictionary)
    
    # If displaying is requested and the molecule was successfully generated
    if show_display and rdkit_mol:
        # Convert the RDKit molecule to an 800x800 pixel image
        img = Draw.MolToImage(rdkit_mol, size=(800, 800))
        display(img)   
    
    # Return the generated RDKit molecule
    return rdkit_mol

#%%
natural_amino_acids = set("ARNDCQEGHILKMFPSTWYVarndcqeghilkmfpstwyv")
print(natural_amino_acids)

# %%
