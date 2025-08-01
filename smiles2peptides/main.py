#%%
import os
import pandas as pd
from rdkit import Chem  
from rdkit.Chem import Draw  
from IPython.display import display 
import re

# Get the path of the amino acid library Excel file
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
EXCEL_PATH = os.path.join(BASE_DIR, "amino_acids_v2.xlsx")


def create_dictionary_from_excel():
    """
    Reads an Excel file and creates a dictionary where the keys are the values
    from the 'Notation' column and the values are tuples with 'Chirality' and 'SMILE'.
    
    Returns:
        dict: Dictionary with the structure {Amino Acid: (Chirality, SMILE)}.
    """
    df = pd.read_excel(EXCEL_PATH)  # Usa la ruta dentro del paquete
    dictionary = {
        str(row['Notation']).strip():
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
        
        # Regex para los posibles inicios
        n_pattern = r"^\[(NH2|NH1|N|NH3):1\]"
        
        if index == 0:
            # Check next character
            if index + 1 < len(characters) and characters[index + 1] in dictionary:
                next_smile = dictionary[characters[index + 1]][1]
                match = re.match(n_pattern, next_smile) #it find the N-terminus
                if match:
                    prefix = match.group(0) #Select what N-terminus is present
                    rest = next_smile[len(prefix):] #it gets the rest of SMILE after the N-terminus string
                    next_smile = prefix + special_smile + rest  #Create the new a.a SMILE with the modificaction
                    next_smile = remove_o_if_needed(next_smile, index + 2, characters) #Remove the H from the N-Terminus
                    return next_smile, index + 1  # Skip next character
                else:
                    raise ValueError(f"The SMILE of the character '{characters[index + 1]}' does not start with '[NH2:1]', '[NH1:1]', '[N:1]' o '[NH3:1]'.")
            else:
                raise ValueError(f"There is no valid character after '{character}'.")
        else:
            # Handle non-first character
            prev_smile = dictionary[characters[index - 1]][1]
            if re.match(n_pattern, prev_smile):
                prev_smile = prev_smile[:1] + special_smile + prev_smile[1:]
                prev_smile = remove_o_if_needed(prev_smile, index, characters)
                return prev_smile, index
            else:
                raise ValueError(f"The SMILE of the character '{characters[index - 1]}' does not start with '[NH2:1]', '[NH1:1]', '[N:1]' o '[NH3:1]'.")
    
    def remove_o_if_needed(smile, index, characters, character=None): # character is optional
        """ Helper function to remove 'O' if the next character is not the last one.
            Also performs nitrogen group replacements and warns if [N:1] is present.
        """
        # Remove "O" if needed
        if index < len(characters) and smile.endswith("O"):
            smile = smile[:-1]
        # Nitrogen group replacements
        if index>1:
            if smile.startswith("[NH2:1]"):
                smile = smile.replace("[NH2:1]", "[NH1:1]", 1)
            elif smile.startswith("[NH1:1]"):
                smile = smile.replace("[NH1:1]", "[N:1]", 1)
            elif smile.startswith("[NH3:1]"):
                smile = smile.replace("[NH3:1]", "[NH2:1]", 1)
            # Warning if starts with [N:1]
            if smile.startswith("[N:1]"):
                if character is not None:
                    print(f"Warning: The amino acid '{character}' cannot form a peptide bond because the nitrogen's bonds are already saturated.")
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
                smile = remove_o_if_needed(smile, i+1, characters)
                concatenated_smile += str(smile).strip()
            else:
                raise ValueError(f"The amino acid '{character}' is not found in the library.")
        
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
        # Usamos el dibujante 2D para poder activar indices
        drawer = Draw.MolDraw2DCairo(800, 800)
        options = drawer.drawOptions()
        options.addAtomIndices = False # Muestra los indices de los átomos
        drawer.DrawMolecule(rdkit_mol)
        drawer.FinishDrawing()
        img = drawer.GetDrawingText()
        
        # Mostrar en Jupyter
        from IPython.display import Image, display
        display(Image(data=img))
    
    # Return the generated RDKit molecule
    return rdkit_mol



peptides = [
    "{PEG2}pr{GlcNAc-T}l{am}",
    "{6-Br-W}{4-hydrox-P}G{me-f}{Pra}nT",
    "v{photo-L}{phospho-Y}l{4-tert-butyl-p}iK",
    "{acm-C}wM{4&5-hydrox-L}{(N-me)-a}{iso-Q}",
    "{4-hydrox-p}{seleno-C}R{trime-L}{3&4-dihydrox-F}G",
    "K{p-carboxyl-F}{photo-L}c{phospho-Y}{me-f}",
    "{nor-r}{photo-M}l{GlcNAc-S}H{me-f}",
    "i{h-s}{(N-me)-A}{Pra}w{acm-C}l",
    "{C-me}d{p-carboxyl-f}{photo-L}{me-Y}K",
    "{photo-l}{p-carboxyl-f}m{3-me-f}w{4-hydrox-p}",
    "a{Pra}{trime-L}y{iso-Q}{photo-M}",
    "{(N-me)-a}H{3-me-f}{photo-M}l{p-carboxyl-F}",
    "M{phospho-Y}{4-hydrox-P}d{acm-c}{Pra}",
    "n{4-hydrox-p}{GlcNAc-S}K{(N-me)-A}{4-tert-butyl-p}",
    "{succinyl-K}{photo-M}{4-hydrox-p}R{me-f}t",
    "G{me-F}{Dip-a}{photo-M}{p-carboxyl-f}{4-hydrox-p}",
    "{seleno-C}v{Pra}{p-carboxyl-F}{GlcNAc-T}w",
    "p{photo-L}{me-f}{trime-L}K{3&4-dihydrox-F}",
    "{4-hydrox-P}y{iso-Q}{photo-M}{h-S}F",
    "W{p-carboxyl-f}M{photo-M}q{C-me}",
    "{phospho-S}{4-hydrox-p}l{succinyl-K}{3-me-f}V",
    "{photo-l}{seleno-C}{(N-me)-A}{4-hydrox-P}dF",
    "c{me-y}{photo-M}{3&4-dihydrox-F}{acm-C}t",
    "{GlcNAc-asn}{photo-M}n{h-S}{4-tert-butyl-p}l",
    "{iso-D}p{Pra}{me-f}H{p-carboxyl-F}",
    "{acm-c}{4-hydrox-P}G{(N-me)-A}R{photo-M}",
    "{4-hydrox-p}{photo-L}{succinyl-K}y{3-me-f}t",
    "L{me-f}{photo-M}{4&5-hydrox-L}m{iso-Q}",
    "{trime-L}H{4-hydrox-p}{photo-l}{C-me}d",
    "{Pra}{seleno-C}{4-hydrox-p}n{phospho-Y}L",
    "t{acm-C}{photo-M}{me-f}{4-tert-butyl-p}k",
    "{4-hydrox-p}V{3-me-f}{iso-Q}{photo-l}y",
    "{4-hydrox-P}{photo-M}{C-me}H{acm-c}m",
    "{succinyl-K}{photo-M}{p-carboxyl-f}{h-s}T",
    "F{3-me-f}{4-hydrox-p}{me-f}{Pra}w",
    "{photo-l}n{4&5-hydrox-L}M{trime-L}y",
    "{photo-M}{me-f}{phospho-Y}{4-hydrox-p}G{C-me}",
    "L{p-carboxyl-F}{photo-L}{Pra}h{iso-Q}",
    "{(N-me)-a}{seleno-C}{photo-M}R{4-hydrox-p}w",
    "{phospho-S}{photo-l}{4-hydrox-P}{3-me-f}yM",
    "{photo-L}F{h-s}{GlcNAc-T}{Pra}c",
    "{trime-L}{4-hydrox-p}y{photo-M}{me-f}k",
    "R{(N-me)-A}{4-hydrox-p}{Pra}{4-tert-butyl-p}l",
    "a{photo-M}{me-F}{4-hydrox-p}{succinyl-K}n",
    "K{4-hydrox-p}{photo-l}q{(N-me)-a}{iso-Q}",
    "{4-hydrox-P}{Pra}{acm-C}{photo-M}Lw",
    "{photo-l}{C-me}G{p-carboxyl-f}{me-f}d",
    "{photo-M}{4-hydrox-p}{GlcNAc-asn}R{3-me-f}v",
    "{me-y}{photo-L}{4-hydrox-P}{succinyl-K}F",
    "{p-carboxyl-f}t{photo-M}K{h-s}{(N-me)-A}"
]

#TODO INCLUIR LAS RESTRICCIONES EN LOS AMINO ACIDOS QUE SOLO PUEDEN ESTAR EN POSICIONES INICIALE O TERMINALES, DEBE SER UNA VERIFACION

for pep in peptides:
    print('----------------------')
    print(pep)
    mol = generating_rdkit_mol(sequence=pep, show_display=False)
    # Obtener átomos con map number 1 y 2
    atoms_map1 = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomMapNum() == 1]
    atoms_map2 = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomMapNum() == 2]
    # Buscar enlaces entre átomos map1 y map2
    bonds_to_highlight = []
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        if (a1 in atoms_map1 and a2 in atoms_map2) or (a1 in atoms_map2 and a2 in atoms_map1):
            bonds_to_highlight.append(bond.GetIdx())
    # Colores opcionales para resaltar (rojo brillante)
    highlight_bond_colors = {idx: (0.0, 1.0, 1.0) for idx in bonds_to_highlight}
    # Dibujar la molécula resaltando los enlaces
    img = Draw.MolToImage(
                        mol,
                        size=(600, 600),
                        highlightBonds=bonds_to_highlight,
                        highlightBondColors=highlight_bond_colors
                        )
    display(img)


# %%
