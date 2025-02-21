#%%
# Import necessary libraries
from rdkit import Chem  
from rdkit.Chem import Draw  
from IPython.display import display 
from functions.processing_seq import process_sequence 
from functions.creating_dic import create_dictionary_from_excel 

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

""" sequence = "{PEG2[CH2CH2]N3}{h-p}{Pra}{photo-M}{me-T}{nap-A}{acm-c}t{am}"  # Change this sequence as needed
rdkit_mol = generating_rdkit_mol(sequence, show_display=True)  """

# %%
