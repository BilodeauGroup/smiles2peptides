import pandas as pd

def create_dictionary_from_excel():
    """
    Reads an Excel file and creates a dictionary where the keys are the values
    from the 'Amino Acid' column and the values are tuples with 'Chirality' and 'SMILE'.
    
    Args:
        file_path (str): Path to the Excel file.
        
    Returns:
        dict: Dictionary with the structure {Amino Acid: (Chirality, SMILE)}.
    """
    # Read the Excel file
    df = pd.read_excel('data.xlsx')
    
    # Create the dictionary
    dictionary = {
                    str(row['Amino Acid']).strip():
                                                    (str(row['Chirality']).strip(),
                                                    str(row['SMILES']).strip()) 
                                                    for _, row in df.iterrows()
                    }
    
    return dictionary