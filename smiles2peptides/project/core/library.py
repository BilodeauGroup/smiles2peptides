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
            dictionary_path_file = os.path.join(base_dir, "amino_acids_library.xlsx")
        
        self.excel_path = dictionary_path_file
        self._dictionary = self._load_dictionary()
    
    def _load_dictionary(self):
        """
        Reads the Excel file and builds an internal dictionary.
        
        Returns:
            dict: Mapping from amino acid notation (str) to a tuple containing (Chirality (str), SMILES (str), Exception (str), Type (str)).
        """
        df = pd.read_excel(self.excel_path, engine='openpyxl')
        
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


        