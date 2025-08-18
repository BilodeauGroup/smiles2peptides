from smiles2peptides.project.utils.amino_acid_utils import  AminoAcidUtils


class AminoAcidBuilder:
    """
    Class for handling peptide bond processing in RDKit molecules.
    Provides methods to identify peptide bonds, highlight them, and optionally fragment the molecule.
    """
    @staticmethod
    def builder_amino_acid_plotting(
                                    mol,
                                    highlight_bonds,
                                    obtain_amino_acids,
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
    def builder_amino_acid_mapping_vector(mol, device):
        """ Returns a tensor mapping each atom to its amino acid fragment."""
        
        peptidic_bonds = AminoAcidUtils.util_peptide_bonds(mol)
        
        return AminoAcidUtils.util_amino_acid_mapping_vector(mol, peptidic_bonds, device)
    
    
    @staticmethod
    def builder_amino_acid_adjacency_matrix(sequence, device, architecture, amino_acid_library):
        """
        Constructs an adjacency matrix for a peptide sequence, representing linear connections between amino acids.
        
        Args:
            sequence (str): Peptide sequence string.
            device (torch.device): Device to store the resulting tensor (e.g., 'cpu' or 'cuda').
            architecture (str): Type of adjacency matrix ('linear' or 'cyclic'). Defaults to 'linear'.
        
        Returns:
            torch.Tensor: Adjacency matrix representing connections between amino acids.
        """
        return AminoAcidUtils.util_amino_acid_adjacency_matrix(sequence, device, architecture, amino_acid_library)
    
    @staticmethod
    def builder_amino_acid_features(sequence, amino_acid_library, device):
        return AminoAcidUtils.util_amino_acid_features_tensor(sequence, amino_acid_library, device)
