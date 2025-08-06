#%%
from smiles2peptides.core.peptide_builder import PeptideBuilder
from smiles2peptides.core.amino_acid_builder import AminoAcidBuilder
from smiles2peptides.core.library import AminoAcidDictionary
from smiles2peptides.utils.peptide_utils import PeptideUtils


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
        """ 
        Constructs RDKit molecule objects for each amino acid in the dictionary.
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
        """
        Returns the atomic node features from the library.
        This method retrieves the atomic node features that were computed from the amino acid molecules.
        The features are stored in a dictionary where keys are feature identifiers and values are feature vectors.
        Returns:
            dict: A dictionary containing atomic node features, where keys are feature identifiers and values are feature vectors.  
        """
        return self.library_atomic_nodes_features
    
    def get_library_atomic_edge_features(self):
        """ 
        Returns the atomic edge features from the library.
        This method retrieves the atomic edge features that were computed from the amino acid molecules.
        The features are stored in a dictionary where keys are feature identifiers and values are feature vectors.
        Returns:
            dict: A dictionary containing atomic edge features, where keys are feature identifiers and values are feature vectors.      
        """
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
    
    
    def get_peptide_tensor_atomic_features(self, mol, device ='cpu'):
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
        
        atom_features_tensors, edge_features_tensors = self.peptide_builder.builder_atomic_features_tensors(
                                                                                                            node_keys_features,
                                                                                                            edge_key_features,
                                                                                                            node_ft_dict,
                                                                                                            edge_ft_dict,
                                                                                                            device
                                                                                                            )
        
        return atom_features_tensors, edge_features_tensors
    
    def describe_peptide_atomic_features(self, mol, device='cpu'):
        """
        Prints the shape and a preview of the atomic and bond feature tensors for a given molecule.
        
        Args:
            mol (rdkit.Chem.Mol): RDKit molecule object.
            device (str, optional): Device to place the tensors on. Default is 'cpu'. Default is 'cpu'.
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
    
    def get_peptide_atomic_adjacency_matrix(self, mol, device='cpu'):
        """
        Constructs an adjacency matrix for the atoms in a peptide molecule.
        
        Args:
            mol (rdkit.Chem.Mol): RDKit molecule object.
            device (str or torch.device): Device to place the tensor on ('cpu' or 'cuda'). Default is 'cpu'.
        
        Returns:
            torch.Tensor: Adjacency matrix of shape [num_atoms, num_atoms].
        """
        return self.peptide_builder.builder_atomic_adjacency_matrix(mol, device)
    
    def get_plot_aminoacids(self, mol, highlight_bonds=False, obtain_amino_acids=False):
        """
        Visualize or fragment the peptide molecule.\
        
        Args:
            mol (Chem.Mol): RDKit molecule.
            highlight_bonds (bool): Highlight peptide bonds.
            obtain_amino_acids (bool): Display amino acid fragments.
        """
        return self.amino_acid_builder.builder_amino_acid_plotting(mol, highlight_bonds, obtain_amino_acids)
    
    def get_amino_acid_atom_mapping(self, mol, device='cpu'):
        """
        Get a tensor mapping each atom to its amino acid fragment.
        
        Args:
            mol (Chem.Mol): RDKit molecule.
        
        Returns:
            torch.Tensor: Mapping vector.
        """
        return self.amino_acid_builder.builder_amino_acid_mapping_vector(mol, device)
    
    def get_amino_acid_adjacency_matrix(self, sequence, device, architecture):
        """
        Get an adjacency matrix for the peptide sequence.
        
        Args:
            sequence (str): Peptide sequence.
            device (torch.device): Device for the resulting tensor. Default is 'cpu'.
            architecture (str): Type of adjacency matrix ('linear' or 'cyclic').
        
        Returns:
            torch.Tensor: Adjacency matrix.
        """
        return self.amino_acid_builder.builder_amino_acid_adjacency_matrix(sequence, device, architecture)
    
    def get_amino_acid_features(self, sequence, device='cpu'):
        """
        Constructs a tensor representing the features of each amino acid in a peptide sequence.
        Args:
            sequence (str): Peptide sequence string.
            device (str or torch.device): Device to store the resulting tensor (e.g., 'cpu' or 'cuda'). Default is 'cpu'.
        """
        return self.amino_acid_builder.builder_amino_acid_features(sequence, self.dictionary, device)


""" peptides = [
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
] """

peptides = [
    "ARNDC", "QEGHI", "LKMFPS", "TWYVAG", "arndc", "qeGhi", "lkmfps", "twyvG",  # L y D individuales
    "AGrdkV", "FPlheT", "GcNqYr", "kMWtvS", "YhGRaC", "LVnpQA", "iGpeKD", "sRqGMF",
    "GfTclWr", "aKGHlp", "MGryta", "wNEKGi", "dvmGqTH", "QRPGvya", "hSLcwmK", "GiFNrya", 
    "eGlVqRM", "aGHkpTV", "CGLsnqA", "ytmwGRP", "sGHvQeM", "KGnywva", "HFctGMp", "GrkLhNa",
    "aGKhdPL", "RMGyVnt", "AQlecvF", "wsiKNrG", "pAGTYwL", "kcGHvmD", "aWErPGt", "qMhGLra",
    "rGnCAyL", "DVKGlta", "sykNpGA", "mQGhRtv", "FhpKGwG", "VnGtkLC", "aiGHydp", "rVLaTGM",
    "GEKhFML", "tGRsaYv"
]


smiles2peptides = Smiles2Peptide()

for pep in peptides:
    print(pep)
    #aqui se instancia la clase para manejar peptidos, se puede usar ademas un diccionario personalizado
    
    mol = smiles2peptides.get_peptide(pep, plot_peptide=True)
    nodes_features, edges_features = smiles2peptides.get_peptide_tensor_atomic_features(mol, device="cuda:0")
    smiles2peptides.get_plot_aminoacids(mol, highlight_bonds=True, obtain_amino_acids=True)
    
    print("Amino Acid mapping Vector:", smiles2peptides.get_amino_acid_atom_mapping(mol))
    print("\nAmino Acid adjacecy Vector:", smiles2peptides.get_amino_acid_adjacency_matrix(pep, device="cuda:0", architecture='linear'))
    print('Nodes: ',nodes_features)
    print('Edges:',edges_features)
    smiles2peptides.describe_peptide_atomic_features(mol, device="cuda:0")
    adjacency_matrix = smiles2peptides.get_peptide_atomic_adjacency_matrix(mol, device="cuda:0")
    print(f"Adjacency Matrix Shape: {adjacency_matrix}")
    print("Amino Acid Features Tensor:", smiles2peptides.get_amino_acid_features(pep, device="cuda:0"))
    print('----------------------\n')



smiles2peptides.get_description_library_atomic_features()
# %%
