#%%
from smiles2peptides.interface import Smiles2Peptide

smiles2pep = Smiles2Peptide()

""" peptides = [
    "{pra}{FITC-Ahx}{Pra}nT{am}",
    "{PEG2}pr{GlcNAc-T}l{am}",
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
]  """

peptides = [
    "ARNDC", "QEGHI", "LKMFPS", "TWYVAG", "arndc", "qeGhi", "lkmfps", "twyvG",  # L y D individuales
    "AGrdkV", "FPlheT", "GcNqYr", "kMWtvS", "YhGRaC", "LVnpQA", "iGpeKD", "sRqGMF",
    "GfTclWr", "aKGHlp", "MGryta", "wNEKGi{am}", "dvmGqTH", "QRPGvya", "hSLcwmK", "GiFNrya", 
    "eGlVqRM", "aGHkpTV", "CGLsnqA", "ytmwGRP", "sGHvQeM", "KGnywva", "HFctGMp", "GrkLhNa",
    "aGKhdPL", "RMGyVnt", "AQlecvF", "wsiKNrG", "pAGTYwL", "kcGHvmD", "aWErPGt", "qMhGLra",
    "rGnCAyL", "DVKGlta", "{ac}sykNpGA", "mQGhRtv", "FhpKGwG", "VnGtkLC", "aiGHydp", "rVLaTGM",
    "GEKhFML", "tGRsaYv"
]


smiles2peptides = Smiles2Peptide()

for pep in peptides:
    print(pep)
    #aqui se instancia la clase para manejar peptidos, se puede usar ademas un diccionario personalizado
    
    mol = smiles2peptides.get_peptide(pep, plot_peptide=True)
    nodes_features, edges_features = smiles2peptides.get_peptide_tensor_atomic_features(mol, device="cpu")
    smiles2peptides.get_plot_aminoacids(mol, highlight_bonds=True, obtain_amino_acids=True)
    
    print("Amino Acid mapping Vector:", smiles2peptides.get_amino_acid_atom_mapping(mol))
    print("\nAmino Acid adjacecy Vector:", smiles2peptides.get_amino_acid_adjacency_matrix(pep, device="cpu", architecture='linear'))
    print('Nodes: ',nodes_features)
    print('Edges:',edges_features)
    smiles2peptides.describe_peptide_atomic_features(mol, device="cpu")
    adjacency_matrix = smiles2peptides.get_peptide_atomic_adjacency_matrix(mol, device="cpu")
    print(f"Adjacency Matrix Shape: {adjacency_matrix}")
    print("Amino Acid Features Tensor:", smiles2peptides.get_amino_acid_features(pep, device="cpu"))
    print('----------------------\n')



smiles2peptides.get_description_library_atomic_features()
# %%
