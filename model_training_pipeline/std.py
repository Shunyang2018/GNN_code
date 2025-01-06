from rdkit import Chem
from molvs import tautomer
from rdkit.Chem.Draw import IPythonConsole #Needed to show molecules
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions #Only needed if modifying defaults
DrawingOptions.bondLineWidth=1.8
DrawingOptions.atomLabelFontSize = 8
def _InitialiseNeutralisationReactions():
    """
    Initialize a list of SMARTS patterns and their corresponding SMILES transformations
    for neutralizing charges in molecules.
    
    Returns:
        list of tuples: Each tuple contains a SMARTS pattern and its corresponding SMILES.
    """    
    patts= (
        # Imidazoles
        ('[n+;H]','n'),
        # Amines
        ('[N+;!H0]','N'),
        # Carboxylic acids and alcohols
        ('[$([O-]);!$([O-][#7])]','O'),
        # Thiols
        ('[S-;X1]','S'),
        # Sulfonamides
        ('[$([N-;X2]S(=O)=O)]','N'),
        # Enamines
        ('[$([N-;X2][C,N]=C)]','N'),
        # Tetrazoles
        ('[n-]','[nH]'),
        # Sulfoxides
        ('[$([S-]=O)]','S'),
        # Amides
        ('[$([N-]C=O)]','N'),
        )
    return [(Chem.MolFromSmarts(x),Chem.MolFromSmiles(y)) for x,y in patts]

_reactions=None
def neutralise_charges(mol, reactions=None):
    """
    Neutralize charges in a molecule using predefined SMARTS patterns.

    Args:
        mol (Mol): RDKit molecule object.
        reactions (list, optional): List of reactions for neutralisation. Defaults to None.

    Returns:
        tuple: (neutralized Mol, boolean indicating if changes were made).
    """
    global _reactions
    if reactions is None:
        if _reactions is None:
            _reactions=_InitialiseNeutralisationReactions()
        reactions=_reactions
    replaced = False
    for i,(reactant, product) in enumerate(reactions):
        while mol.HasSubstructMatch(reactant):
            replaced = True
            rms = Chem.AllChem.ReplaceSubstructs(mol, reactant, product, replaceAll=True)
            mol = rms[0]
    if replaced:
        Chem.SanitizeMol(mol)
        Chem.rdmolops.Kekulize(mol, clearAromaticFlags=True)
        return (mol, True)
    else:
        return (mol, False)

def desalt(mol):
    """
    Remove salts and keep the largest fragment of a molecule.

    Args:
        mol (Mol): RDKit molecule object.

    Returns:
        tuple: (Mol, boolean indicating if desalt was performed).
    """
    d = Chem.rdmolops.GetMolFrags(mol) #these are atom indices
    if len(d) == 1: #If there are fragments or multiple molecules this will be greater than 1
        return mol,False
    my_smiles=Chem.MolToSmiles(mol,True)
    parent_atom_count=0;
    disconnected=my_smiles.split('.')
    #With GetMolFrags, we've already established that there is more than one disconnected structure
    status = False
    for s in disconnected:
        little_mol = Chem.MolFromSmiles(s,sanitize=True)
        #Sanitize=True will fail for choline sulfate.  Can't sanitize the radical.
        if little_mol is not None:
            count = little_mol.GetNumAtoms()
            if count > parent_atom_count:
                parent_atom_count = count
                parent_mol = little_mol
                status = True
    return parent_mol,status

# Given an rdkit mol object, it returns a (rdkit mol object, smiles, and InChI) tuple in standardized form
def mol_to_neutral_desalted_canonical(rdkit_mol):
    """
    Standardize a molecule by neutralizing charges, removing salts, and generating canonical forms.

    Args:
        rdkit_mol (Mol): RDKit molecule object.

    Returns:
        tuple: (Mol, SMILES, InChI) of the standardized molecule, or (False, False, False) on failure.
    """
    canon = tautomer.TautomerCanonicalizer()
    if rdkit_mol is not None:
        mol = rdkit_mol
        if (mol.GetNumAtoms()>0):
            Chem.SanitizeMol(mol)
            Chem.rdmolops.Kekulize(mol, clearAromaticFlags=True)
            mol, status = neutralise_charges(mol)
            mol, status = desalt(mol)
            try:
                mol = standardizer.standardize(mol)
                Chem.SanitizeMol(mol)
                Chem.rdmolops.Kekulize(mol, clearAromaticFlags=True)
            except:
                Chem.SanitizeMol(mol)
                Chem.rdmolops.Kekulize(mol, clearAromaticFlags=True)
            mol = canon.canonicalize(mol)
            Chem.SanitizeMol(mol)
            Chem.rdmolops.Kekulize(mol, clearAromaticFlags=True)
            smiles = Chem.MolToSmiles(mol,isomericSmiles=True)
            inchi = Chem.MolToInchi(mol)
            return (mol, smiles, inchi)
    else:
        return (False, False, False)
def stdsmi(smiles):
    """
    Standardize a SMILES string by converting to a canonical, desalted, and neutralized form.

    Args:
        smiles (str): Input SMILES string.

    Returns:
        str: Standardized SMILES string.
    """
    standard_structure_tuple = mol_to_neutral_desalted_canonical(Chem.MolFromSmiles(smiles))
    #print('Original SMILES: \n' + smiles)
    #print('Standardized SMILES: \n' + standard_structure_tuple[1])
    return standard_structure_tuple[1]