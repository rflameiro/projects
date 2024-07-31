import numpy as np
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from chembl_structure_pipeline import standardizer


def new_standardize_smiles(smiles, return_mol=True):
    """
    Input = SMILES
    return_mol: Return RDKit Mol object (True) or SMILES (False)
    """
    try:
        molecule = Chem.MolFromSmiles(smiles)
    except:
        return
    
    try:
        # Apply standardize_mol
        std1_mol = standardizer.standardize_mol(molecule)
        # Select largest fragment
        desalter = rdMolStandardize.LargestFragmentChooser()
        desalt_mol = desalter.choose(std1_mol)
        # Apply standardize_mol to the largest fragment
        std2_mol = standardizer.standardize_mol(desalt_mol)
        # Select tautomer with highest score
        te = rdMolStandardize.TautomerEnumerator()
        taut_uncharged_parent_clean_mol = te.Canonicalize(std2_mol)
        # Remove stereochemistry information
        taut_smi = Chem.MolToSmiles(taut_uncharged_parent_clean_mol, isomericSmiles=False)  
        # This can still fail conversion to RDKit Mol, therefore, the following is necessary        
        final_smi = Chem.MolToSmiles(Chem.MolFromSmiles(taut_smi))
        
        if not final_smi:
            return
        
        if return_mol:
            return Chem.MolFromSmiles(final_smi)
        else:
            return final_smi
    
    except:
        return


# Check for compounds with non-organic atoms
not_organic_pat = Chem.MolFromSmarts("[!#5;!#6;!#7;!#8;!#16;!#15;!F;!Cl;!Br;!I;!#1]")

def non_organic(smiles):
    # Returns True if a non-organic atom is found 
    # or SMILES failed conversion to RDKit Mol
    try:
        mol = Chem.MolFromSmiles(smiles)
        return bool(mol.GetSubstructMatch(not_organic_pat))
    except:  # failed to convert to rdkit_mol
        return True