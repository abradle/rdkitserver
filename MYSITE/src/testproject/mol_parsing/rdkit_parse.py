#Series of functions to parse molecules based on RDKit
import sys
from rdkit import Chem

def parse_mol_simple(my_type, txt):
    """Function to parse individual mols given a type"""
    if my_type == "mol":
        # Try this way
        mol = Chem.MolFromMolBlock(txt.strip())
        if mol is None:
            mol = Chem.MolFromMolBlock(txt)
        if mol is None:
            mol = Chem.MolFromMolBlock("\n".join(txt.split("\n")[1:]))
        if mol is None:
            mol = Chem.MolFromMolBlock(txt, strictParsing=False)
        if mol is None:
            mol = Chem.MolFromMolBlock(txt.strip(), strictParsing=False)
    elif my_type == "smiles":
        # Assumes that smiles is the first column
        mol = Chem.MolFromSmiles(txt.split()[0])
    elif my_type == "inchi":
        # Assumes that INCHI is the first column
        mol = Chem.MolFromInchi(my_txt.split()[0], my_vals)
    return mol


def parse_mol_json(molobj):
    """Function to get the RDKit mol for a java mol obj"""
    molstr = molobj["source"]
    # Get the format and use this as a starting poitn to work out 
    molformat = molobj["format"]
    # Now parse it with RDKit
    return parse_mol_simple(molformat, molstr)

