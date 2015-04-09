#Series of functions to parse molecules based on RDKit
import sys
from rdkit import Chem
from sanifix import fix_mol


def parse_mol_simple(my_type, txt):
    """Function to parse individual mols given a type"""
    if my_type == "mol":
        # Try this way
        mol = Chem.MolFromMolBlock(txt.strip())
        if mol is None:
            mol = Chem.MolFromMolBlock(txt)
        if mol is None:
            mol = Chem.MolFromMolBlock("\n".join(txt.split("\n")[1:]))
        # Now try to do sanidfix
        if mol is None:
            mol = fix_mol(Chem.MolFromMolBlock(txt, False))
        # Annd again
        if mol is None:
            mol = fix_mol(Chem.MolFromMolBlock(txt.strip(), False))
    elif my_type == "smiles":
        # Assumes that smiles is the first column -> and splits on chemaxon
        mol = Chem.MolFromSmiles(txt.split()[0].split(":")[0])
    elif my_type == "inchi":
        # Assumes that INCHI is the first column
        mol = Chem.MolFromInchi(my_txt.split()[0], my_vals)
    if mol is None:
        print txt
    return mol


def parse_mol_json(molobj):
    """Function to get the RDKit mol for a java mol obj"""
    molstr = molobj["source"]
    # Get the format and use this as a starting poitn to work out 
    molformat = molobj["format"]
    # Now parse it with RDKit
    return parse_mol_simple(molformat, molstr)

