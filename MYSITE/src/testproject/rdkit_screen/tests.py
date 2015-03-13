from rdkit_screen.functions import SimMethods, LibMethods, FPMethods
from rdkit import Chem
import json
import ast
from json_function.json_parse import remove_keys 

def index(request):
    return HttpResponse("WELCOME TO INDEX")


def test_screen():
    """View to take a smiles and then screen against a known library of actives"""
    # Take the smiles in the request object
    scr_mols = [{"RDMOL": Chem.MolFromSmiles(str("CCCCC"))}]  
    # Now run the process
    # Get the library
    libm = LibMethods(ast.literal_eval(open("mols.json").read()))
    mols = libm.get_mols()
    # Get the fps
    fpm = FPMethods("morgan")
    screen_fps = fpm.get_fps(mols)
    # Get the fp for my mol(s)
    mol_fps = fpm.get_fps(scr_mols)
    simm = SimMethods("tanimoto")
    threshold = 0.1
    # Store the result in this
    out_d = {}
    for i, mol in enumerate(mol_fps):
        out_mols = simm.find_sim(mol, screen_fps, threshold)
        # The mol(1) is the smiles of the mol
        out_d["SCREEN "+str(i)] = {"IN_MOL": Chem.MolToSmiles(mol["RDMOL"], isomericSmiles=True), "OUT_MOLS": remove_keys(out_mols)}
    # Now return - preferably as json
    return json.dumps(out_d)
