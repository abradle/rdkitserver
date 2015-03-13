from django.shortcuts import render
from django.http import HttpResponse
from rdkit_screen.functions import SimMethods, LibMethods, FPMethods
from rdkit import Chem
import json
from django.views.decorators.csrf import csrf_exempt
from json_function.json_parse import remove_keys
import ast

def index(request):
    return HttpResponse("WELCOME TO INDEX")


@csrf_exempt
def screen(request):
    """View to take a smiles and then screen against a known library of actives"""
    # Take the smiles in the request object
    my_j = request.body
    my_json = ast.literal_eval(my_j)
    if "SMILES" in my_json:
        smiles = my_json["SMILES"]
        scr_mols = [{"RDMOL": Chem.MolFromSmiles(str(x)} for x in str(smiles).split(".")]
    else:
        return HttpResponse("You must state a SMILES")

    if "FP_METHOD" in my_json:
        fp_method = my_json["FP_METHOD"]
    else:
        fp_method = "morgan"

    if "SIM_METHOD" in my_json:
        sim_method = my_json["SIM_METHOD"]
    else:
        sim_method = "tanimoto"

    if "THRESHOLD" in my_json:
        threshold = float(my_json["THRESHOLD"])
    else:
        threshold = 0.7

    if "SCREEN_LIB" in my_json:
        screen_lib = my_json["SCREEN_LIB"]
    else:
        screen_lib = "default"
     # Now run the process
     # Get the library
    libm = LibMethods(screen_lib)
    mols = libm.get_mols()
    # Get the fps
    fpm = FPMethods(fp_method)
    screen_fps = fpm.get_fps(mols)
    # Get the fp for my mol(s)
    mol_fps = fpm.get_fps(scr_mols)
    simm = SimMethods(sim_method)
    # Store the result in this
    out_d = {}
    for i, mol in enumerate(mol_fps):
        out_mols = simm.find_sim(mol, screen_fps, threshold)
        # The mol(1) is the smiles of the mol
        out_d["SCREEN "+str(i)] = {"IN_MOL": Chem.MolToSmiles(mol["RDMOL"], isomericSmiles=True), "OUT_MOLS": remove_keys(out_mols)}
    # Now return - preferably as json
    return HttpResponse(json.dumps(out_d))
