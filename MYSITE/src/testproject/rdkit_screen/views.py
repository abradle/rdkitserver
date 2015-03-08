from django.shortcuts import render
from django.http import HttpResponse
from rdkit_screen.functions import SimMethods, LibMethods, FPMethods
from rdkit import Chem
import json

def index(request):
    return HttpResponse("WELCOME TO INDEX")

def screen(request):
    """View to take a smiles and then screen against a known library of actives"""
    # Take the smiles in the request object
    if "SMILES" in request.GET:
        smiles = request.GET["SMILES"]
        scr_mols = [Chem.MolFromSmiles(str(smiles))]
    else:
        return HttpResponse("You must state a SMILES")

    if "FP_METHOD" in request.GET:
        fp_method = request.GET["FP_METHOD"]
    else:
        fp_method = "morgan"

    if "SIM_METHOD" in request.GET:
        sim_method = request.GET["SIM_METHOD"]
    else:
        sim_method = "tanimoto"

    if "THRESHOLD" in request.GET:
        threshold = float(request.GET["THRESHOLD"])
    else:
        threshold = 0.7

    if "SCREEN_LIB" in request.GET:
        screen_lib = request.GET["SCREEN_LIB"]
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
    for mol in mol_fps:
        out_mols = simm.find_sim(mol, screen_fps, threshold)
        # The mol(1) is the smiles of the mol
        out_d[mol[1]] = out_mols
    # Now return - preferably as json
    return HttpResponse(json.dumps(out_d))
