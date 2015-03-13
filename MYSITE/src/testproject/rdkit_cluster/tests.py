from django.test import TestCase
from rdkit.ML.Cluster import Butina

from rdkit_screen.functions import SimMethods, LibMethods, FPMethods
from rdkit import Chem
import json
import ast
from json_function.json_parse import remove_keys


def cluster():
    """View to take a file of compounds and cluster the compounds
    Returns an SDFile with the the cluster number as a property"""
    # Read the mols
    # Take the smiles in the request object
    my_json = {}
    libm = LibMethods(ast.literal_eval(open("smis.json").read()))
     # Get the library
    mols = libm.get_mols()
    # Make the fingerprints
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
        threshold = 0.9
    # Now run the process
    # Get the library
    # Get the fps
    fpm = FPMethods(fp_method)
    screen_fps = fpm.get_fps(mols)
    # Get the distance matrix for my mol(s)
    simm = SimMethods(sim_method)
    dists = []
    nfps = len(screen_fps)
    #www now cluster the data:
    for i in range(1,nfps):
        sims = [x["values"]["similarity"] for x in simm.find_sim(screen_fps[i], screen_fps[:i], -1.0)]
        # The mol(1) is the smiles of the mol
        dists.extend([1-x for x in sims])
    # now cluster the data:
    cs = Butina.ClusterData(dists,nfps,threshold,isDistData=True)
    # Out mols
    out_mols = []
    for i, c in enumerate(cs):
        for mol_ind in c:
            my_mol = mols[mol_ind]
            my_mol["values"]["cluster"] = i
            out_mols.append(my_mol)
    return json.dumps(remove_keys(out_mols))
                                                            
