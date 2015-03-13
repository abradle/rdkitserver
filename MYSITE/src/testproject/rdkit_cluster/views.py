# Create your views here.
from django.shortcuts import render
from django.http import HttpResponse
from rdkit_screen.functions import SimMethods, LibMethods, FPMethods
from rdkit_cluster.functions import ClusterMethods
from rdkit import Chem
import json, gzip
from rdkit.ML.Cluster import Butina
from StringIO import StringIO
from django.views.decorators.csrf import csrf_exempt
from json_function.json_parse import remove_keys
import ast

def index(request):
    return HttpResponse("WELCOME TO INDEX")


@csrf_exempt
def cluster(request):
    """View to take a file of compounds and cluster the compounds
    Returns an SDFile with the the cluster number as a property"""
    # Read the mols
    # Take the smiles in the request object
    my_j = request.body
    print my_j
    my_json = ast.literal_eval(my_j)
    if "SCREEN_LIB" in my_json:
        screen_lib = my_json["SCREEN_LIB"]
    else:
        return HttpResponse("MUST SPECIFY COMPOUNDS TO CLUSTER")
     # Get the library
    print "GETTING MOLS"
    libm = LibMethods(screen_lib)
    mols = libm.get_mols()
    print "GOT MOLS"
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
        threshold = 0.7
    # Now run the process
    # Get the library
    # Get the fps
    print "GOT ARGS"
    fpm = FPMethods(fp_method)
    screen_fps = fpm.get_fps(mols)
    print "GOT FPS"
    # Get the distance matrix for my mol(s)
    simm = SimMethods(sim_method)
    dists = []
    nfps = len(screen_fps)
    print nfps, "FPS"
    for i in range(1,nfps):
        sims = [x["values"]["similarity"] for x in simm.find_sim(screen_fps[i], screen_fps[:i], -1.0)]
        # The mol(1) is the smiles of the mol
        dists.extend([1-x for x in sims])
    print "GOT SIMILARITIES"
    # now cluster the data:
    cs = Butina.ClusterData(dists,nfps,threshold,isDistData=True)
    # Out mols
    out_mols = []
    for i, c in enumerate(cs):
        for mol_ind in c:
            my_mol = mols[mol_ind]
            my_mol["values"]["cluster"] = i
            out_mols.append(my_mol)
    return HttpResponse(json.dumps(remove_keys(out_mols)))
