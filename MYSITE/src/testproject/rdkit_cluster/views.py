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
import urllib2

def index(request):
    return HttpResponse("WELCOME TO INDEX")


def process_input(fp_method, sim_method, screen_lib, screen_type, threshold):
    # Now run the process
    # Get the library
    # Get the fps
    libm = LibMethods(screen_lib, screen_type)
    if libm.lib_type is None:
        return HttpResponse("NOT RECOGNISED UPLOAD TYPE" + screen_type)
    mols = libm.get_mols()
    if not mols:
        return HttpResponse("NO VALID MOLECULES!!!")
    fpm = FPMethods(fp_method)
    if not fpm.fp_method:
        return HttpResponse("NOT A REGISTERED FINGERPRINT METHOD -> " + fp_method)
    # Now get the fingerprints
    screen_fps = fpm.get_fps(mols)
    if not screen_fps:
        return HttpResponse("ERROR PRODUCING FINGERPRINTS (FOR SCREENING LIBRARY)")
    # Now get the similarity metric
    simm = SimMethods(sim_method)
    if not simm.sim_meth:
        return HttpResponse("NOT A VALID SIMILARIY METRIC")
    # Now prodcue the distance matric
    dists = []
    nfps = len(screen_fps)
    for i in range(1, nfps):
        sims = [x["values"]["similarity"] for x in simm.find_sim(screen_fps[i], screen_fps[:i], -1.0)]
        # The mol(1) is the smiles of the mol
        dists.extend([1 - x for x in sims])
    # now cluster the data:
    cs = Butina.ClusterData(dists, nfps, threshold, isDistData=True)
    # Out mols is the list for caputring the clusters
    out_mols = []
    # Now loop trhough the clusters outputing the results
    for i, c in enumerate(cs):
        for mol_ind in c:
            my_mol = mols[mol_ind]
            my_mol["values"]["cluster"] = i
            out_mols.append(my_mol)
    # Now retrun the response
    return HttpResponse(json.dumps(remove_keys(out_mols)))

@csrf_exempt
def cluster_mol_body(request):
    """View to take the molecules in the body and params as POST attributes"""
    # Get the mols as the body
    mol_url = request.POST["SCREEN_LIB"]
    # Get the library from the URL specified
    r = urllib2.urlopen(mol_url)
    if "MOL_TYPE" in request.POST:
        mol_type = request.POST["MOL_TYPE"]
    else:
        mol_type = "JSON"
    # Pull it into a JSON
    if mol_type == "JSON":
        screen_lib = ast.literal_eval(r.read())
    else:
        screen_lib = r.read()
    libm = LibMethods(screen_lib, mol_type)
    mols = libm.get_mols()
    # Make the fingerprints
    if "FP_METHOD" in request.POST:
        fp_method = request.POST["FP_METHOD"]
    else:
        fp_method = "morgan"
    # Get the similarity method
    if "SIM_METHOD" in request.POST:
        sim_method = request.POST["SIM_METHOD"]
    else:
        sim_method = "tanimoto"
    # Get the threshold
    if "THRESHOLD" in request.POST:
        threshold = float(request.POST["THRESHOLD"])
    else:
        threshold = 0.7
    # Now do the actual processing
    return process_input(fp_method, sim_method,screen_lib, mol_type, threshold)

@csrf_exempt
def cluster(request):
    """View to take a file of compounds and cluster the compounds
    Returns an SDFile with the the cluster number as a property"""
    # Read the mols
    # Take the smiles in the request object
    my_j = request.body
    my_json = ast.literal_eval(my_j)
    if "mol_type" in my_json:
        mol_tye = my_json["mol_type"]
    else:
        mol_type = "JSON"
    if "SCREEN_LIB" in my_json:
        screen_lib = my_json["SCREEN_LIB"]
    else:
        return HttpResponse("MUST SPECIFY COMPOUNDS TO CLUSTER")
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
    return process_input(fp_method, sim_method, screen_lib, mol_type, threshold)



@csrf_exempt
def cluster_simple(request):
        # Read the mols
    # Take the smiles in the request object
    try:
        screen_lib = dict(request.POST).keys()[0]
    except IndexError:
        return HttpResponse("YOU MUST UPLOAD A LIBRARY")
    # Now get the screening lib
    screen_lib = ast.literal_eval(str(screen_lib))
    # Get the library type
    if "mol_type" in request.GET:
        mol_type = request.GET["mol_type"]
    else:
        mol_type = "JSON"
    # Make the fingerprints
    if "fp_method" in request.GET:
        fp_method = request.GET["fp_method"]
    else:
        fp_method = "morgan"
    if "sim_method" in request.GET:
        sim_method = request.GET["sim_method"]
    else:
        sim_method = "tanimoto"
    if "threshold" in request.GET:
        threshold = float(request.GET["threshold"])
    else:
        threshold = 0.7
    return process_input(fp_method, sim_method, screen_lib, mol_type, threshold)


