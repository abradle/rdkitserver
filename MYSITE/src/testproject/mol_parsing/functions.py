import ast
import urllib,zlib
from rdkit_screen.functions import SimMethods, LibMethods, FPMethods
from rdkit_cluster.functions import ClusterMethods
from rdkit import Chem
import json, gzip
from rdkit.ML.Cluster import Butina
from rdkit.ML.Cluster import Butina
from django.http import HttpResponse
from json_function.json_parse import remove_keys
from StringIO import StringIO

def do_screen(simm, screen_fps, threshold, mol_fps):
    """Function to peform a scren on a library of mols"""
    # Store the result in this
    out_d = {}
    for i, mol in enumerate(mol_fps):
        out_mols = simm.find_sim(mol, screen_fps, threshold)
        # The mol(1) is the smiles of the mol
        out_d = remove_keys(out_mols)
        # Now return - preferably as json
        return HttpResponse(json.dumps(out_d))


def do_clustering(simm, screen_fps, threshold, mols):
    """Function to peform the clustering fro m alibrary"""
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


def process_input(fp_method, sim_method, screen_lib, screen_type, threshold, params=None, scr_mols=None):
    # Now run the process
    # Get the library
    # Get the fps
    if screen_lib is None:
        return HttpResponse(screen_type)
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
    if scr_mols:
        mol_fps = fpm.get_fps(scr_mols)
        if not mol_fps:
            return HttpResponse("ERROR PRODUCING FINGERPRINTS (FOR MOLECULES)")
        else:
            return do_screen(simm, screen_fps, threshold, mol_fps)
    else:
        return do_clustering(simm, screen_fps, threshold, mols)


def find_lib_type(request):
    """Function to handle different types of incoming molecule libraries"""
    if "HTTP_CONTENT_ENCODING" in request.META:
        if request.META["HTTP_CONTENT_ENCODING"] == "gzip":
            in_data = urllib.unquote(request.body)
            compressedFile = StringIO(in_data)
            decompressedFile = gzip.GzipFile(fileobj=compressedFile)
            data = decompressedFile.read()
        else:
            print "NOT VALID INPUT ENCODING"
            return None, "DISALLOWED ENCODING TYPE"
    else:
        data = urllib.unquote(request.body).decode('utf8')
    if request.META["CONTENT_TYPE"] == "application/json":
        try:
            screen_lib = ast.literal_eval(data)
            mol_type = "JSON"
        except:
            print "NOT VALID JSON"
            return None, "NOT VALID JSON"
    elif request.META["CONTENT_TYPE"] == "chemical/x-mdl-sdfile":
        screen_lib = data
        mol_type = "SDF"
    else:
        try:
            screen_lib = dict(request.POST).keys()[0]
            screen_lib = ast.literal_eval(str(screen_lib))
            mol_type = "JSON"
        except IndexError:
            return None, "NONE OR DISALLOWED CONTENT TYPE"
    return screen_lib, mol_type



def request_handler(request):
    """Function to handle a generic request"""
    # Now handle this file upload 
    screen_lib, mol_type = find_lib_type(request)
    # Make the fingerprints
    if "fp_method" in request.GET:
        fp_method = request.GET["fp_method"]
    else:
        fp_method = "morgan"
    if "params" in request.GET:
        params = request.GET["params"].split(",")
    else:
        params = None
    if "sim_method" in request.GET:
        sim_method = request.GET["sim_method"]
    else:
        sim_method = "tanimoto"
    if "threshold" in request.GET:
        threshold = float(request.GET["threshold"])
    else:
        threshold = 0.7
    return mol_type, screen_lib, fp_method, sim_method, threshold, params 
