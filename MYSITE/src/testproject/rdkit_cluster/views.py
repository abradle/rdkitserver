# Create your views here.
from django.shortcuts import render
from django.http import HttpResponse
from rdkit_screen.functions import SimMethods, LibMethods, FPMethods
from rdkit_cluster.functions import ClusterMethods
from rdkit import Chem
import json, gzip
from rdkit.ML.Cluster import Butina
from StringIO import StringIO


def index(request):
    return HttpResponse("WELCOME TO INDEX")


def cluster(request):
    """View to take a file of compounds and cluster the compounds
    Returns an SDFile with the the cluster number as a property"""
    # Read the mols
    #my_text = request.FILES['myfile'].read()
    mols = [x for x in Chem.ForwardSDMolSupplier(request.FILES['myfile']) if x is not None]
    # Make the fingerprints
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
    # Now run the process
    # Get the library
    # Get the fps
    fpm = FPMethods(fp_method)
    screen_fps = fpm.get_fps(mols)
    # Get the distance matrix for my mol(s)
    simm = SimMethods(sim_method)
    dists = []
    nfps = len(screen_fps)
    for i in range(1,nfps):
        sims = [x[1] for x in simm.find_sim(screen_fps[i], screen_fps[:i], -1.0)]
        # The mol(1) is the smiles of the mol
        dists.extend([1-x for x in sims])
    # now cluster the data:
    cs = Butina.ClusterData(dists,nfps,threshold,isDistData=True)
    # Out mols
    out_mols = []

    sio = StringIO()
    w = Chem.SDWriter(sio)
    for i, c in enumerate(cs):
        for mol_ind in c:
            my_mol = mols[mol_ind]
            my_mol.SetProp("cluster_num", str(i))
            w.write(my_mol)
    w.flush()
    return HttpResponse(json.dumps(sio.getvalue()))
