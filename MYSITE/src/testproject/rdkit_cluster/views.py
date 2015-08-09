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
from django.views.decorators.gzip import gzip_page
from json_function.json_parse import remove_keys
import ast
import urllib2
import urllib
from mol_parsing.functions import request_handler, process_input


def index(request):
    return HttpResponse("WELCOME TO INDEX")


@gzip_page
@csrf_exempt
def cluster_simple(request):
        # Read the mols
    # Take the smiles in the request object
    if "dump_out" in request.GET:
        return HttpResponse(json.dumps(str(request))+"\nBODY:" + request.body)
    mol_type, screen_lib, fp_method, sim_method, threshold, params = request_handler(request)
    # Now return the process
    return process_input(fp_method, sim_method, screen_lib, mol_type, threshold, params)
