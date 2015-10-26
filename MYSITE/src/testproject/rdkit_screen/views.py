from django.shortcuts import render
import urllib
from django.http import HttpResponse
from rdkit_screen.functions import SimMethods, LibMethods, FPMethods
from rdkit import Chem
import json
from django.views.decorators.csrf import csrf_exempt
from django.views.decorators.gzip import gzip_page
from json_function.json_parse import remove_keys
import ast
import urllib2
from mol_parsing.functions import request_handler, process_input
import CloseableQueue
import json

def index(request):
    out_d = [
    {
    "id":"rdkit.screening.simple",
    "name":"RDKit screening",
    "description":"RDKit simple descriptor based screening",
    "tags":["virtualscreening", "screening", "moleculardescriptors", "fingerprints", "rdkit"],
    "paths":["/Chemistry/Toolkits/RDKit/Screening","/Chemistry/Screening"],
    "owner":"Tim Dudgeon <tdudgeon@informaticsmatters.com>",
    "layers":["public"],
    "inputClass":"com.im.lac.types.MoleculeObject",
    "outputClass":"com.im.lac.types.MoleculeObject",
    "inputType":"ARRAY",
    "outputType":"ARRAY",
    "accessModes":[
    {
        "id":"asyncHttp",
        "name":"Immediate execution",
        "description":"Execute as an asynchronous REST web service",
        "executionEndpoint":"screen_simple",
        "endpointRelative":True,
        "jobType":"com.im.lac.job.jobdef.AsyncHttpProcessDatasetJobDefinition",
        "parameters":[
            {
            "type":"STRUCTURE",
            "key":"smiles",
            "label":"Query",
            "description":"Query structure (as smiles)"
            },
            {
            "type":"FLOAT",
            "key":"threshold",
            "label":"Similarity Cuttoff",
            "description":"Similarity score cuttoff between 0 and 1 (1 means identical)"
            },
            {
            "type":"STRING",
            "key":"fp_method",
            "label":"Fingerprint",
            "values": ["morgan","maccs","rdkit_topo","atom_pairs"],
            "description":"Fingerprint method"
            },
            {
            "type":"STRING",
            "key":"metric",
            "label":"Metric",
            "values": ["tanimoto","cosine","dice","tversky"],
            "description":"Similarity comparison metric"
            }
        ],
        "adapterClassName":"com.im.lac.services.job.service.adapters.HttpGenericParamsJobAdapter"
    }
    ]
    }]
    return HttpResponse(json.dumps(out_d))


@gzip_page
@csrf_exempt
def screen_simple(request):
    """View to take a smiles and then screen against a known library of actives"""
    import urllib
    # Take the smiles in the request object
    # Now get the library
    if "dump_out" in request.GET:
        return HttpResponse(json.dumps(str(request))+"\nBODY:" + request.body)
    mol_type, screen_lib, fp_method, sim_method, threshold, params = request_handler(request)
    if "smiles" in request.GET:
        smiles = request.GET["smiles"]
        scr_mols = CloseableQueue.CloseableQueue()
        [scr_mols.put({"RDMOL": Chem.MolFromSmiles(str(x))}) for x in str(smiles).split(".")]
        scr_mols.close()
    else:
        return HttpResponse("You must state a SMILES")
    # Now handle this file upload 
    return process_input(fp_method, sim_method, screen_lib, mol_type, threshold, params=None, scr_mols=scr_mols)
