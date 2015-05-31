from docking_runs.models import  *
from docking_runs.functions import *
from docking_runs.tasks import start_docking
from mol_parsing.functions import request_handler
from django.http import HttpResponse
import json
import ast
from django.views.decorators.csrf import csrf_exempt
from rdkit_screen.functions import LibMethods
# View to start a docking
@csrf_exempt
def start_dock(request):
    # Get the library -> either as an ID or as an upload
    if "LIB" in request.GET:
        lib_id = int(request.GET["LIB"])
    else:
        try:
            mol_type, screen_lib, fp_method, sim_method, threshold, params = request_handler(request)
        except:
            screen_lib = None
        if not screen_lib:
            return HttpResponse("NO LIB SPECIFIED")
        else:
            libm = LibMethods(screen_lib, mol_type)
            my_mols = libm.get_mols()
            my_lib = register_lib([x["RDMOL"] for x in my_mols])
            lib_id = my_lib.pk
        # Now process this libary -> save to disk
    # Get the environment
    if "ENV" in request.GET:
        env_id = int(request.GET["ENV"])
        my_env = DockingEnv.objects.filter(pk=env_id)
        if not my_env:
            return HttpResponse("INVALID ENV")
    else:
        return HttpResponse("ERROR DOCKING ENV NOT SPECIFIED")
    my_lib = DockingLib.objects.filter(pk=lib_id)
    if not my_lib:
        return HttpResponse("INVALID LIBARAY")
    run_id = register_docking(my_env[0], my_lib[0])
    # Now set this opff running
    response = start_docking.delay(run_id)
    return HttpResponse(run_id)

# View to check on a docking process
def check_dock(request):
    if "DOCK" in request.GET:
        dock_id = int(request.GET["DOCK"])
    else:
        return HttpResponse("ERROR DOCKING RUN NOT SPECIFIED")
    # Return a JSON of information along`1
    dock = DockingRun.objects.filter(pk=dock_id).values("pk", "completion")
    return HttpResponse(json.dumps(ast.literal_eval(str(dock))))


# View to get a list of the libraries availabel
def libs_avail(request):
    dls =DockingLib.objects.filter().values("pk", "num_mols", "lib_name")
    out_json = json.dumps(ast.literal_eval(str(dls)))
    return HttpResponse(out_json, content_type='application/json')


# View to get a list of the docking runs avaialbele
def docks_avail(request):
    dls = DockingEnv.objects.filter().values("pk", "pdb_code", "uniprot", "smiles", "docking_function")
    return HttpResponse(json.dumps(ast.literal_eval(str((dls)))), content_type='application/json')


def get_output(request):
    if "DOCK" in request.GET:
        dock_id = int(request.GET["DOCK"])
    else:
        return HttpResponse("ERROR DOCKING RUN NOT SPECIFIED")
    dock = DockingRun.objects.filter(pk=dock_id).values("pk", "completion", "output_text")  
    return HttpResponse(json.dumps(ast.literal_eval(str((dock)))), content_type='application/json')
