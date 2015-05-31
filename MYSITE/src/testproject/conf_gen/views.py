from django.shortcuts import render
import urllib
from django.http import HttpResponse
from rdkit import Chem
import json
from django.views.decorators.csrf import csrf_exempt
from django.views.decorators.gzip import gzip_page
from json_function.json_parse import remove_keys
import ast
import urllib2
from mol_parsing.functions import request_handler, process_input
from conf_gen.functions import *
from rdkit_screen.functions import LibMethods
from conf_gen.tasks import *


@csrf_exempt
def gen_confs(request):
    """View to take a library of input molecules and generate conformations from them"""
    mol_type, screen_lib, fp_method, sim_method, threshold, params = request_handler(request)
    # Now handle this file upload 
    libm = LibMethods(screen_lib, mol_type)
    my_mols = libm.get_mols()
    my_lib = register_lib([x["RDMOL"] for x in my_mols])
    my_pk = setup_conf_run(my_lib, num_mols=None)    
    response = start_conf_gen.delay(my_pk)
    return HttpResponse(str(my_pk))


@csrf_exempt
def gen_moments(request):
    """View to take a library of molecules and return them annotated with USRCAT moments"""
    if "CONF" in request.GET:
        conf_id = int(request.GET["CONF"])
    else:
        return HttpResponse("ERROR CONF RUN NOT SPECIFIED")
    # Get the library
    lib = ConformerRun.objects.get(pk=conf_id).lib_id
    # Now make set of this celery task
    my_pk = setup_moment_run(lib)
    response = start_moment_gen.delay(my_pk)
    return HttpResponse(str(my_pk))   


def get_confs(request):
    """View to get the progress for a given conformer generation job"""
    if "CONF" in request.GET:
        conf_id = int(request.GET["CONF"])
    else:
        return HttpResponse("ERROR CONF RUN NOT SPECIFIED")
    # Return a JSON of information along`1
    conf = ConformerRun.objects.filter(pk=conf_id).values("output_text", "pk", "completion")
    return HttpResponse(json.dumps(ast.literal_eval(str((conf)))), content_type='application/json')

def get_moments(request):
    """View to get the progress for a given conformer generation job"""
    if "CONF" in request.GET:
	conf_id = int(request.GET["CONF"])
    else:
	return HttpResponse("ERROR CONF RUN NOT SPECIFIED")
    # Return a JSON of information along`1
    conf = MomentRun.objects.filter(pk=conf_id).values("output_text", "pk", "completion")
    return HttpResponse(json.dumps(ast.literal_eval(str((conf)))), content_type='application/json')


def get_progress(request):
    """View to get the progress for a given conformer generation job"""
    if "CONF" in request.GET:
        conf_id = int(request.GET["CONF"])
    else:
        return HttpResponse("ERROR CONF RUN NOT SPECIFIED")
    # Return a JSON of information along`1
    conf = ConformerRun.objects.filter(pk=conf_id).values("pk", "completion")
    return HttpResponse(json.dumps(ast.literal_eval(str(conf))))
