from __future__ import absolute_import

from celery import shared_task
from docking_runs.models import DockingRun
from rdkit import Chem
import json

@shared_task
def start_docking(run_id):
    """Function to start a docking"""
    dr = DockingRun.objects.get(pk=run_id)
    if dr.env_id.docking_function == "VINA":
        try:
            from VINADOCKING.functions import get_lib_ligs, run_docking
        except ImportError:
            print "VINA docking suite not available"
    else:
        print "Undefined function"
        return None    
    # Get the ligs available for this docking
    my_ligs = get_lib_ligs(dr.lib_id.uuid)
    tot = len(my_ligs)
    old = -1
    print "PROCESSED",tot ," LIGANDS"
    print "CARRYING OUT DOCKING" 
    dr.output_text  = ""
    out_l = []
    for i, lig in enumerate(my_ligs):
        # Flag to stop, so jobs can be cancelled
        if dr.stop_flag == True:
            break
        if float(i) / float(tot) < dr.completion:
            continue
        # Carry out the docking
        my_results = run_docking("/DockingSetup/"+dr.env_id.uuid, lig, "/DockingLibs/"+dr.lib_id.uuid, "/DockingRuns/"+dr.uuid)
        # Append the result text
        for mol in my_results:
             out_l.append({"source": Chem.MolToMolBlock(mol), "format": "mol","values": {"score": mol.GetProp("SCORE")}})
        dr.output_text = json.dumps(out_l)
        # Assign the completion to a value
        dr.completion = float(i) / float(tot)
        print dr.completion*100, "% COMPLETE"
        # Now save this as the result
        dr.save()
