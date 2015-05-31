from __future__ import absolute_import
from celery import shared_task
from rdkit import Chem
import json
from conf_gen.models import *  
from conf_gen.functions import *

def store_confs(cmpd, confs):
    """Function to store all the conformers for a given compound"""  
    # Loop through the conformations
    out_confs = []
    for conf in confs:
        new_conf = Conformation()
        new_conf.cmpd_id = cmpd
        new_conf.sdf_text = confs[conf]["source"]
        new_conf.rmsd = confs[conf]["values"]["RMSD"]
        new_conf.energy = confs[conf]["values"]["ENERGY"]
        new_conf.diff_energy = confs[conf]["values"]["ENERGY_DIFF"]
        new_conf.save()


@shared_task
def start_conf_gen(run_id):
    """Function to start a conformer generation"""
    cr = ConformerRun.objects.get(pk=run_id)
    # Get the ligs available for this docking
    my_ligs = get_lib_ligs(cr.lib_id)
    tot = len(my_ligs)
    old = -1
    cr.output_text  = "["
    for i, lig in enumerate(my_ligs):
        # Flag to stop, so jobs can be cancelled
        if cr.stop_flag == True:
            break
        if float(i) / float(tot) < cr.completion:
            continue
        # Get the RDMol
        rdmol = Chem.MolFromSmiles(str(lig.smiles))
        # Now make the conformations
        confs = generate_confs(rdmol, cr.num_mols, cr.min_rmsd, cr.max_energy, cr.max_energy_diff)
        # NOW STORE THE CONFORMATIONS
        if confs:
            store_confs(lig, confs)
        cr.output_text += str(confs)
        # Assign the completion to a value
        cr.completion = float(i) / float(tot)
        print cr.completion*100, "% COMPLETE"
        # Now save this as the result
        cr.save()
    

def store_moments(conf, moments, desc_type):
    """Function to store the moments generated in the above analysis"""
    this_desc = Descriptor()
    this_desc.conf_id = conf
    this_desc.desc_type = desc_type
    try:
        this_desc.validate_unique()
        this_desc.save()
    except ValidationError:
        return 
    # Now cycle through the moments and store them in the database
    for i, moment in enumerate(moments[0]):
        dv = DescVal()
        dv.bit_ind = i
        dv.value = moment
        dv.descriptor_id = this_desc
        # Now need to validate unique because there should be no duplicates
        dv.save()


@shared_task
def start_moment_gen(run_id):
    cr = MomentRun.objects.get(pk=run_id)
    # Get the ligs available for this docking
    my_ligs = get_lib_confs(cr.lib_id)
    tot = len(my_ligs)
    old = -1
    print my_ligs
    cr.output_text  = ""
    for i, lig in enumerate(my_ligs):
        # Flag to stop, so jobs can be cancelled
        if cr.stop_flag == True:
            break
        if float(i) / float(tot) < cr.completion:
            continue
        rdmol = Chem.MolFromMolBlock(str(lig.sdf_text))
        # Get the RDMol from the 3D information available
        moments = generate_moments(rdmol)
        # NOW ASSIGN THE DESCRIPTORS TO THIS CONFORMATION
        store_moments(lig, moments, "USRCAT")
        cr.output_text += str({"format": "mol", "source": lig.sdf_text, "values": {"moments": moments[0]} })
        # Assign the completion to a value
        cr.completion = float(i) / float(tot)
        print cr.completion*100, "% COMPLETE"
        cr.save()
