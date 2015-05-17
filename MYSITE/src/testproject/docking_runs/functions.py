from docking_runs.models import DockingRun, DockingEnv, DockingLib 
import uuid
from rdkit import Chem
import os


def register_docking(env_id, lib_id):
    """Function to register a docking"""
    new_dock = DockingRun()
    new_dock.lib_id = lib_id
    new_dock.env_id = env_id
    new_dock.uuid = uuid.uuid4().get_hex()
    new_dock.save()
    return new_dock.pk


def register_env(dock_type, pdb_code, uuid=None, name=None, prot_name=None, uniprot=None, smiles=None, meta_data=None):
    """Function to register a docking env"""
    # Create this new docking env
    de = DockingEnv()
    # Process the function
    if dock_type == "VINA":
        try:
            from VINADOCKING.functions import prepare_env
        except ImportError:
            print "VINA docking suite not available"
    if not uuid:
        conf_path = prepare_env(pdb_code, smiles)
        print conf_path
        uuid = os.path.split(conf_path)[0].split("/")[-1]
        print uuid
  #.split("/")[-1]
        print "UUID created: ", uuid
    else:
        pass
    # Now set the details
    de.docking_function = dock_type
    # The uuid, this finds the appropriate directtoy
    de.uuid = uuid
    # The name of this entry for public consumption
    if name:
        de.entry_name = name
    # PDB code
    if pdb_code:
        de.pdb_code = pdb_code
    # Protein name, e.g. CDK2
    if prot_name:
        de.prot_name = prot_name
    # Protein uniprot
    if uniprot:
        de.uniprot = uniprot
    # Ligand site smiles
    if smiles:
        de.smiles = smiles
    # Meta_Data
    if meta_data:
        de.meta_data = meta_data
    de.save()

    

def register_lib(in_mols, lib_name=None):
    """Function to register a lib"""
    print "Creating lib"
    new_lib = DockingLib()
    new_lib.uuid = uuid.uuid4().get_hex()
    os.mkdir("/DockingLibs/"+new_lib.uuid)
    my_sd = "/DockingLibs/"+new_lib.uuid+"/mols.sdf"
    out_f = Chem.SDWriter(my_sd)
    [out_f.write(x) for x in in_mols if x]
    out_f.close()
    new_lib.sdf_info = open(my_sd).read()
    new_lib.num_mols = len([x for x in in_mols if x])
    if lib_name:
        new_lib.lib_name = lib_name
    # Save the SD file in the directory
    new_lib.save()
    print "Created lib"
    print new_lib.pk
    return new_lib
