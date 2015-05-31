from rdkit import Chem
from usrcat.toolkits.rd import generate_moments
from rdkit.Chem import AllChem
from django.core.exceptions import ValidationError
from conf_gen.models import *
import uuid

def get_lib_ligs(lib_id):
    """Function to get the ligands for this libary"""
    # Get the compounds within it 
    cmpds = lib_id.cmpd_id.all()
    return cmpds


def get_lib_confs(lib_id):
    """Function to get the conformers for this library"""
    # Get this library object
    # Get the compounds within it 
    confs = lib_id.cmpd_id.all().values_list("conformation", flat=True)
    return Conformation.objects.filter(pk__in=confs)


#### WRIITEN IN PYTHON -> WOULD BE MUCH FASTER IN C++
def generate_conformations(mol, num_confs=10, opt="MMFF"):
    """Function to generate a conformation for a given molecule and return the RMSD and the energy"""
    out_energy = []
    for i in range(int(num_confs)):
        # clean up the conformation
        if opt is "MMFF":
            ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=i)
        else:
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=i)
        if not ff:
            print "Error generating forcefield"
            return None, None, None
        ff.Initialize()
        n = 4
        more = ff.Minimize()
        while more and n:
            more = ff.Minimize()
            n -= 1
        # Get the RMSD
        out_energy.append(ff.CalcEnergy())
    # Now do all against all RMSD
    out_rmsd = AllChem.GetConformerRMSMatrix(mol, atomIds=None, prealigned=False)
    return (mol, out_energy, out_rmsd)


def generate_confs(rdmol, num_mols, min_rmsd=0.5, max_energy=None, max_energy_diff=None):
    """Function to generate diverse, low energy, low energy spread conformations for a given molecule"""
    # Take these molecules and embed molecules
    AllChem.EmbedMultipleConfs(rdmol, numConfs=int(num_mols), clearConfs=True, useRandomCoords=True, forceTol=0.00000000000000001)
    # Now minimise them
    mol, energy, rmsd = generate_conformations(rdmol, num_confs=num_mols, opt="MMFF")
    if not mol:
        return None
    min_energy = min(energy)
    out_d = {}
    # Now circle through this
    for i in range(int(num_mols)):
        if i != 0:
            my_rmsd = min(rmsd[i-1:i+i-1])
            print my_rmsd
            if min_rmsd:
                if my_rmsd < min_rmsd:
                    continue
        else:
            my_rmsd = -1.0
        if max_energy:
            if energy[i] > max_energy:
                continue
        out_d[i] = {"format": "mol", "values": {"ENERGY_DIFF": energy[i]-min_energy, "RMSD": my_rmsd, "ENERGY": energy[i]}, "source": Chem.MolToMolBlock(rdmol, confId=i)}
    return out_d


# Function to set up a conformer generation run
def setup_conf_run(input_lib, num_mols=None):
    """Function to set up a conformer generation run 
    Takes a library as input"""
    conf_run = ConformerRun()
    # The input lib to use
    conf_run.lib_id = input_lib
    # The uuid for this one
    conf_run.uuid = uuid.uuid4().get_hex() 
    # Completion
    conf_run.completion = 0.0
    if num_mols:
        conf_run.num_mols = num_mols 
    conf_run.save()
    return conf_run.pk


def register_lib(in_mols, lib_name=None):
    """Function to register a lib"""
    print "Creating lib"
    new_lib = InputLib()
    new_lib.uuid = uuid.uuid4().get_hex()
    new_lib.save()
    # The SDF blcok SDF info
    for my_mol in in_mols:
        smiles = Chem.MolToSmiles(my_mol, isomericSmiles=True) 
        inchi = Chem.MolToInchi(my_mol) 
        my_cmpd = Compound.objects.get_or_create(smiles=smiles, inchi=inchi)[0]
        new_lib.cmpd_id.add(my_cmpd)
    new_lib.num_mols = len(in_mols)
    # Lib name
    if lib_name:
        new_lib.lib_name = lib_name
    # Save the SD file in the directory
    new_lib.save()
    print "Created lib"
    print new_lib.pk
    return new_lib


# Function to set up a moment generation run
def setup_moment_run(input_lib):
    """Function to set up a moment generation run
    Takes a library as input"""
    mom_run = MomentRun()
    # The input lib to use
    mom_run.lib_id = input_lib
    # The uuid for this one
    mom_run.uuid = uuid.uuid4().get_hex() 
    # Completion
    mom_run.completion = 0.0
    mom_run.save()
    return mom_run.pk
