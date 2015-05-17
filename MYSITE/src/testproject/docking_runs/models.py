from django.db import models

# Docking environments-> e.g. a pre set docking run
class DockingEnv(models.Model):
    """A model to hold information about a given docking environment"""
    # The program used -> not essential
    docking_function = models.TextField(null=True)
    # The uuid, this finds the appropriate directtoy
    uuid = models.TextField(unique=True)
    # The name of this entry for public consumption
    entry_name = models.TextField(null=True)
    # PDB code
    pdb_code = models.TextField()
    # Protein name, e.g. CDK2
    prot_name = models.TextField(null=True)
    # Protein uniprot
    uniprot = models.TextField(null=True)
    # Ligand site smiles
    smiles = models.TextField(null=True)
    # Meta_Data
    meta_data = models.TextField(null=True)


# Docking libs
class DockingLib(models.Model):
    """A model to contain information about a docking lib"""
    # Uuid 
    uuid = models.TextField(unique=True)
    # The SDF blcok SDF info
    sdf_info = models.TextField(null=True)
    # Error mols 
    errors = models.FloatField(default=0.0)
    # Error text
    error_log = models.TextField(null=True)
    # Num mols
    num_mols = models.IntegerField(default=0)
    # Lib name
    lib_name = models.TextField(null=True)


# Docking runs -> e.g. a docking run set up
class DockingRun(models.Model):
    """A model to hold information about a given docking run"""
    # The docking env it relates to
    env_id = models.ForeignKey(DockingEnv)
    # The docking lib to use
    lib_id = models.ForeignKey(DockingLib)
    # The uuid for this one
    uuid = models.TextField(unique=True)
    # Completion
    completion = models.FloatField(default=0.0)
    stop_flag = models.BooleanField(default=False)
    # Errors
    errros = models.FloatField(default=0.0)
    # Output result
    output_text = models.TextField(null=True)
    # Error log
    error_log = models.TextField(null=True)
