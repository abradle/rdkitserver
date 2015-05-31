from django.db import models


class Compound(models.Model):
    """Unique molecular entities - so that we don't generat conformers multiple times unnecesarily"""
    # Go uique on smiles - because this just to save time, so I only ever want totally identical molecules to show up
    smiles = models.TextField(db_index=True, unique=True)
    inchi = models.TextField(db_index=True)


class Conformation(models.Model):
    """A conformation for a compound - stores the 3D coords and how it was made"""
    cmpd_id = models.ForeignKey(Compound)
    sdf_text = models.TextField()
    rmsd = models.FloatField()
    energy = models.FloatField()
    diff_energy = models.FloatField()


class Descriptor(models.Model):
    """A descriptor for a given molecule - e.g. a USRCAT fingerprint"""
    conf_id = models.ForeignKey(Conformation)
    desc_type = models.TextField()

    class Meta:
        unique_together = ("conf_id", "desc_type")


class DescVal(models.Model):
    """A class to store the individual descriptor bits"""
    value = models.FloatField()
    bit_ind = models.FloatField()
    descriptor_id = models.ForeignKey(Descriptor)

    class Meta:
        unique_together = ("bit_ind", "descriptor_id")


class InputLib(models.Model):
    """Place to store an input library of molecules for processing"""
    # Uuid 
    uuid = models.TextField(unique=True)
    # The compounds to which we are referring
    cmpd_id = models.ManyToManyField(Compound)
    # Error mols 
    errors = models.FloatField(default=0.0)
    # Error text
    error_log = models.TextField(null=True)
    # Num mols
    num_mols = models.IntegerField(default=0)
    # Lib name
    lib_name = models.TextField(null=True)


class MomentRun(models.Model):
    """A model to hold information about a process generating moments"""
    # The input lib to use
    lib_id = models.ForeignKey(InputLib)
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


class ConformerRun(models.Model):
    """A model to hold information about a process generating coformers"""
    # The input lib to use
    lib_id = models.ForeignKey(InputLib)
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
    # Parameters for the run
    num_mols = models.FloatField(default=50.0) 
    min_rmsd = models.FloatField(null=True)
    max_energy = models.FloatField(null=True)
    max_energy_diff = models.FloatField(null=True)
