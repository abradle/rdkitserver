from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from mol_parsing import rdkit_parse

def morgan(m):
    return AllChem.GetMorganFingerprintAsBitVect(m,2)

class FPMethods():
    def __init__(self, fp_method):
        self.fp_method = fp_method
# Now define the functions
        self.f_dict = { "morgan": morgan}
        if self.fp_method not in self.f_dict:
            self.fp_method = None
    def get_fps(self, mols):
        error_counter = 0
        # Now set the FPs by iterating through the list -> get exceptions here
        for i, m in enumerate(mols):
            try:
                my_fp = self.f_dict[self.fp_method](m["RDMOL"])
            except:
                error_counter +=1
                my_fp = None
            mols[i]["FP"] = my_fp
        if error_counter:
            print "ERROR CREATING", str(i), "FINGERPRINTS"
        if error_counter == i+1:
            # If they have all failed then fail
            return None
        else:
            return mols

def parse_json_mols(mols):
    """Function to parse json mol objs"""
    out_mols = []
    for i, m in enumerate(mols):
        rdmol = rdkit_parse.parse_mol_json(m)
        if rdmol is None:
            print "NONE MOL"
            continue
        m["RDMOL"] = rdmol
        out_mols.append(m)
    return out_mols

def add_values_dict(my_mols):
    """Function to add a values dict to the JSON if it doesn't exist"""
    out_mols = []
    for mol in my_mols:
        if "values" in mol:
            pass
        else:
            mol["values"] = {}
        out_mols.append(mol)
    return out_mols


def sdf_mols_to_json(mols):
    """Function to parse SDF mols and return a JSON back"""
    out_mols = []
    
    for i,m in enumerate(mols.split("$$$$")):
        print m
        out_mols.append({"source": m, "format": "mol", "values": {}})
    return out_mols


def parse_sdf_mols(mols):
    """Function to get a json of molecules from an SDF"""
    my_json = sdf_mols_to_json(mols)
    return parse_json_mols(my_json)    


class LibMethods():
    def __init__(self, lib_mols, lib_type="JSON"):
        self.lib_mols = lib_mols
        self.f_d = {"JSON": parse_json_mols, "SDF": parse_sdf_mols}
        if lib_type in self.f_d:
            self.lib_type = lib_type
        else:
            self.lib_type = None

    def get_mols(self):
        # Get the molecules
        my_mols = self.f_d[self.lib_type](self.lib_mols)
        my_mols = add_values_dict(my_mols)
        return my_mols


def tanimoto(mol, lib_in):
    return DataStructs.BulkTanimotoSimilarity(mol,lib_in)

class SimMethods():
    def __init__(self, sim_meth):
        self.sim_meth = sim_meth
        self.f_dict = {"tanimoto": tanimoto}
        if self.sim_meth not in self.f_dict:
            self.sim_meth = None
    def find_sim(self, mol, lib_in, threshold):
        my_sim = self.f_dict[self.sim_meth](mol["FP"],[x["FP"] for x in lib_in])
        # Now make the out list to return
        out_ans = []
        for i,sim in enumerate(my_sim):
            if sim > threshold:
                my_ele = lib_in[i]
                my_ele["values"]["similarity"] = sim
                out_ans.append(my_ele)
        return out_ans
