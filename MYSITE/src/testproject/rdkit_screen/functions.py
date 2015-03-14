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
    def get_fps(self, mols):
        # Now set the FP
        for i,m in enumerate(mols):
            my_fp = self.f_dict[self.fp_method](m["RDMOL"])
            mols[i]["FP"] = my_fp
        return mols

def parse_json_mols(mols):
    """Function to parse json mol objs"""
    out_mols = []
    for i,m in enumerate(mols):
        rdmol = rdkit_parse.parse_mol_json(m)
        if m is None:
            print "NONE MOL"
            continue
        m["RDMOL"] = rdmol
        out_mols.append(m)
    return out_mols


class LibMethods():
    def __init__(self, lib_mols, lib_type="JSON"):
        self.lib_mols = lib_mols
        self.lib_type = lib_type
    def get_mols(self):
        my_mols = parse_json_mols(self.lib_mols)
        return my_mols


def tanimoto(mol, lib_in):
    return DataStructs.BulkTanimotoSimilarity(mol,lib_in)

class SimMethods():
    def __init__(self, sim_meth):
        self.sim_meth = sim_meth
        self.f_dict = {"tanimoto": tanimoto}
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
