from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs


def morgan(m):
    return AllChem.GetMorganFingerprintAsBitVect(m,2)

class FPMethods():
    def __init__(self, fp_method):
        self.fp_method = fp_method
# Now define the functions
        self.f_dict = { "morgan": morgan}
    def get_fps(self, mols):
        my_fps = []
        for m in mols:
            my_fps.append([self.f_dict[self.fp_method](m), Chem.MolToSmiles(m, isomericSmiles=True)])
        return my_fps

class LibMethods():
    def __init__(self, lib_name):
        self.lib_name = lib_name
    def get_mols(self):
        if self.lib_name == "default":
            my_mols = Chem.SmilesMolSupplier("/MYSITE/src/testproject/data/nci1000.smiles", delimiter="\t")#[Chem.MolFromSmiles("CCC"), Chem.MolFromSmiles("CC")] 
        else:
            print "UNDEFINED LIBRARY"
            my_mols = None
        return [x for x in my_mols if x != None] 

def tanimoto(mol, lib_in):
    return DataStructs.BulkTanimotoSimilarity(mol,lib_in)

class SimMethods():
    def __init__(self, sim_meth):
        self.sim_meth = sim_meth
        self.f_dict = {"tanimoto": tanimoto}
    def find_sim(self, mol, lib_in, threshold):
        my_sim = self.f_dict[self.sim_meth](mol[0],[x[0] for x in lib_in])
        # Now make the out list to return
        out_ans = []
        for i,sim in enumerate(my_sim):
            if sim > threshold:
                out_ans.append([lib_in[i][1], sim])
        return out_ans
