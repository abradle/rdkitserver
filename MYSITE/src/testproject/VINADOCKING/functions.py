from rdkit import Chem
import subprocess, urllib
from rdkit.Chem import AllChem
import tempfile
import uuid
import os
import sys


def get_pdb(my_pdbs):
    """Function to get all the PDBs for a given list -> and return a list of their paths on the file system"""
    print "Downloading PDB data"
    # Now get the SMILES and take out the most likely one if any available
    import xml.etree.ElementTree as ET
    import urllib, urllib2, sys
    # The dictionary to write stuff to
    out_d = {}
    # First get all the PDBs for this uniprot
    tot = len(my_pdbs)
    old = -1 
    for i, my_item in enumerate(my_pdbs):
        # Print the progress
        myuuid = uuid.uuid4().get_hex()
        my_dir = "/DockingSetup/"+myuuid
        if i * 100 / tot != old:
            old = i * 100 / tot
            sys.stdout.write("\rLoading structural data from PDB %d%% complete..." % old)
            sys.stdout.flush()
        pdb = my_item["pdb_code"]
        if "smiles" in my_item:
            smiles = my_item["smiles"]
        else:
            smiles = None
        out_d[pdb] = {}
        # Now get the pdb file data
        my_pdb = urllib.urlopen("http://www.rcsb.org/pdb/files/" + pdb.lower() + ".pdb").read()
        os.mkdir(my_dir)
        out_p = my_dir+"/"+pdb+".pdb"
        out_f = open(out_p, "w")
        out_f.write(my_pdb)
        out_f.close()
        # Write it to a file
        out_d[pdb]["path"] = out_p
        # Now get the smiles data - and pick the best oe
        my_s = urllib.urlopen("http://www.rcsb.org/pdb/rest/ligandInfo?structureId=" + pdb).read()
        e = ET.fromstring(my_s)
        my_mol = None
        for item in e[0]:
            # Now pick the most druglike one??? Most no H A - not very good
            try:
                try_mol = Chem.MolFromSmiles(str(item[4].text.strip()))
            except IndexError:
                continue
            if not try_mol:
                continue
            if my_mol and try_mol.GetNumHeavyAtoms() < my_mol.GetNumHeavyAtoms():
                continue
            my_mol = try_mol
        if my_mol and not smiles:
            out_d[pdb]["smiles"] = Chem.MolToSmiles(my_mol, isomericSmiles=True)
        elif smiles:
            out_d[pdb]["smiles"] = smiles
        else:
            out_d[pdb]["smiles"] = None
    return out_d

# Extract ligand
def assign_temp(block, smiles):
    """Helper function to create an RDMol from PDB information"""
    tmp = Chem.MolFromPDBBlock(block)
    template = Chem.MolFromSmiles(smiles)
    if not template:
        print smiles, " not recognised"
        return None
    AllChem.AddHs(template)
    try:
        mol = AllChem.AssignBondOrdersFromTemplate(template, tmp)
    except ValueError:
        print "DOESN'T FIT", smiles
        mol = None
    except:
        print "UNSPECIFED ERRROR"
        mol = None
    return mol


def get_ligs(lines):
    """Helper function to extract LIGAND data from a PDB file"""
    old_id = "NOTME"
    old_chain = ""
    out_d = []
    for line in lines:
        if line[:6] != "HETATM":
            continue
        res = line[17:20]
        if res in ["HOH", "EDO", "ACT"]:
                continue
        if old_id == "NOTME":
            line_list = [line]
            old_id = res
        elif res != old_id or int(line[22:26]) != old_chain:
            out_d.append("".join(line_list))
            line_list = [line]
            old_id = res
        else:
            line_list.append(line)
        old_chain = int(line[22:26])
    if old_id == "NOTME":
        return []
    out_d.append("".join(line_list))
    return out_d  


def remove_hetatm(lines):
    return "".join([line for line in lines if line[:6] not in ["HETATM", "CONECT"] ])


def extract_ligand(file_lines, smiles):
    mols = [assign_temp(block, smiles) for block in get_ligs(file_lines) if assign_temp(block, smiles) is not None]
    return mols

# Prepare protein
# Prepare the ligands
def prepare_mol(pdb_file, my_type):
    """Function to prepare an RDKit molecule using ADT"""
    pythonsh = '/mgltools_x86_64Linux2_1.5.6/bin/pythonsh'
    prepare_lig = '/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py'
    prepare_prot = '/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py'
    revert_pdbqt = '/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/pdbqt_to_pdb.py'
    # Set the out file
    if my_type=="lig":
       out_file = pdb_file.replace(".pdb",".pdbqt") 
       subprocess.call([pythonsh,prepare_lig,"-l",pdb_file,"-o", out_file])
    elif my_type=="prot":
       out_file = pdb_file.replace(".pdb",".pdbqt") 
       subprocess.call([pythonsh,prepare_prot,"-v","-r",pdb_file,"-o",out_file,"-A",'bonds_hydrogens'])
    elif my_type=="result":
       out_file = pdb_file.replace(".pdbqt",".pdb")
       subprocess.call([pythonsh, revert_pdbqt, "-f", pdb_file, "-o", out_file])
    return out_file


# Run docking
def run_vina(conf_file, result, my_lig=None):
    """Function to call a vina docking"""
    vina_path = "/usr/bin/vina"
    if my_lig:
        subprocess.call([vina_path, "--out", result, "--config", conf_file, "--ligand", my_lig])
    else:
        subprocess.call([vina_path, "--out", result, "--config", conf_file])

def centre_of_mass_from_rdmol(rdmol):
    """Function to return the unweighted Centre of Mass from an SD block.
    Takes an RDKit mol
    Returns three floats, the x,y and z coordinate for the centre of mass"""
    # Gives the atoms
    atoms = rdmol.GetAtoms()
    conf = rdmol.GetConformer()
    x_coord = y_coord = z_coord = 0.0
    numatoms = 0.0
    # Assume all heavy atoms have the same mass
    for atom in atoms:
        if atom.GetAtomicNum() == 1 or atom.GetSmarts() == "[*]":
            continue
        numatoms += 1.0
        coords = conf.GetAtomPosition(atom.GetIdx())
        if coords.x == 0.0 and coords.y == 0.0 and coords.z == 0.0:
            print sdf_text
            sys.exit()
        x_coord += float(coords.x)
        y_coord += float(coords.y)
        z_coord += float(coords.z)
    # Now we have all the coords -> we want to loop through
    if numatoms == 0:
        sys.stderr.write("Zero atoms error...")
        return None, None, None
    return int(x_coord / numatoms), int(y_coord / numatoms), int(z_coord / numatoms)


def parse_result(in_pdbqt): 
    """Function to parse a pdbqt file and parse the result
    Returns a list of mols with property set"""
    out_sd = in_pdbqt.replace(".pdbqt",".sdf")
    print "WRITING TO: ", out_sd
    out_mols = []
    # First split on model
    my_mols = open(in_pdbqt).read().split("MODEL")   
    for i, mol in enumerate(my_mols):
        if i == 0:
            continue
        my_file = tempfile.NamedTemporaryFile(mode='w', suffix='.pdbqt', delete=False)
        my_file.write(mol)
        my_name = my_file.name
        my_file.close()
        pdb_file = prepare_mol(my_file.name, "result")
        rdmol = Chem.MolFromPDBFile(pdb_file)
        if not rdmol:
            continue
        score = [line for line in mol.split("\n") if "REMARK VINA RESULT" in line][0].split()[3]
        rdmol.SetProp("SCORE", score)
        out_mols.append(rdmol)
    return out_mols 

def prepare_lib(in_mols):
    """Function to prepare a library for docking"""
    lib_id = "/DockingLibs/"+uuid.uuid4().get_hex()
    os.mkdir(lib_id)


def save_lib(in_mols, lib_id):
    """Function to save the library"""
    out_lib = lib_id+"/mols.sdf"
    if os.path.isfile(out_lib):
        out_sd = None
    else:
        print "WRITING TO: ", out_lib
        out_sd = Chem.SDWriter(out_lib)
    # Loop through the mols
    for i, mol in enumerate(in_mols):
        if mol is None:
            print "NONE MOL"
            continue
        if out_sd:
            out_sd.write(mol)
        out_file = lib_id+"/"+str(i)+"MOL.pdbqt"
        if os.path.isfile(out_file):
            continue
        my_file = tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False)
        my_name = my_file.name
        my_file.close()
       # Convert to PDB
        Chem.MolToPDBFile(mol, my_name)
        # Convert to PDBQT
        tmp_file = prepare_mol(my_name,"lig")
        if os.path.isfile(tmp_file):
            out_pdbqt = open(out_file, "w")
            my_text = open(tmp_file).read()
            out_pdbqt.write(my_text)
        else:
            print "ERROR IN MAKING MOL"
   # Save as pdbqt using the number
    # Return the uuid of the library
    return lib_id


def get_lib_ligs(lib_uuid):
    """Function to retrieve the ligands from a unique id"""
    my_lib = "/DockingLibs/"+lib_uuid
    # First make the library
    in_mols = [mol for mol in Chem.SDMolSupplier(str(my_lib+"/mols.sdf"))]
    for mol in in_mols:
        AllChem.EmbedMolecule(mol)
    save_lib(in_mols, my_lib)
    print my_lib
    out_ligs = [x for x in os.listdir(my_lib) if x.endswith(".pdbqt")]
    return sorted(out_ligs)


def run_docking(conf_dir, ligand, lig_dir, run_dir):
    if not os.path.exists(run_dir):
        os.mkdir(run_dir)
    my_result = run_dir + "/"+ligand.replace(".pdbqt", "_out.pdbqt")
    my_lig = lig_dir + "/" + ligand
    conf_path = os.path.join(conf_dir,[my_file for my_file in os.listdir(conf_dir) if my_file.endswith(".conf")][0])
    run_vina(conf_path, my_result, my_lig)
    # Parses the output and returns an SDFile for the composite result
    if not os.path.isfile(my_result):
        return []
    else:
        return parse_result(my_result)


def prepare_env(pdb_code, smiles=""):
    """Function to prepare a given docking environment"""
    # Get the pdb code
    me_d = get_pdb([{"pdb_code": pdb_code,"smiles":smiles}])
    # Get the files
    me = me_d.keys()[0]
    file_lines = open(me_d[me]["path"]).readlines()
    out_lig = extract_ligand(file_lines, me_d[me]["smiles"])
    apo_path = me_d[me]["path"].replace(".pdb", "APO.pdb")
    me_d[me]["apo_path"] = apo_path
    out_f = open(apo_path,"w")
    out_f.write(remove_hetatm(file_lines))
    for i, lig in enumerate(out_lig):
        lig_path = me_d[me]["path"].replace(".pdb", "LIG"+str(i))+".pdb"
        Chem.MolToPDBFile(lig, lig_path)
        if i == 0:
            me_d[me]["docking_lig"] = lig_path
    me_d[me]["dock_prot"] = prepare_mol(me_d[me]["apo_path"],"prot")
    me_d[me]["dock_lig"] = prepare_mol(me_d[me]["docking_lig"],"lig")
    # Now get the centre of mass 
    me_d[me]["x_com"], me_d[me]["y_com"], me_d[me]["z_com"] = centre_of_mass_from_rdmol(lig)
    # Print me_d
    conf_path = me_d[me]["path"].replace(".pdb", "VINA.conf")
    out_f = open(conf_path, "w")
    out_text = "\n".join(["receptor = "+me_d[me]["dock_prot"], "ligand = "+me_d[me]["dock_lig"], "center_x =  " + str(me_d[me]["x_com"]),  "center_y =  " + str(me_d[me]["y_com"]), "center_z =  " + str(me_d[me]["z_com"]), "size_x = 25", "size_y = 25", "size_z = 25"])
    out_f.write(out_text)
    out_f.close()
    return conf_path


if __name__ == "__main__":
    me_d = get_pdb([{"pdb_code":sys.argv[1],"smiles":""}])
    me = me_d.keys()[0]
    file_lines = open(me_d[me]["path"]).readlines()
    out_lig = extract_ligand(file_lines, me_d[me]["smiles"])
    apo_path = me_d[me]["path"].replace(".pdb", "APO.pdb")
    me_d[me]["apo_path"] = apo_path
    out_f = open(apo_path,"w")
    out_f.write(remove_hetatm(file_lines))
    for i, lig in enumerate(out_lig):
        lig_path = me_d[me]["path"].replace(".pdb", "LIG"+str(i))+".pdb"
        Chem.MolToPDBFile(lig, lig_path)
        if i == 0:
            me_d[me]["docking_lig"] = lig_path
    me_d[me]["dock_prot"] = prepare_mol(me_d[me]["apo_path"],"prot")
    me_d[me]["dock_lig"] = prepare_mol(me_d[me]["docking_lig"],"lig")
    # Now get the centre of mass 
    me_d[me]["x_com"], me_d[me]["y_com"], me_d[me]["z_com"] = centre_of_mass_from_rdmol(lig)
    # Print me_d
    conf_path = me_d[me]["path"].replace(".pdb", "VINA.conf")
    out_f = open(conf_path, "w")
    out_text = "\n".join(["receptor = "+me_d[me]["dock_prot"], "ligand = "+me_d[me]["dock_lig"], "center_x =  " + str(me_d[me]["x_com"]),  "center_y =  " + str(me_d[me]["y_com"]), "center_z =  " + str(me_d[me]["z_com"]), "size_x = 25", "size_y = 25", "size_z = 25"])
    out_f.write(out_text)
    out_f.close()
    run_docking(os.path.split(conf_path)[0], me_d[me]["dock_lig"])
