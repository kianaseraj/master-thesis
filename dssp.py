from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import numpy as np
import os
from tqdm import tqdm
from os import walk

#extracting relative ASA with DSSP
p = PDBParser()
#pdb_dir = "/path/to/pdb files directory"
#out_dir = "/Path/t/output directory"

for (dirpath, dirnames, filenames) in walk(pdb_dir):
    break 
for word in filenames[:]:
  if not word.endswith("pdb"):
    filenames.remove(word)

for file in tqdm(filenames):
  # DSSP data is accessed by a tuple (chain_id, res_id)
  """
outputs provided by dssp are:
(dssp index, amino acid, secondary structure, relative ASA, phi, psi,
NH_O_1_relidx, NH_O_1_energy, O_NH_1_relidx, O_NH_1_energy,
NH_O_2_relidx, NH_O_2_energy, O_NH_2_relidx, O_NH_2_energy)
  """
  try:
    id = os.path.splitext(file)[0]
    structure = p.get_structure(f"{id}", f"{file}")
    model = structure[0]
    dssp = DSSP(model, f"{file}", dssp = "mkdssp")


    dssp_list = []
    pdb_list = []
    for a_key in dssp.keys():
        a = dssp[a_key]
        dssp_list.append(a)
    np.save(out_dir + f"/{id}", str(dssp_list))
   
  except Exception as e:
    pdb_list.append(f"{file}")
    print(f"Error parsing {file}: {e}")
    continue