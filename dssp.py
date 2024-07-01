from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import numpy as np
import os
from tqdm import tqdm
from os import walk

"""extracting relative accessible surface area(ASA) of each residue using Bio.PDB.DSSP module"""


#pdb_dir = "/path/to/pdb files directory which their ASA will be extracted" 
#out_dir = "/Path/to/output directory to save the result"



#only keeping files with pdb exatension for further works.
for (dirpath, dirnames, filenames) in walk(pdb_dir):
    break 
for word in filenames[:]:
  if not word.endswith("pdb"):
    filenames.remove(word)



p = PDBParser() #instantiating a Parser for PDB files!

#generate and save the secondary structure and accessibility values for each pdb file with DSSP program!
for file in tqdm(filenames): 
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