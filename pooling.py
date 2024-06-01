import torch
import torch.nn as nn
from constants import *
from torch.autograd import Variable
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB import PDBParser
from dataclasses import dataclass
from abc import abstractmethod, ABC
import numpy as np
import os
import ast

#Amino acids having less than 30% of relative ASA will be eliminated from feature embeddings.

"""
based on the dssp putput:
(dssp index, amino acid, secondary structure, relative ASA, phi, psi,
NH_O_1_relidx, NH_O_1_energy, O_NH_1_relidx, O_NH_1_energy,
NH_O_2_relidx, NH_O_2_energy, O_NH_2_relidx, O_NH_2_energy)
"""
dir_path = "/Users/kianaseraj/Desktop/github-kianaseraj/master-thesis/dataset/dssp_go"
os.chdir(dir_path)
def dssp_get(dssp_file):
    prot_id = os.path.splitext(dssp_file)[0]
    dssp_output = np.load(dssp_file)
    dssp_output = dssp_output.item()
    dssp_output = ast.literal_eval(dssp_output)
    return dssp_output, prot_id

def accessible_residues(dssp_output):
    
    surface_accessible_aa = []
    for residue_index, dssp_val in enumerate(dssp_output):
        if dssp_val[3] >= 0.3:
            surface_accessible_aa.append(residue_index)
    return surface_accessible_aa

#storing all sequences with their surface-accessible residues in a dict format
ASA_dict = {}
for file in os.listdir(dir_path):
    dssp_output, prot_id = dssp_get(file)
    surface_accessible_aa = accessible_residues(dssp_output)
    ASA_dict[prot_id] = surface_accessible_aa
    

    

        
    


    