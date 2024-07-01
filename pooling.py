
import torch
from typing import List
from abc import abstractmethod, ABC
import numpy as np
import ast


"""Define pooling object consisting of pooling methods : mean, max, ASA, etc."""

class pooling(ABC):
  def __init__(self):
    pass
  
  @abstractmethod
  def embedding_pooling(self):
    pass
  
    
def dssp_get(dssp_file):
    """read dssp file"""
    prot_id = dssp_file.split("/")[-1].split(".")[0]
    dssp_output = np.load(dssp_file)
    dssp_output = dssp_output.item()
    dssp_output = ast.literal_eval(dssp_output)
    return dssp_output, prot_id

def accessible_residues(dssp_output) -> List:
    """based on the dssp values, residues having less than 0.3 relative accessible surface area will be eliminated from the list of indices"""
    surface_accessible_aa = []
    for residue_index, dssp_val in enumerate(dssp_output):
        if float(dssp_val[3]) >= 0.3:
            surface_accessible_aa.append(residue_index)
    return surface_accessible_aa

def sequence_embedding(embedding_file, sequence_file):
  """the esm token representation inserts 0 at the beginning and 2 at the end of each sequence.Additionally, to pad sequences to the max sequence length, number 1 will be added at the end after token id 2!
  example : "ATIDPK"=[0, 5, 11, 12, 13, 14, 15, 2, 1, 1, 1, 1]
  Hence, check if the embedding's length matches the length of its sequence, otherwise, eliminated the padding tokens!""" 
  embedding = torch.load(embedding_file)
  with open(sequence_file, "r") as f:
    sequence = f.read()
  embedding = embedding[0][1:len(sequence)+1]
  return embedding


class ASA_pooling():
  def __init__(self):
    super().__init__()
    
  def embedding_pooling(dssp_file, embedding_file : str, sequence_file : str, save_path : str) ->torch.Tensor:
    """
    Amino acids having less than 30% of relative ASA will be eliminated from feature embeddings.
    """
    dssp_output, prot_id = dssp_get(dssp_file)
    surface_accessible_aa = accessible_residues(dssp_output)
    embedding = sequence_embedding(embedding_file, sequence_file)
    embedding = embedding[surface_accessible_aa]
    torch.save(embedding, save_path+"/"+f"{prot_id}.pt")

    