#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 14:14:22 2022

@author: Alexander G. Lucaci

@Description: AlignmentProfiler provides quick summary statistics for multiple
                 sequence alignments.
 
 **Current Features**
 
 *   Get nucleotide frequencies
 
 **Future development**
 *   Get site likelihood statistics
 *   Get GC Content
 *   Get Amino acid relative frequencies
 *   Get essential versus non-essential amino acid statistics
 *   Get ungapped nucleotide content
 *   Bootstrap an alignment
 *   Test for ambiguous characters
 *   Get gap character analysis
 *   Get codon diversity
 *   Get protein diversity
 *   Generate Hamming distance plot
 *   Get Jukes-Cantor distances
 *   Get Relative Synonymous Codon Usage (RSCU)


"""

# =============================================================================
# Imports
# =============================================================================

import pandas as pd
from pandas.plotting import scatter_matrix
import random
from matplotlib import pyplot as plt
import itertools
import numpy as np
from matplotlib.pyplot import figure
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np 
from pandas import DataFrame
import seaborn as sns
from tqdm import tqdm
import argparse
import os
import sys
import pathlib 

# =============================================================================
# Declares
# =============================================================================

parser = argparse.ArgumentParser(description='Process some integers.')

parser.add_argument('--input', metavar='I', type=str,
                    help='Input multiple sequence alignment')

args = parser.parse_args()

# =============================================================================
# Classes
# =============================================================================

class AlignmentProfiler(): # Loads in a fasta dictionary
  def __init__(self, _Alignment, MoleculeType="DNA"):
    
    File_Extension = pathlib.Path(_Alignment).suffix
    if File_Extension == ".nex":
        File_Type = "nexus"
    else:
        File_Type = "fasta"
    #end if  
    
    #self.records = self.process_sites(FASTA_FILE)
    #self.CodonTable = CodonTable
    self.NumSites = self.Get_NumSites(_Alignment, File_Type)
    self.NumSequences = self.Get_NumSequences(_Alignment, File_Type)
    #self.NumSequences = self.Get_NumSequences()
    #self.alignment_type = "DNA" # Default
    self.MoleculeType = MoleculeType
    self.GapCharacter = "-"                    # Default
    self.DNA_Characters = ["T", "C", "G", "A"] # Default
    
    self.NumInvariantSites = Get_NumInvariantSites(_Alignment, File_Type)
  #end method

  def Get_NumSites(self, _Alignment, File_Type):
      sites = 0
      with open(_Alignment, "r") as handle:
          for record in SeqIO.parse(handle, File_Type):
              if len(str(record.seq)) % 3 == 0:
                  sites = int(len(str(record.seq)) / 3)
                  return sites
              else:
                  print("# ERROR: Number of sites is not a multiple of 3")
            #end if
          #end for
          
      #end with
      #return sites, sequences  
  #end method

  def Get_NumSequences(self, _Alignment, File_Type):
      sequences = 0
      with open(_Alignment, "r") as handle:
          for record in SeqIO.parse(handle, File_Type):
              sequences += 1
          #end for
      #end with
      return sequences  
  #end method
  
  """
  def Get_NumInvariantSites(self, _Alignment, File_Type):
      InvariantSites = 0
      
      for i in range(self.NumSites):
          column_data = []
          
          with open(_Alignment, "r") as handle:
              for record in SeqIO.parse(handle, File_Type):
                  if len(str(record.seq)) % 3 == 0:
                      sites = int(len(str(record.seq)) / 3)
                      return sites
                  else:
                      print("# ERROR: Number of sites is not a multiple of 3")
                #end if
              #end for
          #end with
    
  #end method
  """


#end class

# =============================================================================
# Helper functions
# =============================================================================

# =============================================================================
# Main
# =============================================================================

print("# Starting to profile multiple sequence alignment:", args.input)

if os.path.exists(args.input) and os.path.getsize(args.input) > 0:
    _MSA = AlignmentProfiler(args.input)
    print("# Loading complete on an alignment with", _MSA.NumSequences, "sequences, and", _MSA.NumSites, "sites")
else:
    print("# ERROR: the file is either empty or does not exist")
#end if

# Invariant Sites





"""
# Summary stats
print("# Loaded an alignment from:", input_file)
print("# This alignment contains", TEST.num_sequences, "sequences")
print("# This alignment contains", TEST.num_sites, "sites")
print("# Number of invariant sites in alignment:", len(TEST.invariant_sites()))
print("# Fraction of invariant sites in alignment:", len(TEST.invariant_sites())/TEST.num_sites)
gap_list = TEST.gaps_distribution()
avg_gap = sum(gap_list)/ len(gap_list)
print("# Average measure of gappiness in alignment is:", avg_gap)
print("# Maximum measure of gappiness in alignment at a particular site is:", max(gap_list))
print("# Minimum measure of gappiness in alignment at a particular site is:", min(gap_list))

"""


# =============================================================================
# End of file
# =============================================================================
