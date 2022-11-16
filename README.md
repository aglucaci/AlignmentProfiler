# AlignmentProfiler

A package to get summary statistics and diagnostics for a multiple sequence alignment


## Test run

Example code 

```python3 src/AlignmentProfiler/AlignmentProfiler.py --input data/Benchmark/adh.nex --output results/adh.nex.csv```

Example output 

```
# ==============================================================================
# Starting to profile multiple sequence alignment: data/Benchmark/adh.nex
# ==============================================================================

# Checking alignment for invariant sites
100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 254/254 [00:01<00:00, 129.00it/s]
# Checking NT Frequencies
23it [00:00, 1686.99it/s]

# Loading complete on an alignment with 23 sequences, and 254 codon sites

# ==============================================================================
# Reporting alignment statistics
# ==============================================================================

# The alignment has 27 invariant sites ( 10.62992125984252 % ) 
# The alignment the following nucleotide frequencies...

Nucleotide            Frequency
------------------  -----------
Adenine (A)            0.233253
Thymine (T)            0.227091
Guanine (G)            0.253509
Cytosine (C)           0.286146
Any Nucleotide (N)     0

```




