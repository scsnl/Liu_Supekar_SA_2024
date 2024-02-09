# please install abagen in terminal first
# by running the following:
# pip install brainsmash
# 

import numpy
import pandas as pd
from brainsmash.mapgen.base import Base

# input data
# brain map  a one-dimensional scalar vector (Nx1, N = number of regions)
brain_map_file = "Y:\\projects\\jinliu5\\2021_Longt_math_gene\\results\\smri\\vbm\\Stanford_cohort\\CCA_GMV_math_N219\\CCA_math_brainmap_N219.txt"  # use absolute paths if necessary
# brain_map_file = "//Volumes//menon//projects//jinliu5//2021_Longt_math_gene//results//smri//vbm//Stanford_cohort//CCA_GMV_math_N219//CCA_math_brainmap_N219.txt"  # use absolute paths if necessary
# distance matrix containing a measure of distance between each pair of elements in the brain map (NxN)
dist_mat_file = "C:\\Users\\jinliu5\\Box\\Jin Liu\\2021 Longt math gene\\0. ELN\\1. log of results\\data, scripts, and results\\scripts\\BrainSMASH\\Distance_BN246.txt"
num_permutation = 1000  # modified the number of permutation times (M)!

# generating surrogate maps (MxN)
base = Base(x=brain_map_file, D=dist_mat_file)
surrogates = base(n=num_permutation)
print(surrogates)
numpy.savetxt("surrogatebrainmap.txt",surrogates) 

# from brainsmash.mapgen.eval import base_fit
# from brainsmash.utils.eval import sampled_fit  analogous function for Sampled class
# base_fit(brain_map_file, dist_mat_file, nsurr=100)