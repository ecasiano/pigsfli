# Takes average <S2> from many random seeds
# and combines them into one file

import os
import numpy as np

# Save all the file names in the path as strings to a list
path = 'Data/'
filenames_all = os.listdir(path)

# Only keep the files that contain the string 'SWAP' in their name
string_to_keep = 'SWAP'  # The string that kept files will contain
files_SWAP = []
for filename in filenames_all:
    if string_to_keep in filename: files_SWAP.append(filename)
        
# Sort SWAP files in ascending order of random seed
files_SWAP.sort()

# SHOULD I COMPUTE <S2> before combining files?
# or SHOULD I WRITE total number of counts for each SWAP-sector?

# Will continue later, after asking Adrian and Chris....

# Initialize arrays to save seed and average S2
N_elements = np.shape(files_egs)[0]
v_list = np.ones(N_elements)
egs_list = np.ones(N_elements)