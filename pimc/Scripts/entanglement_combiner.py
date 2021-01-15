# Takes average <S2> from many random seeds
# and combines them into one file

import os
import numpy as np

# Save all the file names in the path as strings to a list
path = '../'
filenames_all = os.listdir(path)

# Only keep the files that contain the string 'SWAP' in their name
string_to_keep = '_SWAP_'  # The string that kept files will contain
U_to_keep = '10.000000'
files_SWAP = []
for filename in filenames_all:
    if (string_to_keep in filename) and (U_to_keep in filename):
        if os.stat(path+filename).st_size > 0:
            files_SWAP.append(filename)
        
# Sort SWAP files in ascending order of random seed
files_SWAP.sort()

# Get total number of seeds and total columns per data file
reference_file = np.loadtxt(path+files_SWAP[0])
number_of_seeds = len(files_SWAP)
columns_per_file = reference_file.shape[1]

combined_SWAP_data = np.zeros((number_of_seeds,columns_per_file))
for i,filename in enumerate(files_SWAP):
    data = np.loadtxt(path+filename)
    data_mean = np.mean(data,axis=0)
    combined_SWAP_data[i] = data_mean

# Get mean and std dev of each column (swapped region size)
SWAP_mean = np.mean(combined_SWAP_data,axis=0)
SWAP_var = np.var(combined_SWAP_data,axis=0)

# Calculate S2
SWAP_0 = SWAP_mean[0]
S2 = -np.log(SWAP_mean / SWAP_0)

# Determine error propagation of S2

# Print out <S2> +/- error
print("<S2(l)>: ",S2)
