# Takes average <S2> from many random seeds
# and combines them into one file

import os
import numpy as np

# Save all the file names in the path as strings to a list
path = '../Data/'
filenames_all = os.listdir(path)

# Only keep the files that contain the string 'SWAP' in their name,
# the desired "U" in their name, and that were able to collect all bins
string_to_keep = '_SWAP_'  # The string that kept files will contain
U_to_keep = '3.300000'
bins_wanted = 1000
files_SWAP = []
for filename in filenames_all:
    if (string_to_keep in filename) and (U_to_keep in filename):
        if os.stat(path+filename).st_size > 0:
            with open(path+filename) as f:
               count = sum(1 for _ in f)
            if count == bins_wanted:
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

# Calculate S2 of each bin
SWAP_0 = combined_SWAP_data[:,0]
S2_data = -np.log(combined_SWAP_data / SWAP_0[:,None])

# Get mean and std dev,err of S2
S2_mean = np.mean(S2_data,axis=0)
S2_stderr = np.std(S2_data,axis=0)/np.sqrt(number_of_seeds)

# Print out <S2> +/- error
for l in range(columns_per_file):
    print(f"<S2(l={l})> = {S2_mean[l]} +/- {S2_stderr[l]:0.8f}")
    
print(number_of_seeds)
