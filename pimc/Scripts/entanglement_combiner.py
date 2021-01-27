# Takes average <S2> from many random seeds
# and combines them into one file

import os
import numpy as np

incomplete_seeds = []
seeds_list = list(range(1000))
seeds_measured = []

# Save all the file names in the path as strings to a list
path = path="/Users/ecasiano/Desktop/PrototypeScripts/Data/"
filenames_all = os.listdir(path)

# Saves the files relevant to the Renyi Entanglement Entropy calculation
files_SWAP = []

# Iterate over all filenames
for filename in filenames_all:
    
    # Extract parameter information from file name
    parameters = filename.split("_")
    
    if parameters[0]=='1D' or parameters[0]=='2D' or parameters[0]=='3D':
    
        D = int((parameters[0])[0]) # hypercube dimension
        L = int(parameters[1]) # hypercube linear size
        N = int(parameters[2]) # total particles
        l = int(parameters[3]) # subsystem linear size
        U = float(parameters[4]) # interaction potential
        t = float(parameters[5]) # tunneling parameter
        beta = float(parameters[6]) # imaginary time length (K_B*T)**(-1)
        bins_wanted = int(parameters[7]) # number of bins saved in file
        filetype = parameters[8] # identifies the data stored in file
        seed = int(parameters[9].split(".")[0]) # random seed used

        if filetype=='SWAP':
            # Set parameters of simulations from differenet seeds we want to evaluate [D,L,N,l,U,t,beta,bins,type]
            parameters_to_evaluate = [1,
                                      4,
                                      4,
                                      2,
                                      3.300000,
                                      1.000000,
                                      4.000000,
                                      10000,
                                      'SWAP']

            if [D,L,N,l,U,t,beta,bins_wanted,filetype] == parameters_to_evaluate:
                if os.stat(path+filename).st_size > 0:
                    with open(path+filename) as f:
                       count = sum(1 for _ in f)
                    if count == bins_wanted: # only consider files that managed to save number of bins wanted
                        files_SWAP.append(filename)
                        seeds_measured.append(seed)
                    else:
                        incomplete_seeds.append(seed)
        
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
    print(f"<S2(l={l})> = {S2_mean[l]:0.8f} +/- {S2_stderr[l]:0.8f}")
    
print(number_of_seeds)
print("incomplete seeds: ",[int(i) for i in incomplete_seeds])

#seeds_measured.sort()
#print([x for x in range(seeds_measured[0],seeds_measured[-1]+1) if x not in seeds_measured])
