# Takes average <S2> from many random seeds
# and combines them into one file

# What we're doing on one file:
# 1) Loading P(n) for desired swapped sector mA
# 2) Load SWAP(n) files
# 3) Taking column sum for P(n) file for desired mA 
# 4) Taking column sum for SWAP(n) file
# 5) Dividing column sums of P(n) over column sums of SWAP(n)
# 6) Taking negative log of above result
# 7) Taking the average over seeds of above results
# 8) Take the standard error of above average

import os
import numpy as np

x = np.array([0,1,-1])
y = np.array([0,0,0])
z = x/y
print(z)
print(np.nan_to_num(z,nan=0,posinf=1000,neginf=-1000))
print(z)

incomplete_seeds = [] 
seeds_list = list(range(1000))
seeds_measured = []

# Save all the file names in the path as strings to a list
path = path="/Users/ecasiano/Desktop/PrototypeScripts/TestSystematicErrorv3/"
filenames_all = os.listdir(path)

# Set desired total number of particles
L_want = 6
N_want = 6
l_want = 6
beta_want = 2.000000
bin_size_want = 10000
bins_want = 1000
D_want = 1
U_want = 3.300000
t_want = 1.000000

mA_sector_wanted = l_want**D_want//2

# Saves the files relevant to P(n) & S2(n) calculation
files_PnSquared = []
files_SWAPn = []
files_Pn = []

seed_min = 285
seed_max = 285
# Iterate over all filenames in path
for filename in filenames_all:

    # Extract parameter information from file name
    parameters = filename.split("_")
    
    if parameters[0]=='1D' or parameters[0]=='2D' or parameters[0]=='3D':
    
        D = int((parameters[0])[0]) # hypercube dimension
        L = int(parameters[1]) # hypercube linear size
        N = int(parameters[2]) # total particles
        l = int(parameters[3]) # subsystem linear size (actually l_max)
        U = float(parameters[4]) # interaction potential
        t = float(parameters[5]) # tunneling parameter
        beta = float(parameters[6]) # imaginary time length (K_B*T)**(-1)
        bin_size = int(parameters[7])
        bins_wanted = int(parameters[8]) # number of bins saved in file
        filetype = (parameters[9]).split("-mA") # identifies the data stored in file
        seed = int(parameters[10].split(".")[0]) # random seed used
                    
        if filetype[0]=='PnSquared' and int(filetype[1])==mA_sector_wanted:
            # Set parameters of simulations from differenet seeds we want to evaluate [D,L,N,l,U,t,beta,bins,type]
            parameters_to_evaluate = [D_want,
                                      L_want,
                                      N_want,
                                      l_want,
                                      U_want,
                                      t_want,
                                      beta_want,
                                      bin_size_want,
                                      bins_want,
                                      'PnSquared']

            if [D,L,N,l,U,t,beta,bin_size,bins_wanted,filetype[0]] == parameters_to_evaluate:
                if os.stat(path+filename).st_size > 0:
                    with open(path+filename) as f:
                       count = sum(1 for _ in f)
                    if count == bins_wanted: # only consider files that managed to save number of bins wanted
                        files_PnSquared.append(filename)
                        
                        filename_splitted = filename.split('_')
                        filename_splitted[9] = 'Pn-mA'+str(mA_sector_wanted)
                        filename_Pn = "_".join(filename_splitted)
                        files_Pn.append(filename_Pn)
                        
                        filename_splitted = filename.split('_')
                        filename_splitted[9] = 'SWAPn-mA'+str(mA_sector_wanted)
                        filename_SWAPn = "_".join(filename_splitted)
                        files_SWAPn.append(filename_SWAPn)
                            
                        seeds_measured.append(seed)
                            
                    else: 
                        incomplete_seeds.append(seed)
                        
# Get total number of seeds 
number_of_seeds = len(files_PnSquared)//1

# Get column sum of P(n) files for each seed
print('Initial number of seeds: ',number_of_seeds)
PnSquared_col_sums = np.zeros((number_of_seeds,N_want+1)).astype(int)
Pn_col_sums = np.zeros((number_of_seeds,N_want+1)).astype(int)
SWAP_col_sums = np.zeros((number_of_seeds,N_want+1)).astype(int)
for s in range(number_of_seeds):

    data = np.loadtxt(path+files_PnSquared[s])
    PnSquared_col_sums[s] = np.sum(data,axis=0)
 
    data = np.loadtxt(path+files_Pn[s])
    Pn_col_sums[s] = np.sum(data,axis=0)
    
    data = np.loadtxt(path+files_SWAPn[s])
    SWAP_col_sums[s] = np.sum(data,axis=0)
#     print(np.mean(data,axis=1))
    
# print("\n\n")
# print(SWAP_col_sums.astype(int))
# print("----------------------------------")
# print(PnSquared_col_sums.astype(int))

# print(SWAP_col_sums[:,0])
# print(SWAP_col_sums[:,-1])


# Get invalid row indices where there is a zero.
# This avoids division by zero and/or log(0)=inf errors.
# invalid_rows = []
# for i in range(number_of_seeds):
#     if 0 in SWAP_col_sums[i]: 
#         invalid_rows.append(i)
#         print("!!!!!")
#     if 0 in PnSquared_col_sums[i]: 
#         invalid_rows.append(i)
#         print("?????") # IF THIS PRINTS, THE IF ABOVE SHOULD HAVE PRINTED TOO
        
#!!! Actually, the statement above might not necessarily be true. Think about it.

# Try doing a slightly longer run, where we add to the measurement counter only after
# the histograms corresponding to \ell=2 have a certain number of bins (try with 1000)

# invalid_rows = np.unique(np.array(invalid_rows)).astype(int)

# Eliminate invalid rows
# SWAP_col_sums = np.delete(SWAP_col_sums,invalid_rows,axis=0)
# PnSquared_col_sums = np.delete(PnSquared_col_sums,invalid_rows,axis=0)


# Calculate expectation value of SWAP operator for each n-sector & seed
SWAP = SWAP_col_sums/PnSquared_col_sums

# The ratio above calculated some 0/0 = nan elements and some num/0 elements
SWAP = np.nan_to_num(SWAP,nan=1.0,posinf=1.0)

# Calculate Second Renyi Entanglement Entropy for each n-sector & seed
S2 = -np.log(SWAP)

# The log above calculated some log(0)=-inf elements
S2 = np.nan_to_num(S2,neginf=0,posinf=0)

# Get Pn
Pn = Pn_col_sums / np.sum(Pn_col_sums,axis=1)[:,None]

# Some seeds might have been thrown away
number_of_seeds = SWAP.shape[0]

print("\nFinal num_seeds: ",number_of_seeds)

# Get mean of S2 and std errors
S2_mean = np.mean(S2,axis=0)
S2_err = np.std(S2,axis=0) / np.sqrt(number_of_seeds)

# Get mean of Pn and std errors
Pn_mean = np.mean(Pn,axis=0)
Pn_err = np.std(Pn,axis=0) / np.sqrt(number_of_seeds)

# Print P(n),S2(n) for each sector to screen
for i in range(len(S2_mean)):
    print("P(n=%d) = %.6f +/- %.6f"%(i,Pn_mean[i],Pn_err[i]))
print("\n")
for i in range(len(S2_mean)):
    print("S2(n=%d) = %.6f +/- %.6f"%(i,S2_mean[i],S2_err[i]))

# print("\n")
# for i in range(len(S2_mean)):
#     print("Numerator(n=%d) = %.6f"%(i,np.mean(SWAP_col_sums[i],axis=0)))

# print("\n")
# for i in range(len(S2_mean)):
#     print("Denominator(n=%d) = %.6f"%(i,np.mean(PnSquared_col_sums[i],axis=0)))
    
print("incomplete seeds: ",[int(i) for i in incomplete_seeds])

# print("Seeds not included for some reason: ")

# for i in seeds_list:
#     if not(i in seeds_measured):
#         print(i,end=",")
