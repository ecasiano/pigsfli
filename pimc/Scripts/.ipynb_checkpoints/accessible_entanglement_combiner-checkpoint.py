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

# Save all the file names in the path as strings to a list
path = path="/Users/ecasiano/Desktop/PrototypeScripts/Data/"
filenames_all = os.listdir(path)

# Set desired total number of particles
L_want = 4
N_want = 4
l_want = 2
beta_want = 4.000000
bins_want = 10000

# Saves the files relevant to P(n) & S2(n) calculation
files_PnSquared = []
files_SWAPn = []
files_Pn = []
for i in range(N_want+1):
    files_SWAPn.append([])

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
        bins_wanted = int(parameters[7]) # number of bins saved in file
        filetype = (parameters[8]).split("-mA") # identifies the data stored in file
        seed = int(parameters[9].split(".")[0]) # random seed used
        
        mA_sector_wanted = 2
        
        if filetype[0]=='PnSquared' and int(filetype[1])==mA_sector_wanted:
            
            # Set parameters of simulations from differenet seeds we want to evaluate [D,L,N,l,U,t,beta,bins,type]
            parameters_to_evaluate = [1,
                                      L_want,
                                      N_want,
                                      l_want,
                                      3.300000,
                                      1.000000,
                                      beta_want,
                                      bins_want,
                                      'PnSquared']

            if [D,L,N,l,U,t,beta,bins_wanted,filetype[0]] == parameters_to_evaluate:
                if os.stat(path+filename).st_size > 0:
                    with open(path+filename) as f:
                       count = sum(1 for _ in f)
                    if count == bins_wanted: # only consider files that managed to save number of bins wanted
                        files_PnSquared.append(filename)
                        
                        for i in range(N+1):
                            filename_splitted = filename.split('_')
                            filename_splitted[8] = 'SWAP-n'+str(i)
                            filename_SWAPn = "_".join(filename_splitted)
                            
                            files_SWAPn[i].append(filename_SWAPn)
                        
                        filename_splitted = filename.split('_')
                        filename_splitted[8] = 'Pn-mA'+str(mA_sector_wanted)
                        filename_Pn = "_".join(filename_splitted)
                        files_Pn.append(filename_Pn)

# Get total number of seeds 
number_of_seeds = len(files_PnSquared)
            
# Get column sum of P(n) files for each seed
print('number of seeds: ',number_of_seeds)
Pn_squared_l_col_sums = np.zeros((number_of_seeds,N_want+1))
Pn_l_col_sums = np.zeros((number_of_seeds,N_want+1))
for s in range(number_of_seeds):
    data = np.loadtxt(path+files_PnSquared[s])
    Pn_squared_l_col_sums[s] = np.mean(data,axis=0)
    
    data = np.loadtxt(path+files_Pn[s])
    Pn_l_col_sums[s] = np.mean(data,axis=0)
       
SWAP_col_sums = np.zeros((number_of_seeds,N_want+1))
for s in range(number_of_seeds):
    for n in range(N_want+1):
        data = np.loadtxt(path+files_SWAPn[n][s])
        data = data[:,mA_sector_wanted]
        SWAP_col_sums[s][n] = (np.mean(data))

# Calculate expectation value of SWAP operator for each n-sector & seed
SWAP = SWAP_col_sums/Pn_squared_l_col_sums
invalid_rows = np.where(SWAP==0)[0] # Get rid of rows with SWAP(n)=0 to avoid log(0)=inf error
SWAP = np.delete(SWAP, invalid_rows, axis=0)

# Calculate Second Renyi Entanglement Entropy for each n-sector & seed
S2 = -np.log(SWAP)

# Get Pn
Pn = Pn_l_col_sums / np.sum(Pn_l_col_sums,axis=1)[:,None]

# Some seeds might have been thrown away
number_of_seeds = SWAP.shape[0]

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
    
print("Final num_seeds: ",SWAP.shape)
