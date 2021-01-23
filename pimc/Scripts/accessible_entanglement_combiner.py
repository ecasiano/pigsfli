# Takes average <S2> from many random seeds
# and combines them into one file

import os
import numpy as np

# 1) Get file names of P(n) data corresponding to desired m_A file
# 2) While you're at it, get file names of corresponding SWAPn files
# 3) Maybe sort them such that random seeds match

# 4) Finish the above steps first...

# Save all the file names in the path as strings to a list
path = '/Users/ecasiano/XCode/pimc/pimc/Data/'
filenames_all = os.listdir(path)

N = 4
# Saves the files relevant to the Renyi Entanglement Entropy calculation
files_Pn = []
files_SWAPn = []
for i in range(N+1):
    files_SWAPn.append([])

# Iterate over all filenames
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
                                      4,
                                      N,
                                      2,
                                      3.300000,
                                      1.000000,
                                      4.000000,
                                      1000,
                                      'PnSquared']

            if [D,L,N,l,U,t,beta,bins_wanted,filetype[0]] == parameters_to_evaluate:
                if os.stat(path+filename).st_size > 0:
                    with open(path+filename) as f:
                       count = sum(1 for _ in f)
                    if count == bins_wanted: # only consider files that managed to save number of bins wanted
                        files_Pn.append(filename)
                        
                        for i in range(N+1):
                            filename_splitted = filename.split('_')
                            filename_splitted[8] = 'SWAP-n'+str(i)
                            filename_SWAPn = "_".join(filename_splitted)
                            
                            files_SWAPn[i].append(filename_SWAPn)

# Get total number of seeds and total columns per data file
number_of_seeds = len(files_Pn)

# Stores the mean Pn for each Pn-mA file
Pn_raw = np.zeros((number_of_seeds,N+1))
SWAPn_raw = np.zeros((number_of_seeds,N+1))
for i,filename in enumerate(files_Pn):
    data_Pn = np.loadtxt(path+filename)
    
#     data_Pn_norms = np.sum(data_Pn,axis=1) # bin norms
        
#     # Normalize (some bins might be all zeros. Don't divide these)
#     data_Pn /= data_Pn_norms[:,None]
#     if 0 in data_Pn_norms:
#         x = np.where(data_Pn_norms==0)
#         data_Pn[x] = 0
        
    Pn_raw[i] = np.mean(data_Pn,axis=0)
    Pn_raw[i] /= np.sum(Pn_raw[i]) 
    
    for j in range(len(files_SWAPn)): # iterates over N+1 list elements
        for k,file_SWAPn in enumerate(files_SWAPn[j]): # iterates over number_of_seeds list elements
            
            data_SWAPn_raw = np.loadtxt(path+file_SWAPn)
            data_SWAPn_mean = np.mean(data_SWAPn_raw,axis=0) # now we only have one 1D array
            
            data_SWAPn_unswapped_prob = data_SWAPn_mean[0] # mean prob of no swap
            
            if data_SWAPn_unswapped_prob != 0:
                SWAPn_mean = -np.log(data_SWAPn_mean/data_SWAPn_unswapped_prob)
            else:
                SWAPn_mean = np.zeros(mA_sector_wanted+1)
            
            SWAPn_raw[k][j] = SWAPn_mean[mA_sector_wanted]
            
# combined_SWAP_data = np.zeros((number_of_seeds,columns_per_file))
# for i,filename in enumerate(files_SWAP):
#     data = np.loadtxt(path+filename)
#     data_mean = np.mean(data,axis=0)
#     combined_SWAP_data[i] = data_mean

# # Calculate S2 of each bin
# SWAP_0 = combined_SWAP_data[:,0]
# S2_data = -np.log(combined_SWAP_data / SWAP_0[:,None])

# # Get mean and std dev,err of S2
# S2_mean = np.mean(S2_data,axis=0)
# S2_stderr = np.std(S2_data,axis=0)/np.sqrt(number_of_seeds)

# # Print out <S2> +/- error
# for l in range(columns_per_file):
#     print(f"<S2(l={l})> = {S2_mean[l]:0.8f} +/- {S2_stderr[l]:0.8f}")
    
# print(number_of_seeds)
