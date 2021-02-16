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

incomplete_seeds = [] 
seeds_list = list(range(1000))
seeds_measured = []

# Save all the file names in the path as strings to a list
path = path="/Users/ecasiano/Desktop/PrototypeScripts/TestSystematicErrorv4/"
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
for mA_sector_wanted in range(l_want**D_want//2,l_want**D_want//2+1):

    # Saves the files relevant to P(n) & S2(n) calculation
    files_PnSquared = []
    files_SWAPn = []
    files_Pn = []

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
    print('Initial number of seeds: ',number_of_seeds,"\n")
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
        
    # Jacknife 
    SWAP_col_seed_sum = np.sum(SWAP_col_sums,axis=0)
    PnSquared_col_seed_sum = np.sum(PnSquared_col_sums,axis=0)
    Pn_col_seed_sum = np.sum(Pn_col_sums,axis=0)
    S2_i = np.zeros(SWAP_col_sums.shape)
    S2_acc_i = np.zeros(SWAP_col_sums.shape[0])
    for i in range(SWAP_col_sums.shape[0]):
        
        SWAP_jacknifed_sum = SWAP_col_seed_sum - SWAP_col_sums[i]
        PnSquared_jacknifed_sum = PnSquared_col_seed_sum - PnSquared_col_sums[i]
        Pn_jacknifed_sum = Pn_col_seed_sum - Pn_col_sums[i]

        S2_i[i] = -np.log ( SWAP_jacknifed_sum / PnSquared_jacknifed_sum ) 
        S2_acc_i[i] = np.sum(Pn_jacknifed_sum * S2_i[i]) / np.sum(Pn_jacknifed_sum)

    S2 = S2_i 
    
    # Get Pn
    Pn = Pn_col_sums / np.sum(Pn_col_sums,axis=1)[:,None]

    # Some seeds might have been thrown away
#     number_of_seeds = SWAP.shape[0]

    print("\nFinal number of seeds: ",number_of_seeds)

    S2_mean = np.mean(S2,axis=0)
    S2_err = np.std(S2,axis=0) * np.sqrt(S2.shape[0])

    S2_acc_mean = np.mean(S2_acc_i,axis=0)
    S2_acc_err = np.std(S2_acc_i,axis=0) * np.sqrt(S2_acc_i.shape[0])

    # Get mean of Pn and std errors
    Pn_mean = np.mean(Pn,axis=0)
    Pn_err = np.std(Pn,axis=0) / np.sqrt(number_of_seeds)

    print("\npartition size: ", mA_sector_wanted,"\n")
    # Print P(n),S2(n) for each sector to screen
    for i in range(len(S2_mean)):
        print("P(n=%d) = %.4f +/- %.4f"%(i,Pn_mean[i],Pn_err[i]))
    print("\n")
    for i in range(len(S2_mean)):
        print("S2(n=%d) = %.4f +/- %.4f"%(i,S2_mean[i],S2_err[i]))
        
    print("\nS2_acc = %.4f +/- %.4f"%(S2_acc_mean,S2_acc_err))
    
print("\nincomplete seeds: ",[int(i) for i in incomplete_seeds])

# print("Seeds not included for some reason: ")

# for i in seeds_list:
#     if not(i in seeds_measured):
#         print(i,end=",")