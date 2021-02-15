# Takes average <S2> from many random seeds
# and combines them into one file

import os
import numpy as np

incomplete_seeds = []
seeds_list = list(range(1000))
seeds_measured = []

# Append sweep results to same list so we can copy paste to plotting script
all_S2_results = []
all_S2_errors = []

# Save all the file names in the path as strings to a list
path = path="/Users/ecasiano/Desktop/PrototypeScripts/BetaScaling/"
filenames_all = os.listdir(path)

for BETA in [0.6 , 0.7, 0.8, 0.9, 1., 1.15, 1.3,1.5, 1.75, 2.  , 3.  , 4.  , 5.  , 6. ]:
    print("BETA=",BETA)
    # Set desired total number of particles
    L_want = 4
    N_want = 4
    l_want = 2
    beta_want = BETA
    bin_size_want = 10000
    bins_want = 1000
    D_want = 1
    U_want = 10.0 # For U=0.1, there is one beta that I forgot to simulate. Figure out which one it was then run it.
    t_want = 1.000000

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
            l = int(parameters[3]) # subsystem linear size (actually l_max)
            U = float(parameters[4]) # interaction potential
            t = float(parameters[5]) # tunneling parameter
            beta = float(parameters[6]) # imaginary time length (K_B*T)**(-1)
            bin_size = int(parameters[7])
            bins_wanted = int(parameters[8]) # number of bins saved in file
            filetype = (parameters[9]) # identifies the data stored in file
            seed = int(parameters[10].split(".")[0]) # random seed used

            if filetype=='SWAP':

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
                                          'SWAP']

                if [D,L,N,l,U,t,beta,bin_size,bins_wanted,filetype] == parameters_to_evaluate:
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
        
    # Save the l=2 results
    all_S2_results.append(S2_mean[l_want])
    all_S2_errors.append(S2_stderr[l_want])

    print(number_of_seeds)
    print("incomplete seeds: ",[int(i) for i in incomplete_seeds])

print("\n\n")
for result in all_S2_results:
    print(f"{result:0.8f}",end=",")
print("\n\n")
for error in all_S2_errors:
    print(f"{error:0.8f}",end=",")

#seeds_measured.sort()
#print([x for x in range(seeds_measured[0],seeds_measured[-1]+1) if x not in seeds_measured])
