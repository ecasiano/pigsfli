# Takes average <S2> from many random seeds
# and combines them into one file

import os
import numpy as np

# Set values of U sweep
# U_list = np.round(np.geomspace(0.01,100,20),4)
U_list = np.array([3.3])
# U_list = np.array([0.500000,
# 0.730000,
# 1.065800,
# 1.556100, 
# 2.272000, 
# 3.300000, 
# 4.843100, 
# 7.071100,   
# 10.323900,
# 16.666667, 
# 22.007100, 
# 32.130800,
# 46.911700,
# 68.492100,
# 100.000000])

beta_list = [0.6,0.7,0.8,0.9,1.0,1.15,1.30,1.50,
             1.75,2.0,2.25,2.50,2.75,3.0,3.25,
             3.50,3.75,4.0,5.0,6.0]
beta_list = [0.6,0.7,0.8,0.9,1.0,1.15,1.30,1.50,
             1.75,2.0,2.5,3.0,3.5,
             4.0,4.5,5.0,5.5,6.0]
# beta_list = [0.6,0.7,0.8,0.9,1.0,1.15,1.30,1.50,
#              1.75,2.0,3.0,4.0]
# beta_list = [7.0]


# Append sweep results to same list so we can copy paste to plotting script
K_plot = []
K_err_plot = []
for U in U_list:
    for beta in beta_list:
        for mA_sector_wanted in [4]:
        
            incomplete_seeds = [] 
            seeds_list = list(range(1000))
            seeds_measured = []

            # Set desired parameters of calculation
            L = "%d"%(2*mA_sector_wanted)
            N = "%d"%(2*mA_sector_wanted)
            l_max = "%d"%mA_sector_wanted
            beta = "%.6f"%beta
            bin_size = "10000"
            bins_wanted = "1000"
            D = "1"
            U = "%.6f"%(U)
            t = "1.000000"

            # Get path where raw data for the simulation is stored
            path = D+"D_"+L+"_"+N+"_"+l_max+"_"+U+"_"+\
                   t+"_"+beta+"_"+bin_size+"_"+bins_wanted+"/"
            path = "/Users/ecasiano/Desktop/Data/"+path

            # Stores all files in the directory
            filenames_all = os.listdir(path)

            # Convert strings back to numbers
            L = int(L)
            N = int(N)
            l_max = int(l_max)
            beta = float(beta)
            bin_size = int(bin_size)
            bins_wanted = int(bins_wanted)
            D = int(D)
            U = float(U)
            t = float(t)

            # Saves the files relevant to the Renyi Entanglement Entropy calculation
            files_K = []

            # Iterate over all filenames
            for filename in filenames_all:

                # Extract parameter information from file name
                parameters = filename.split("_")

                if parameters[0]=='1D' or parameters[0]=='2D' or parameters[0]=='3D':

                    D_want = int((parameters[0])[0]) # hypercube dimension
                    L_want = int(parameters[1]) # hypercube linear size
                    N_want = int(parameters[2]) # total particles
                    l_want = int(parameters[3]) # subsystem linear size (actually l_max)
                    U_want = float(parameters[4]) # interaction potential
                    t_want = float(parameters[5]) # tunneling parameter
                    beta_want = float(parameters[6]) # imaginary time length (K_B*T)**(-1)
                    bin_size_want = int(parameters[7])
                    bins_want = int(parameters[8]) # number of bins saved in file
                    filetype = (parameters[9]) # identifies the data stored in file
                    seed = int(parameters[10].split(".")[0]) # random seed used

                    if filetype=='K':

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
                                                  'K']

                        if [D,L,N,l_max,U,t,beta,bin_size,bins_wanted,filetype] == parameters_to_evaluate:
                            if os.stat(path+filename).st_size > 0:
                                with open(path+filename) as f:
                                   count = sum(1 for _ in f)
                                if count == bins_wanted: # only consider files that managed to save number of bins wanted
                                    files_K.append(filename)
                                    seeds_measured.append(seed)
                                else:
                                    incomplete_seeds.append(seed)

            # Sort SWAP files in ascending order of random seed
            files_K.sort()

            # Get total number of seeds and total columns per data file
            reference_file = np.loadtxt(path+files_K[0])
            number_of_seeds = len(files_K)
            
            print("\nN: ",N)
            print("L: ",L)
            print("partition size: ", mA_sector_wanted)
            print("U: ",U)
            print("beta: ",beta)

            combined_K_data = np.zeros((number_of_seeds))
            for i,filename in enumerate(files_K):
                data = np.loadtxt(path+filename)[:]
                data_mean = np.mean(data,axis=0)
                combined_K_data[i] = data_mean
        
            # Get mean and std dev,err of S2
            K_mean = np.mean(combined_K_data,axis=0)
            K_stderr = np.std(combined_K_data,axis=0)/np.sqrt(number_of_seeds)

#             # Print out <S2> +/- error
#             for l in range(columns_per_file):
#                 print(f"<S2(l={l})> = {S2_mean[l]:0.8f} +/- {S2_stderr[l]:0.8f}")
            print("<K> = %.4f +/- %.4f"%(K_mean,K_stderr))

            # Save the l=2 results
            K_plot.append(K_mean)
            K_err_plot.append(K_stderr)

print(number_of_seeds)
print("incomplete seeds: ",[int(i) for i in incomplete_seeds])

print("\n\n")
for result in K_plot:
    print(f"{result:0.8f}",end=",")
print("\n\n")
for error in K_err_plot:
    print(f"{error:0.8f}",end=",")

print("\n")
print("<S2>=",K_plot)
print("S2_err=",K_err_plot)
print("beta=",beta_list)
#seeds_measured.sort()
#print([x for x in range(seeds_measured[0],seeds_measured[-1]+1) if x not in seeds_measured])
