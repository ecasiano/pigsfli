# Searches directory for all Bose-Hubbard ground state energy (egs) file
# at fixed interaction strength but varying variational parameter

import os
import numpy as np

# Save all the file names in the path as strings to a list
path = './'
files_all = os.listdir(path)

# Only keep the files that contain the string 'egs_' in their name
keep_string = '-sector'  # The string that kept files will contain
n_sector_files = []
for file in files_all:
    if keep_string in file: 
        print("KAKAROT")
        n_sector_files.append(file)
        
# Sort egs files in ascending order n-sector
n_sector_files.sort()

# Initialize arrays to save v and average egs
N_elements = np.shape(files_egs)[0]
v_list = np.ones(N_elements)
egs_list = np.ones(N_elements)

print(files_egs)
# # Set from where to start the data depending on equilibration time
# equil_time = 0.2 #Set the percentage of unequilibrated data to throw away
# data = np.loadtxt('Data/egs_4_4_9.9763_1.2500_150000.dat')
# begin_data = int(np.shape(data)[0]*equil_time) # Index of starting element after equilibration

# for i,file_name in enumerate(files_egs):
    
#     file_name = path+file_name # Pre-pend the directory where the data file is located
#     # Load equilibrated data file for each variational parameter v
#     data = np.loadtxt(file_name)[begin_data:]

#     # Extract BH parameters from file name
#     L,N,U,v,M = file_name.split("_")[1:]
#     M = M.split(".")[0]
#     L,N,U,v,M = int(L),int(N),float(U),float(v),int(M) #Promote from str to int OR float
    
#     # Save v and average egs to array
#     v_list[i] = v
#     egs_list[i] = np.mean(data)

#     print("%.4f %%"%((i+1)/N_elements * 100))
    
# # Reshape arrays as columns
# v_list = np.reshape(v_list,(N_elements,1))
# egs_list = np.reshape(egs_list,(N_elements,1))
# cols = np.hstack((v_list,egs_list))
# with open("vsweep_%i_%i_%.4f_%.4f_%i.dat"%(L,N,U,v,M),"w+") as data:
#         np.savetxt(data,cols,delimiter=" ",fmt="%.16f",header="MC_step Egs // BH Parameters: L=%d,N=%d,U=%.14f,v=%.14f,MC_steps=%i"%(L,N,U,v,M))