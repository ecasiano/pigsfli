//
//  pimc.cpp
//  pimc
//
//  Created by Emanuel Casiano-Diaz on 8/22/20.
//  Copyright © 2020 Emanuel Casiano-Díaz. All rights reserved.
//

#include "pimc.hpp"

// Main
int main(){
    
    // Initialize a Mersenne Twister RNG for each replica
    int seed_A=17,seed_B=42;
//    vector<boost::random::mt19937> rng;
//    boost::random::mt19937 rng_A(seed_A);
//    boost::random::mt19937 rng_B(seed_B);
//    rng.push_back(rng_A);
//    rng.push_back(rng_B);
    
    boost::random::mt19937 rng(seed_A);

    
    // Create a uniform distribution with support: [0.0,1.0)
    boost::random::uniform_real_distribution<double> rnum(0.0,1.0);
    
    // Create integer distribution with support: [0,14]
    boost::random::uniform_int_distribution<> updates(0,14);
    
    // Bose-Hubbard parameters
    int L,D,M,N;
    double t,U,mu;
    string boundary_condition;
    vector<int> initial_fock_state;
    
    // Simulation parameters
    double eta,beta;
    bool canonical;
    unsigned long long int sweeps_pre,sweeps,sweep;
    int label; // random update label;
    
    // Adjacency matrix
    vector<vector<int>> adjacency_matrix;
    int total_nn;
    
    // Declare the data structure
    vector<Kink> kinks_vector;
    vector<vector<Kink>> paths;
    
    // Trackers
//    int num_kinks,head_idx,tail_idx,N_zero,N_beta;
//    double N_tracker;
//    vector<int> last_kinks;
//    unsigned long long int Z_ctr,measurement_attempts;
//    double Z_frac;
    bool not_equilibrated;
    
    // Replicated trackers
    vector<int> num_kinks,head_idx,tail_idx,N_zero,N_beta,bin_ctr;
    vector<double> N_tracker,Z_frac;
    vector<vector<int>> last_kinks;
    vector<unsigned long long int> Z_ctr,measurement_attempts;
    
    // Observables
//    double N_sum,diagonal_energy,kinetic_energy;
//    vector<double> tr_kinetic_energy,tr_diagonal_energy;
    
    // Replicated observables
    vector<double> N_sum,diagonal_energy,kinetic_energy;
    vector<vector<double>> tr_kinetic_energy,tr_diagonal_energy;

    // Measurement settings
    double measurement_center,measurement_plus_minus;
    int measurement_frequency,bin_size;
//    int bin_ctr;
    vector<double> measurement_centers;
    vector<int> fock_state_at_slice;
    
    // Declare data files
//    ofstream kinetic_energy_file,diagonal_energy_file,total_energy_file,
//    tr_kinetic_energy_file,tr_diagonal_energy_file;
    
    vector<ofstream> kinetic_energy_file,diagonal_energy_file,total_energy_file,
    tr_kinetic_energy_file,tr_diagonal_energy_file;
    
    // mu-calibration variables
    double mu_initial,N_hist_sum,P_N_peak,mu_right,mu_left;
    bool N_target_in_bins;
    vector<int> N_data,N_hist,N_bins;
    vector<double> P_N;
    int N_min,N_max,peak_idx,N_idx;
    unsigned long long int  dummy_counter;
    
    // SWAP
    int num_replicas;
    
    // Attempt/Acceptance counters
    unsigned long long int  insert_worm_attempts=0,insert_worm_accepts=0;
    unsigned long long int  delete_worm_attempts=0,delete_worm_accepts=0;

    unsigned long long int  insert_anti_attempts=0,insert_anti_accepts=0;
    unsigned long long int  delete_anti_attempts=0,delete_anti_accepts=0;
    
    unsigned long long int  insertZero_worm_attempts=0,insertZero_worm_accepts=0;
    unsigned long long int  deleteZero_worm_attempts=0,deleteZero_worm_accepts=0;

    unsigned long long int  insertZero_anti_attempts=0,insertZero_anti_accepts=0;
    unsigned long long int  deleteZero_anti_attempts=0,deleteZero_anti_accepts=0;
    
    unsigned long long int  insertBeta_worm_attempts=0,insertBeta_worm_accepts=0;
    unsigned long long int  deleteBeta_worm_attempts=0,deleteBeta_worm_accepts=0;

    unsigned long long int  insertBeta_anti_attempts=0,insertBeta_anti_accepts=0;
    unsigned long long int  deleteBeta_anti_attempts=0,deleteBeta_anti_accepts=0;
    
    unsigned long long int  advance_head_attempts=0,advance_head_accepts=0;
    unsigned long long int  recede_head_attempts=0,recede_head_accepts=0;
    
    unsigned long long int  advance_tail_attempts=0,advance_tail_accepts=0;
    unsigned long long int  recede_tail_attempts=0,recede_tail_accepts=0;
    
    unsigned long long int  ikbh_attempts=0,ikbh_accepts=0;
    unsigned long long int  dkbh_attempts=0,dkbh_accepts=0;
    
    unsigned long long int  ikah_attempts=0,ikah_accepts=0;
    unsigned long long int  dkah_attempts=0,dkah_accepts=0;
    
    unsigned long long int  ikbt_attempts=0,ikbt_accepts=0;
    unsigned long long int  dkbt_attempts=0,dkbt_accepts=0;
    
    unsigned long long int  ikat_attempts=0,ikat_accepts=0;
    unsigned long long int  dkat_attempts=0,dkat_accepts=0;
    
/*------------------------- Initialize variables -----------------------------*/

    // SWAP
    num_replicas=5;
    
    // Bose-Hubbard parameters
    L=2;
    D=1;
    M=pow(L,D);
    N=M;
    t=1.0;
    U=1.0;
    mu=-2.63596;
    boundary_condition="pbc";
    
    // Initialize Fock State
    initial_fock_state = random_boson_config(M,N,rng);
    
    // Simulation parameters
    eta=0.566325;
    beta=1.0;
    canonical=true;
    sweeps=10000000;
    sweeps_pre=1000000;
    sweep=beta*M;
    if (sweep==0){sweep=M;}
    
    // Adjacency matrix
    build_hypercube_adjacency_matrix(L,D,boundary_condition,adjacency_matrix);
    total_nn=0;
    for (int i=0;i<adjacency_matrix[0].size();i++){total_nn+=1;}
            
    // Trackers
//    num_kinks=M;
//    N_tracker=N;
//    head_idx=-1;
//    tail_idx=-1;
//    N_zero=N;
//    N_beta=N;
    
    // Replicated trackers
    for (int i=0;i<num_replicas;i++){
        num_kinks.push_back(M);
        N_tracker.push_back(M);
        head_idx.push_back(-1);
        tail_idx.push_back(-1);
        N_zero.push_back(N);
        N_beta.push_back(N);
        measurement_attempts.push_back(0);
        Z_frac.push_back(0);
        bin_ctr.push_back(0);
    }
    
    // Initialize vector containing indices of last kinks at each site
    for (int r=0;r<num_replicas;r++){
        last_kinks.push_back(vector<int> (M,-1));
        for (int i=0;i<M;i++){last_kinks[r][i]=i;}
    }

//    kinks_vector=create_kinks_vector(initial_fock_state,M);
    paths.push_back(create_kinks_vector(initial_fock_state,M));
    paths.push_back(create_kinks_vector(initial_fock_state,M));
    
    for (int r=0;r<num_replicas;r++){
        paths.push_back(create_kinks_vector(initial_fock_state,M));
//        paths.push_back(create_kinks_vector(initial_fock_state,M));
    }

    // Observables and other measurements
//    N_sum=0.0;
//    diagonal_energy=0.0;
//    kinetic_energy=0.0;
//    Z_ctr=0;
//    measurement_attempts=0;
    
    // Replicated observables and other measurements
    for (int i=0;i<num_replicas;i++){
        N_sum.push_back(0);
        diagonal_energy.push_back(0);
        kinetic_energy.push_back(0);
        Z_ctr.push_back(0);
        measurement_attempts.push_back(0);
//        N_beta.push_back(N);
    }
    
    // Measurement settings
    measurement_center=beta/2.0;
    measurement_plus_minus=0.10*beta;
    measurement_frequency=1;
    bin_size=500;
//    bin_ctr=0;
    measurement_centers=get_measurement_centers(beta);
    for (int i=0;i<M;i++){fock_state_at_slice.push_back(0);}
    
/*------------------- Try drawing a pretty welcome message -------------------*/

    cout << R"(
  _           _   _   _           ______ _____ _____  _____
 | |         | | | | (_)          | ___ \_   _|  __ \/  ___|
 | |     __ _| |_| |_ _  ___ ___  | |_/ / | | | |  \/\ `--.
 | |    / _` | __| __| |/ __/ _ \ |  __/  | | | | __  `--. \
 | |___| (_| | |_| |_| | (_|  __/ | |    _| |_| |_\ \/\__/ /
 \_____/\__,_|\__|\__|_|\___\___| \_|    \___/ \____/\____/
                                                               )";
    
//    cout << R"(
//                                 _
//                                ( `.
//                  _,--------.__  ))\`.
//              _,-"   ,::::.    `(( (  \___
//            ,'      .:::::'      \`-\ |   `-.
//          ,'     ___                         \
//         /     -'   `-.               .    ;; \
//        :::          : \    :         ~)       )-._
//        ;::          :: |   .      .  :       /::._)
//       ( `:          ;  :  /: .    (  :__ .,~/_.-'
//       /__ :        .__/_ (:' ,--.  `./o `.|'
//      ((_\`.    `:.      `-.._    `.__`._o )
// -hrr- `-'  `""`.____,-.___/`_\______ """"`.
//                           `-`       `-. ,\_\
//                                        `-')";
//
    cout << endl << endl;

/*------------------- Pre-equilibration 1: mu calibration --------------------*/

    not_equilibrated=true;
//    mu_initial=mu;
//    dummy_counter=0;
//
//    if (beta>=1.0){sweeps_pre*=(beta*M);}
//    else {sweeps_pre*=M;}
//
//    cout << "Stage (1/4): Determining mu..." << endl << endl;
//
//    // Iterate until particle distribution P(N) is peaked at target N
//    while (true){
//
//        if (!canonical){break;}
//
//        // Restart the data structure and trackers
//        kinks_vector=create_kinks_vector(initial_fock_state,M);
//        for (int i=0;i<M;i++){last_kinks[i] = i;}
//        num_kinks=M;
//        N_tracker=N*1.0;
//        head_idx=-1;
//        tail_idx=-1;
//        N_zero=N;
//        N_beta=N;
//        dummy_counter=0;
//
//        N_data.clear();
//        N_hist.clear();
//        P_N.clear();
//        N_bins.clear();
//        N_hist_sum=0.0;
//        N_min=-1;
//        N_max=-1;
//
//        boost::random::uniform_int_distribution<> updates(0, 14);
//
//        for (unsigned long long int m=0;m<sweeps_pre;m++){
//
//              label = updates(rng);
//
//              if (label==0){     // worm_insert
//                  insert_worm(kinks_vector,num_kinks,head_idx,tail_idx,
//                              M,N,U,mu,t,beta,eta,canonical,N_tracker,
//                              N_zero,N_beta,last_kinks,
//                              dummy_counter,dummy_counter,
//                              dummy_counter,dummy_counter,rng);
//              }
//              else if (label==1){ // worm_delete
//                  delete_worm(kinks_vector,num_kinks,head_idx,tail_idx,
//                              M,N,U,mu,t,beta,eta,canonical,N_tracker,
//                              N_zero,N_beta,last_kinks,
//                              dummy_counter,dummy_counter,
//                              dummy_counter,dummy_counter,rng);
//              }
//              else if (label==2){ // insertZero
//                  insertZero(kinks_vector,num_kinks,head_idx,tail_idx,
//                             M,N,U,mu,t,beta,eta,canonical,N_tracker,
//                             N_zero,N_beta,last_kinks,
//                             dummy_counter,dummy_counter,
//                             dummy_counter,dummy_counter,rng);
//
//              }
//              else if (label==3){ // deleteZero
//                  deleteZero(kinks_vector,num_kinks,head_idx,tail_idx,
//                             M,N,U,mu,t,beta,eta,canonical,N_tracker,
//                             N_zero,N_beta,last_kinks,
//                             dummy_counter,dummy_counter,
//                             dummy_counter,dummy_counter,rng);
//              }
//              else if (label==4){ // insertBeta
//                  insertBeta(kinks_vector,num_kinks,head_idx,tail_idx,
//                             M,N,U,mu,t,beta,eta,canonical,N_tracker,
//                             N_zero,N_beta,last_kinks,
//                             dummy_counter,dummy_counter,
//                             dummy_counter,dummy_counter,rng);
//              }
//              else if (label==5){ // deleteBeta
//                  deleteBeta(kinks_vector,num_kinks,head_idx,tail_idx,
//                             M,N,U,mu,t,beta,eta,canonical,N_tracker,
//                             N_zero,N_beta,last_kinks,
//                             dummy_counter,dummy_counter,
//                             dummy_counter,dummy_counter,rng);
//              }
//              else if (label==6){ // timeshift
//                  timeshift(kinks_vector,num_kinks,head_idx,tail_idx,
//                             M,N,U,mu,t,beta,eta,canonical,N_tracker,
//                             N_zero,N_beta,last_kinks,
//                             dummy_counter,dummy_counter,
//                             dummy_counter,dummy_counter,
//                             dummy_counter,dummy_counter,
//                             dummy_counter,dummy_counter,rng);
//              }
//              else if (label==7){ // insert kink before head
//                  insert_kink_before_head(kinks_vector,num_kinks,head_idx,tail_idx,
//                             M,N,U,mu,t,adjacency_matrix,total_nn,
//                             beta,eta,canonical,N_tracker,
//                             N_zero,N_beta,last_kinks,
//                             dummy_counter,dummy_counter,rng);
//              }
//              else if (label==8){ // delete kink before head
//                  delete_kink_before_head(kinks_vector,num_kinks,head_idx,tail_idx,
//                             M,N,U,mu,t,adjacency_matrix,total_nn,
//                             beta,eta,canonical,N_tracker,
//                             N_zero,N_beta,last_kinks,
//                             dummy_counter,dummy_counter,rng);
//              }
//              else if (label==9){ // insert kink after head
//                  insert_kink_after_head(kinks_vector,num_kinks,head_idx,tail_idx,
//                             M,N,U,mu,t,adjacency_matrix,total_nn,
//                             beta,eta,canonical,N_tracker,
//                             N_zero,N_beta,last_kinks,
//                             dummy_counter,dummy_counter,rng);
//              }
//              else if (label==10){ // delete kink after head
//                  delete_kink_after_head(kinks_vector,num_kinks,head_idx,tail_idx,
//                             M,N,U,mu,t,adjacency_matrix,total_nn,
//                             beta,eta,canonical,N_tracker,
//                             N_zero,N_beta,last_kinks,
//                             dummy_counter,dummy_counter,rng);
//                      }
//              else if (label==11){ // insert kink before tail
//                  insert_kink_before_tail(kinks_vector,num_kinks,head_idx,tail_idx,
//                             M,N,U,mu,t,adjacency_matrix,total_nn,
//                             beta,eta,canonical,N_tracker,
//                             N_zero,N_beta,last_kinks,
//                             dummy_counter,dummy_counter,rng);
//              }
//              else if (label==12){ // delete kink before tail
//                  delete_kink_before_tail(kinks_vector,num_kinks,head_idx,tail_idx,
//                             M,N,U,mu,t,adjacency_matrix,total_nn,
//                             beta,eta,canonical,N_tracker,
//                             N_zero,N_beta,last_kinks,
//                             dummy_counter,dummy_counter,rng);
//              }
//              else if (label==13){ // insert kink after tail
//                   insert_kink_after_tail(kinks_vector,num_kinks,head_idx,tail_idx,
//                              M,N,U,mu,t,adjacency_matrix,total_nn,
//                              beta,eta,canonical,N_tracker,
//                              N_zero,N_beta,last_kinks,
//                              dummy_counter,dummy_counter,rng);
//               }
//               else if (label==14){ // delete kink after tail
//                   delete_kink_after_tail(kinks_vector,num_kinks,head_idx,tail_idx,
//                              M,N,U,mu,t,adjacency_matrix,total_nn,
//                              beta,eta,canonical,N_tracker,
//                              N_zero,N_beta,last_kinks,
//                              dummy_counter,dummy_counter,rng);
//               }
//              else{
//                  // lol
//              }
//
//            // Measure the total number of particles
//            if (head_idx==-1 && tail_idx==-1 &&
//            m%(sweep*measurement_frequency)==0 && m>=0.25*sweeps_pre){
//                N_data.push_back(N_beta);
//                // If we did not collect data, decrease eta and try again.
//                if (m>=0.99*sweeps_pre&&N_data.size()<5){eta*=0.5;break;}
//            }
//        }
//
//        // Find the minimum and maximum number of particles measured
//        N_min=*min_element(N_data.begin(), N_data.end());
//        N_max=*max_element(N_data.begin(), N_data.end());
//
//        // Generate the support of the distribution & initialize the histogram
//        N_target_in_bins=false;
//        for (int i=N_min;i<=N_max;i++){
//            N_bins.push_back(i);
//            N_hist.push_back(0);
//            P_N.push_back(0);
//            if (i==N){N_target_in_bins=true;}
//        }
//
//        // Get the index of the target N in the support of the distribution
//        N_idx = N-N_min;
//
//        // Fill out the histogram
//        for (int i=0;i<N_data.size();i++){
//            N_hist[N_data[i]-N_min]+=1;
//            N_hist_sum+=1.0;
//        }
//
//        // Build the normalized probability distribution P(N) & find its peak
//        peak_idx=0;
//        P_N_peak=P_N[peak_idx];
//        for (int i=0;i<P_N.size();i++){
//            P_N[i]=N_hist[i]/N_hist_sum;
//            if (P_N[i]>P_N_peak){
//                peak_idx=i;
//                P_N_peak=P_N[peak_idx];
//            }
//        }
//
//        // Print out current mu and probability distribution
//        cout << "mu: " << mu << endl;
//        cout << "N     P(N)"<<endl;
//        for (int i=0;i<N_bins.size();i++){
//            cout << setw(6) << left << N_bins[i];
//            for (int j=0;j<=static_cast<int>(100*P_N[i]);j++){
//                cout<<"*";
//            }
//            cout<<endl;
//        }
//        cout << endl;
//
//        if (N_target_in_bins){
//            // Stop the loop if the peak is at P(N)
//            if (peak_idx==N_idx){break;}
//
//            else{
//                // Estimate mu via Eq. 15 in:https://arxiv.org/pdf/1312.6177.pdf
//                if (std::count(N_bins.begin(), N_bins.end(), N-1) &&
//                    std::count(N_bins.begin(), N_bins.end(), N+1)){
//                    mu_right=mu-(1/beta)*log(P_N[N_idx+1]/P_N[N_idx]);
//                    mu_left=mu-(1/beta)*log(P_N[N_idx]/P_N[N_idx-1]);
//                    mu=0.5*(mu_left+mu_right);
//                }
//                else if (std::count(N_bins.begin(), N_bins.end(), N+1)){
//                    mu_right=mu-(1/beta)*log(P_N[N_idx+1]/P_N[N_idx]);
//                    mu=mu_right;
//                }
//                else{
//                    mu_left=mu-(1/beta)*log(P_N[N_idx]/P_N[N_idx-1]);
//                    mu=mu_left;
//                }
//            }
//        }
//        else{ // Target N not in P_N
//            if (N_bins[peak_idx]>N){mu-=0.5;}
//            else{mu+=1;}
//        }
//    }
//
///*------------------ Pre-equilibration 2: eta calibration --------------------*/
//
//    cout << "Stage (2/4): Determining eta..." << endl << endl;
//    cout << setw(16) <<"eta";
//    cout << "| Z-fraction (%)" << endl;
//    cout << "--------------------------------"<<endl;
//
//    // Iterate until particle distribution P(N) is peaked at target N
//      while (true){
//
//          if (!canonical){break;}
//
//          // Restart the data structure and trackers
//          kinks_vector=create_kinks_vector(initial_fock_state,M);
//          for (int i=0;i<M;i++){last_kinks[i] = i;}
//          num_kinks=M;
//          N_tracker=N*1.0;
//          head_idx=-1;
//          tail_idx=-1;
//          N_zero=N;
//          N_beta=N;
//          dummy_counter=0;
//
//          N_data.clear();
//          N_hist.clear();
//          P_N.clear();
//          N_bins.clear();
//          N_hist_sum=0.0;
//          N_min=-1;
//          N_max=-1;
//
//          Z_frac=0.0;
//          measurement_attempts=0;
//
//          boost::random::uniform_int_distribution<> updates(0, 14);
//
//          for (unsigned long long int m=0;m<sweeps_pre;m++){
//
//                label = updates(rng);
//
//                if (label==0){     // worm_insert
//                    insert_worm(kinks_vector,num_kinks,head_idx,tail_idx,
//                                M,N,U,mu,t,beta,eta,canonical,N_tracker,
//                                N_zero,N_beta,last_kinks,
//                                dummy_counter,dummy_counter,
//                                dummy_counter,dummy_counter,rng);
//                }
//                else if (label==1){ // worm_delete
//                    delete_worm(kinks_vector,num_kinks,head_idx,tail_idx,
//                                M,N,U,mu,t,beta,eta,canonical,N_tracker,
//                                N_zero,N_beta,last_kinks,
//                                dummy_counter,dummy_counter,
//                                dummy_counter,dummy_counter,rng);
//                }
//                else if (label==2){ // insertZero
//                    insertZero(kinks_vector,num_kinks,head_idx,tail_idx,
//                               M,N,U,mu,t,beta,eta,canonical,N_tracker,
//                               N_zero,N_beta,last_kinks,
//                               dummy_counter,dummy_counter,
//                               dummy_counter,dummy_counter,rng);
//
//                }
//                else if (label==3){ // deleteZero
//                    deleteZero(kinks_vector,num_kinks,head_idx,tail_idx,
//                               M,N,U,mu,t,beta,eta,canonical,N_tracker,
//                               N_zero,N_beta,last_kinks,
//                               dummy_counter,dummy_counter,
//                               dummy_counter,dummy_counter,rng);
//                }
//                else if (label==4){ // insertBeta
//                    insertBeta(kinks_vector,num_kinks,head_idx,tail_idx,
//                               M,N,U,mu,t,beta,eta,canonical,N_tracker,
//                               N_zero,N_beta,last_kinks,
//                               dummy_counter,dummy_counter,
//                               dummy_counter,dummy_counter,rng);
//                }
//                else if (label==5){ // deleteBeta
//                    deleteBeta(kinks_vector,num_kinks,head_idx,tail_idx,
//                               M,N,U,mu,t,beta,eta,canonical,N_tracker,
//                               N_zero,N_beta,last_kinks,
//                               dummy_counter,dummy_counter,
//                               dummy_counter,dummy_counter,rng);
//                }
//                else if (label==6){ // timeshift
//                    timeshift(kinks_vector,num_kinks,head_idx,tail_idx,
//                               M,N,U,mu,t,beta,eta,canonical,N_tracker,
//                               N_zero,N_beta,last_kinks,
//                               dummy_counter,dummy_counter,
//                               dummy_counter,dummy_counter,
//                               dummy_counter,dummy_counter,
//                               dummy_counter,dummy_counter,rng);
//                }
//                else if (label==7){ // insert kink before head
//                    insert_kink_before_head(kinks_vector,num_kinks,head_idx,tail_idx,
//                               M,N,U,mu,t,adjacency_matrix,total_nn,
//                               beta,eta,canonical,N_tracker,
//                               N_zero,N_beta,last_kinks,
//                               dummy_counter,dummy_counter,rng);
//                }
//                else if (label==8){ // delete kink before head
//                    delete_kink_before_head(kinks_vector,num_kinks,head_idx,tail_idx,
//                               M,N,U,mu,t,adjacency_matrix,total_nn,
//                               beta,eta,canonical,N_tracker,
//                               N_zero,N_beta,last_kinks,
//                               dummy_counter,dummy_counter,rng);
//                }
//                else if (label==9){ // insert kink after head
//                    insert_kink_after_head(kinks_vector,num_kinks,head_idx,tail_idx,
//                               M,N,U,mu,t,adjacency_matrix,total_nn,
//                               beta,eta,canonical,N_tracker,
//                               N_zero,N_beta,last_kinks,
//                               dummy_counter,dummy_counter,rng);
//                }
//                else if (label==10){ // delete kink after head
//                    delete_kink_after_head(kinks_vector,num_kinks,head_idx,tail_idx,
//                               M,N,U,mu,t,adjacency_matrix,total_nn,
//                               beta,eta,canonical,N_tracker,
//                               N_zero,N_beta,last_kinks,
//                               dummy_counter,dummy_counter,rng);
//                        }
//                else if (label==11){ // insert kink before tail
//                    insert_kink_before_tail(kinks_vector,num_kinks,head_idx,tail_idx,
//                               M,N,U,mu,t,adjacency_matrix,total_nn,
//                               beta,eta,canonical,N_tracker,
//                               N_zero,N_beta,last_kinks,
//                               dummy_counter,dummy_counter,rng);
//                }
//                else if (label==12){ // delete kink before tail
//                    delete_kink_before_tail(kinks_vector,num_kinks,head_idx,tail_idx,
//                               M,N,U,mu,t,adjacency_matrix,total_nn,
//                               beta,eta,canonical,N_tracker,
//                               N_zero,N_beta,last_kinks,
//                               dummy_counter,dummy_counter,rng);
//                }
//                else if (label==13){ // insert kink after tail
//                     insert_kink_after_tail(kinks_vector,num_kinks,head_idx,tail_idx,
//                                M,N,U,mu,t,adjacency_matrix,total_nn,
//                                beta,eta,canonical,N_tracker,
//                                N_zero,N_beta,last_kinks,
//                                dummy_counter,dummy_counter,rng);
//                 }
//                 else if (label==14){ // delete kink after tail
//                     delete_kink_after_tail(kinks_vector,num_kinks,head_idx,tail_idx,
//                                M,N,U,mu,t,adjacency_matrix,total_nn,
//                                beta,eta,canonical,N_tracker,
//                                N_zero,N_beta,last_kinks,
//                                dummy_counter,dummy_counter,rng);
//                 }
//                else{
//                    // lol
//                }
//
//              // Measure the total number of particles
//              if (m%(sweep*measurement_frequency)==0 && m>=0.25*sweeps_pre){
//                  measurement_attempts+=1;
//                  if (head_idx==-1 and tail_idx==-1){Z_frac += 1;}
//              }
//          }
//
//          Z_frac/=measurement_attempts;
//          cout << setw(16) << eta;
//          cout << "| " << Z_frac*100.0 << endl;
//
//
//          // Modify eta if necessary
//          if (Z_frac > 0.13 &&
//              Z_frac < 0.17){break;}
//          else{
//              if (Z_frac < 0.17){eta *= 0.5;}
//              else{eta *= 1.5;}
//          }
//      }
    
/*---------------------------- Open files ------------------------------------*/
    
    for (int r=0;r<num_replicas;r++){
        
        ofstream K_out,V_out,tr_K_out,tr_V_out;
        string K_name,V_name,tr_K_name,tr_V_name,rep;
        
        rep=r+65; // rep=65='A'...rep=66='B'...rep=67='C'...
        
        K_name=to_string(L)+"_"+to_string(M)+"_"+
        to_string(U)+"_"+to_string(mu)+"_"+
        to_string(t)+"_"+to_string(beta)+"_"+
        to_string(sweeps)+"_"+"seed_"+to_string(D)+"D_"+
        "can_"+"K_"+"rep"+rep+"_.dat";
        
        V_name=to_string(L)+"_"+to_string(M)+"_"+
        to_string(U)+"_"+to_string(mu)+"_"+
        to_string(t)+"_"+to_string(beta)+"_"+
        to_string(sweeps)+"_"+"seed_"+to_string(D)+"D_"+
        "can_"+"V_"+"rep"+rep+"_.dat";
        
        tr_K_name=to_string(L)+"_"+to_string(M)+"_"+
        to_string(U)+"_"+to_string(mu)+"_"+
        to_string(t)+"_"+to_string(beta)+"_"+
        to_string(sweeps)+"_"+"seed_"+to_string(D)+"D_"+
        "can_"+"tauResolvedK_"+"rep"+rep+".dat";
        
        tr_V_name=to_string(L)+"_"+to_string(M)+"_"+
        to_string(U)+"_"+to_string(mu)+"_"+
        to_string(t)+"_"+to_string(beta)+"_"+
        to_string(sweeps)+"_"+"seed_"+to_string(D)+"D_"+
        "can_"+"tauResolvedV_"+"rep"+rep+".dat";
        
        K_out.open(K_name);
        V_out.open(V_name);
        tr_K_out.open(tr_K_name);
        tr_V_out.open(tr_V_name);
        
        kinetic_energy_file.push_back(std::move(K_out));
        diagonal_energy_file.push_back(std::move(V_out));
        tr_kinetic_energy_file.push_back(std::move(tr_K_out));
        tr_diagonal_energy_file.push_back(std::move(tr_V_out));

//    kinetic_energy_file.open(to_string(L)+"_"+to_string(M)+"_"+
//                             to_string(U)+"_"+to_string(mu)+"_"+
//                             to_string(t)+"_"+to_string(beta)+"_"+
//                             to_string(sweeps)+"_"+"seed_"+to_string(D)+"D_"+
//                             "can"+"_K.dat",fstream::out);
//    if( !kinetic_energy_file ) { // file couldn't be opened
//       cerr << "Error: kinetic energy file could not be opened" << endl;
//       exit(1);
//    }
//
//    diagonal_energy_file.open(to_string(L)+"_"+to_string(M)+"_"+
//                             to_string(U)+"_"+to_string(mu)+"_"+
//                             to_string(t)+"_"+to_string(beta)+"_"+
//                             to_string(sweeps)+"_"+"seed_"+to_string(D)+"D_"+
//                             "can"+"_V.dat",fstream::out);
//    if( !diagonal_energy_file ) { // file couldn't be opened
//       cerr << "Error: diagonal energy file could not be opened" << endl;
//       exit(1);
//    }
//
//    tr_kinetic_energy_file.open(to_string(L)+"_"+to_string(M)+"_"+
//                             to_string(U)+"_"+to_string(mu)+"_"+
//                             to_string(t)+"_"+to_string(beta)+"_"+
//                             to_string(sweeps)+"_"+"seed_"+to_string(D)+"D_"+
//                             "can"+"_tauResolvedK.dat",fstream::out);
//    if( !tr_kinetic_energy_file ) { // file couldn't be opened
//       cerr << "Error: tr kinetic energy file could not be opened" << endl;
//       exit(1);
//    }
//
//    tr_diagonal_energy_file.open(to_string(L)+"_"+to_string(M)+"_"+
//                             to_string(U)+"_"+to_string(mu)+"_"+
//                             to_string(t)+"_"+to_string(beta)+"_"+
//                             to_string(sweeps)+"_"+"seed_"+to_string(D)+"D_"+
//                             "can"+"_tauResolvedV.dat",fstream::out);
//    if( !tr_diagonal_energy_file ) { // file couldn't be opened
//       cerr << "Error: tr diagonal energy file could not be opened" << endl;
//       exit(1);
//    }
//
//    tr_kinetic_energy_file_B.open(to_string(L)+"_"+to_string(M)+"_"+
//                             to_string(U)+"_"+to_string(mu)+"_"+
//                             to_string(t)+"_"+to_string(beta)+"_"+
//                             to_string(sweeps)+"_"+"seed_"+to_string(D)+"D_"+
//                             "can"+"_tauResolvedK_B.dat",fstream::out);
//    if( !tr_kinetic_energy_file_B ) { // file couldn't be opened
//       cerr << "Error: tr kinetic energy file B could not be opened" << endl;
//       exit(1);
//    }
//
//    tr_diagonal_energy_file_B.open(to_string(L)+"_"+to_string(M)+"_"+
//                             to_string(U)+"_"+to_string(mu)+"_"+
//                             to_string(t)+"_"+to_string(beta)+"_"+
//                             to_string(sweeps)+"_"+"seed_"+to_string(D)+"D_"+
//                             "can"+"_tauResolvedV_B.dat",fstream::out);
//    if( !tr_diagonal_energy_file_B ) { // file couldn't be opened
//       cerr << "Error: tr diagonal energy B file could not be opened" << endl;
//       exit(1);
//    }
    }
    
/*---------------------------- Monte Carlo -----------------------------------*/
    
    // Time main function execution
    auto start = high_resolution_clock::now();
    
    // Restart the data structure and trackers
//    kinks_vector=create_kinks_vector(initial_fock_state,M);


//    for (int i=0; i<M; i++){last_kinks[i] = i;}
//    num_kinks=M;
//    N_tracker=N*1.0;
//    head_idx=-1;
//    tail_idx=-1;
//    N_zero=N;
//    N_beta=N;
    
//    Z_frac=0.0;
//    measurement_attempts=0;
    std::fill(Z_frac.begin(),Z_frac.end(),0);
    std::fill(measurement_attempts.begin(),measurement_attempts.end(),0);
    
    if (beta>=1.0){sweeps*=(beta*M);}
    else {sweeps*=M;}
    
//    // Initialize tau resolved estimators
//    for (int i=0;i<measurement_centers.size();i++){
//        tr_kinetic_energy.push_back(0.0);
//        tr_diagonal_energy.push_back(0.0);
//    }
        // Initialize tau resolved estimators
        for (int r=0;r<num_replicas;r++){
            tr_kinetic_energy.push_back(vector<double>
                                        (measurement_centers.size(),0.0));
            tr_diagonal_energy.push_back(vector<double>
                                         (measurement_centers.size(),0.0));
        }
    
    
    cout << endl << "Stage (3/4): Equilibrating..." << endl << endl;
    for (int r=0;r<num_replicas;r++){
    for (unsigned long long int m=0; m < sweeps; m++){
        label = updates(rng);
        
        if (label==0){     // worm_insert
            insert_worm(paths[r],num_kinks[r],head_idx[r],tail_idx[r],
                        M,N,U,mu,t,beta,eta,canonical,N_tracker[r],
                        N_zero[r],N_beta[r],last_kinks[r],
                        insert_worm_attempts,insert_worm_accepts,
                        insert_anti_attempts,insert_anti_accepts,rng);
        }
        else if (label==1){ // worm_delete
            delete_worm(paths[r],num_kinks[r],head_idx[r],tail_idx[r],
                        M,N,U,mu,t,beta,eta,canonical,N_tracker[r],
                        N_zero[r],N_beta[r],last_kinks[r],
                        delete_worm_attempts,delete_worm_accepts,
                        delete_anti_attempts,delete_anti_accepts,rng);
        }
        else if (label==2){ // insertZero
            insertZero(paths[r],num_kinks[r],head_idx[r],tail_idx[r],
                       M,N,U,mu,t,beta,eta,canonical,N_tracker[r],
                       N_zero[r],N_beta[r],last_kinks[r],
                       insertZero_worm_attempts,insertZero_worm_accepts,
                       insertZero_anti_attempts,insertZero_anti_accepts,rng);
            
        }
        else if (label==3){ // deleteZero
            deleteZero(paths[r],num_kinks[r],head_idx[r],tail_idx[r],
                       M,N,U,mu,t,beta,eta,canonical,N_tracker[r],
                       N_zero[r],N_beta[r],last_kinks[r],
                       deleteZero_worm_attempts,deleteZero_worm_accepts,
                       deleteZero_anti_attempts,deleteZero_anti_accepts,rng);
        }
        else if (label==4){ // insertBeta
            insertBeta(paths[r],num_kinks[r],head_idx[r],tail_idx[r],
                       M,N,U,mu,t,beta,eta,canonical,N_tracker[r],
                       N_zero[r],N_beta[r],last_kinks[r],
                       insertBeta_worm_attempts,insertBeta_worm_accepts,
                       insertBeta_anti_attempts,insertBeta_anti_accepts,rng);
        }
        else if (label==5){ // deleteBeta
            deleteBeta(paths[r],num_kinks[r],head_idx[r],tail_idx[r],
                       M,N,U,mu,t,beta,eta,canonical,N_tracker[r],
                       N_zero[r],N_beta[r],last_kinks[r],
                       deleteBeta_worm_attempts,deleteBeta_worm_accepts,
                       deleteBeta_anti_attempts,deleteBeta_anti_accepts,rng);
        }
        else if (label==6){ // timeshift
            timeshift(paths[r],num_kinks[r],head_idx[r],tail_idx[r],
                       M,N,U,mu,t,beta,eta,canonical,N_tracker[r],
                       N_zero[r],N_beta[r],last_kinks[r],
                       advance_head_attempts,advance_head_accepts,
                       recede_head_attempts,recede_head_accepts,
                       advance_tail_attempts,advance_tail_accepts,
                       recede_tail_attempts,recede_tail_accepts,rng);
        }
        else if (label==7){ // insert kink before head
            insert_kink_before_head(paths[r],num_kinks[r],head_idx[r],tail_idx[r],
                       M,N,U,mu,t,adjacency_matrix,total_nn,
                       beta,eta,canonical,N_tracker[r],
                       N_zero[r],N_beta[r],last_kinks[r],
                       ikbh_attempts,ikbh_accepts,rng);
        }
        else if (label==8){ // delete kink before head
            delete_kink_before_head(paths[r],num_kinks[r],head_idx[r],tail_idx[r],
                       M,N,U,mu,t,adjacency_matrix,total_nn,
                       beta,eta,canonical,N_tracker[r],
                       N_zero[r],N_beta[r],last_kinks[r],
                       dkbh_attempts,dkbh_accepts,rng);
        }
        else if (label==9){ // insert kink after head
            insert_kink_after_head(paths[r],num_kinks[r],head_idx[r],tail_idx[r],
                       M,N,U,mu,t,adjacency_matrix,total_nn,
                       beta,eta,canonical,N_tracker[r],
                       N_zero[r],N_beta[r],last_kinks[r],
                       ikah_attempts,ikah_accepts,rng);
        }
        else if (label==10){ // delete kink after head
            delete_kink_after_head(paths[r],num_kinks[r],head_idx[r],tail_idx[r],
                       M,N,U,mu,t,adjacency_matrix,total_nn,
                       beta,eta,canonical,N_tracker[r],
                       N_zero[r],N_beta[r],last_kinks[r],
                       dkah_attempts,dkah_accepts,rng);
                }
        else if (label==11){ // insert kink before tail
            insert_kink_before_tail(paths[r],num_kinks[r],head_idx[r],tail_idx[r],
                       M,N,U,mu,t,adjacency_matrix,total_nn,
                       beta,eta,canonical,N_tracker[r],
                       N_zero[r],N_beta[r],last_kinks[r],
                       ikbt_attempts,ikbt_accepts,rng);
        }
        else if (label==12){ // delete kink before tail
            delete_kink_before_tail(paths[r],num_kinks[r],head_idx[r],tail_idx[r],
                       M,N,U,mu,t,adjacency_matrix,total_nn,
                       beta,eta,canonical,N_tracker[r],
                       N_zero[r],N_beta[r],last_kinks[r],
                       dkbt_attempts,dkbt_accepts,rng);
        }
        else if (label==13){ // insert kink after tail
             insert_kink_after_tail(paths[r],num_kinks[r],head_idx[r],tail_idx[r],
                        M,N,U,mu,t,adjacency_matrix,total_nn,
                        beta,eta,canonical,N_tracker[r],
                        N_zero[r],N_beta[r],last_kinks[r],
                        ikat_attempts,ikat_accepts,rng);
         }
         else if (label==14){ // delete kink after tail
             delete_kink_after_tail(paths[r],num_kinks[r],head_idx[r],tail_idx[r],
                        M,N,U,mu,t,adjacency_matrix,total_nn,
                        beta,eta,canonical,N_tracker[r],
                        N_zero[r],N_beta[r],last_kinks[r],
                        dkat_attempts,dkat_accepts,rng);
         }
        else{
            // lol
        }
        
        
/*------------------------- Unit Tests (kind of) -----------------------------*/
//
//        // Unit test #1: No last kink indices should be repeated
//        if (last_kinks[0]==last_kinks[1]
//           || last_kinks[0]==last_kinks[2]
//           || last_kinks[0]==last_kinks[3]
//           || last_kinks[1]==last_kinks[2]
//           || last_kinks[1]==last_kinks[3]
//           || last_kinks[2]==last_kinks[3]){
//            cout << "ERROR: Every site should have different last indices";
//            cout << "label " << label << " m " << m << endl;
//            // Print out the indices of each sites last kink
//            cout << "Last kiNk indices before update: ";
//            for (int i=0; i<M ; i++){
//                cout << last_kinks[i] << " ";
//            }
//            cout << endl;
//            break;
//        }
//
//
//        // Unit test #2: src and dest of worm tail kink should be the same
//        if (tail_idx!=-1){
//            if (kinks_vector[tail_idx].src!=kinks_vector[tail_idx].dest){
//                cout << "ERROR: src,dest of worm tail not the same ";
//                cout << "label: " << label << " m " << m <<
//                 " num_kinks: " << num_kinks << " head_idx: " << head_idx <<
//                 " tail_idx: " << tail_idx << endl;
//                // Print out the indices of each sites last kink
//                cout << "Structure after worm end idx error: " << endl;
//                for (int i=0; i<num_kinks+5 ; i++){
//                    cout << i << " " << kinks_vector[i] << endl;
//                }
//                cout << endl;
//                break;
//            }
//        }
//
//        // Unit test 3: src and dest of worm head kink should be the same
//        if (head_idx!=-1){
//            if (kinks_vector[head_idx].src!=kinks_vector[head_idx].dest){
//                cout << "ERROR: src,dest of worm head not the same ";
//                cout << "label: " << label << " m " << m <<
//                 " num_kinks: " << num_kinks << " head_idx: " << head_idx <<
//                 " tail_idx: " << tail_idx << endl;
//                // Print out the indices of each sites last kink
//                cout << "Structure after worm end idx error: " << endl;
//                for (int i=0; i<num_kinks+5 ; i++){
//                    cout << i << " " << kinks_vector[i] << endl;
//                }
//                cout << endl;
//                break;
//            }
//        }
//
//        // Unit test 4: Last indices should have .next equal to -1
//        if (kinks_vector[last_kinks[0]].next!=-1
//            || kinks_vector[last_kinks[1]].next!=-1
//            || kinks_vector[last_kinks[2]].next!=-1
//            || kinks_vector[last_kinks[3]].next!=-1){
//            cout << "ERROR: Last indices should have .next equal to -1 ";
//            cout << "label: " << label << " m " << m <<
//             " num_kinks: " << num_kinks << " head_idx: " << head_idx <<
//             " tail_idx: " << tail_idx << endl;
//            // Print out the indices of each sites last kink
//            cout << "Last indices: " << endl;
//            for (int i=0; i<M; i++){
//                cout << last_kinks[i] << " ";
//            }
//            cout << endl;
//            cout << ".next of each of the last kinks:"<<endl;
//            for (int i=0; i<M; i++){
//                cout << kinks_vector[last_kinks[i]].next << " ";
//            }
//            cout << endl;
//            cout << "Structure after worm end idx error: " << endl;
//            for (int i=0; i<num_kinks+5 ; i++){
//                cout << i << " " << kinks_vector[i] << endl;
//            }
//            cout << endl;
//            break;
//        }
//
//        // Unit test 5: Conservation of N_zero
//        if (head_idx==-1 && tail_idx==-1 && canonical){
//            if (N_tracker < N-1 || N_tracker > N+1){
//                cout << "ERROR: Total particle number N not conserved" << endl;
//                cout << "label: " << label << " m " << m <<
//                 " num_kinks: " << num_kinks << " head_idx: " << head_idx <<
//                 " tail_idx: " << tail_idx << endl;
//                // Print out the indices of each sites last kink
//                cout << "Structure after worm end idx error: " << endl;
//                for (int i=0; i<num_kinks+5 ; i++){
//                    cout << i << " " << kinks_vector[i] << endl;
//                }
//                break;
//            }
//        }
//
//        // Unit Test 6: Sum the particles at the beggining when no worm ends
//        if (head_idx==-1 && tail_idx==-1 && canonical){
//            int N_total_sum_zero=0;
//            for (int i=0; i<M; i++){
//                N_total_sum_zero+=kinks_vector[i].n;
//            }
//            if (N_total_sum_zero < N-1 || N_total_sum_zero > N+1){
//                cout << "ERROR: Total particle number N at zero not conserved" << endl;
//                cout << "label: " << label << " m " << m <<
//                 " num_kinks: " << num_kinks << " head_idx: " << head_idx <<
//                 " tail_idx: " << tail_idx << endl;
//                cout << "N_zero (unit test 6): " << N_total_sum_zero << endl;
//                // Print out the indices of each sites last kink
//                cout << "Structure after worm end idx error: " << endl;
//                for (int i=0; i<num_kinks+5 ; i++){
//                    cout << i << " " << kinks_vector[i] << endl;
//                }
//                break;
//            }
//        }
//
//        // Unit Test 7: Sum the particles at the end when no worm ends
//        if (head_idx==-1 && tail_idx==-1 && canonical){
//            int N_total_sum_beta=0;
//            for (int i=0; i<M; i++){
//                N_total_sum_beta+=kinks_vector[last_kinks[i]].n;
//            }
//            if (N_total_sum_beta < N-1 || N_total_sum_beta > N+1){
//                cout << "ERROR: Total particle number at beta not conserved" << endl;
//                cout << "label: " << label << " m " << m <<
//                 " num_kinks: " << num_kinks << " head_idx: " << head_idx <<
//                 " tail_idx: " << tail_idx << " N_tracker: " << N_tracker <<endl;
//                cout << "N_beta (unit test 7): " << N_total_sum_beta << endl;
//                // Print out the indices of each sites last kink
//                cout << "Structure after worm end idx error: " << endl;
//                for (int i=0; i<num_kinks+5 ; i++){
//                    cout << i << " " << kinks_vector[i] << endl;
//                }
//                break;
//            }
//        }

//        // Unit test 8: Check that the prev,next attributes are different to kink index
//        for (int i=0; i<num_kinks; i++){
//            if (i==kinks_vector[i].prev
//                || i==kinks_vector[i].next){
//                cout << "ERROR: the kink with index " <<
//                i << " has the same prev or next" << endl;
//                break;}
//        }
        
/*----------------------------- Measurements ---------------------------------*/

        if (m%(sweep*measurement_frequency)==0 && m>=sweeps*0.25){
            
            if (not_equilibrated){
                not_equilibrated=false;
                cout << "Stage (4/4): Main Monte Carlo loop..." << endl;
            }
            
            measurement_attempts[r]+=1;
            if (head_idx[r]==-1 and tail_idx[r]==-1){
                N_sum[r] += N_tracker[r];
                Z_ctr[r] += 1;
                           
                if (N_beta[r]==N){ // canonical measurement
                    
                // Get fock state at desired measurement center
                get_fock_state(measurement_center,M,fock_state_at_slice,
                               paths[r]);
                    
                // Measure and accumulate <K>
                kinetic_energy[r]+=pimc_kinetic_energy(paths[r],num_kinks[r],
                            measurement_center,measurement_plus_minus,M,t,beta);
                    
                // Measure and accumulate <V>
                diagonal_energy[r]+=pimc_diagonal_energy(fock_state_at_slice,
                                                      M,canonical,U,mu);
                    
                tau_resolved_kinetic_energy(paths[r],num_kinks[r],M,t,beta,
                                            measurement_centers,
                                            tr_kinetic_energy[r]);
                    
                tau_resolved_diagonal_energy(paths[r],num_kinks[r],
                                            M,canonical,U,mu,beta,
                                            measurement_centers,
                                            tr_diagonal_energy[r]);
                    
                bin_ctr[r]+=1;
                // Take binned averages and write to disk
                if (bin_ctr[r]==bin_size){
                    kinetic_energy_file[r]<<fixed<<setprecision(17)<<
                    kinetic_energy[r]/bin_size<<endl;
                    diagonal_energy_file[r]<<fixed<<setprecision(17)<<
                    diagonal_energy[r]/bin_size<<endl;
                    
                    // Save tau resolved estimators
                    for (int i=0; i<measurement_centers.size(); i++){
                        tr_kinetic_energy_file[r]<<fixed<<setprecision(17)<<
                        tr_kinetic_energy[r][i]/bin_size << " ";
                        
                        tr_diagonal_energy_file[r]<<fixed<<setprecision(17)<<
                        tr_diagonal_energy[r][i]/bin_size << " ";
                    }
                    tr_kinetic_energy_file[r]<<endl;
                    tr_diagonal_energy_file[r]<<endl;
                    
                    bin_ctr[r]=0;
                    kinetic_energy[r]=0.0;
                    diagonal_energy[r]=0.0;
                    std::fill(tr_kinetic_energy[r].begin(),
                              tr_kinetic_energy[r].end(),0);
                    std::fill(tr_diagonal_energy[r].begin(),
                              tr_diagonal_energy[r].end(),0);
                    }
                }
            }
        }
    } // end of replica loop
    } // end of sweep loop
    
    // Close data files
    for (int r=0;r<num_replicas;r++){
        kinetic_energy_file[r].close();
        diagonal_energy_file[r].close();
        tr_kinetic_energy_file[r].close();
        tr_diagonal_energy_file[r].close();
    }
/*--------------------------------- FIN --------------------------------------*/

    cout << endl << "-------- Detailed Balance --------" << endl;

    cout<< endl << "Insert Worm: "<<insert_worm_accepts<<"/"<<
                                    insert_worm_attempts<<endl;
    cout<<         "Delete Worm: "<<delete_worm_accepts<<"/"<<
                                    delete_worm_attempts<<endl;
    
    cout<< endl <<"Insert Anti: "<<insert_anti_accepts<<"/"<<
                           insert_anti_attempts<<endl;
    cout<<"Delete Anti: "<<delete_anti_accepts<<"/"<<
                           delete_anti_attempts<<endl;
    
    cout<< endl <<"InsertZero Worm: "<<insertZero_worm_accepts<<"/"<<
                               insertZero_worm_attempts<<endl;
    cout<<"DeleteZero Worm: "<<deleteZero_worm_accepts<<"/"<<
                               deleteZero_worm_attempts<<endl;
    
    cout<< endl <<"InsertZero Anti: "<<insertZero_anti_accepts<<"/"<<
                               insertZero_anti_attempts<<endl;
    cout<<"DeleteZero Anti: "<<deleteZero_anti_accepts<<"/"<<
                               deleteZero_anti_attempts<<endl;
    
    cout<< endl <<"InsertBeta Worm: "<<insertBeta_worm_accepts<<"/"<<
                               insertBeta_worm_attempts<<endl;
    cout<<"DeleteBeta Worm: "<<deleteBeta_worm_accepts<<"/"<<
                               deleteBeta_worm_attempts<<endl;
    
    cout<< endl <<"InsertBeta Anti: "<<insertBeta_anti_accepts<<"/"<<
                               insertBeta_anti_attempts<<endl;
    cout<<"DeleteBeta Anti: "<<deleteBeta_anti_accepts<<"/"<<
                               deleteBeta_anti_attempts<<endl;
    
    cout<< endl <<"Advance Head: "<<advance_head_accepts<<"/"<<
                               advance_head_attempts<<endl;
    cout<<"Recede  Head: "<< recede_head_accepts<<"/"<<
                               recede_head_attempts<<endl;
    
    cout<< endl <<"Advance Tail: "<<advance_tail_accepts<<"/"<<
                               advance_tail_attempts<<endl;
    cout<<"Recede  Tail: "<<recede_tail_accepts<<"/"<<
                               recede_tail_attempts<<endl;
    
    cout<< endl <<"IKBH: "<<ikbh_accepts<<"/"<<
                               ikbh_attempts<<endl;
    cout <<"DKBH: "<<dkbh_accepts<<"/"<<
                               dkbh_attempts<<endl;
    
    cout<< endl <<"IKAH: "<<ikah_accepts<<"/"<<
                               ikah_attempts<<endl;
    cout <<"DKAH: "<<dkah_accepts<<"/"<<
                               dkah_attempts<<endl;
    
    cout<< endl <<"IKBT: "<<ikbt_accepts<<"/"<<
                               ikbt_attempts<<endl;
    cout <<"DKBT: "<<dkbt_accepts<<"/"<<
                               dkbt_attempts<<endl;
    
    cout<< endl <<"IKAT: "<<ikat_accepts<<"/"<<
                               ikat_attempts<<endl;
    cout <<"DKAT: "<<dkat_accepts<<"/"<<
                               dkat_attempts<<endl;
    
    auto end = high_resolution_clock::now();

    auto elapsed_time = duration_cast<nanoseconds>(end - start);
    double duration = elapsed_time.count() * 1e-9;
    
    cout << endl << "beta: " << beta << endl;
    cout << endl << "sweeps: " << sweeps/(beta*M) << endl;
    cout << "Z_ctr: " << Z_ctr[0] << endl;
    cout<<"Z_frac: "<<Z_ctr[0]*100.0/measurement_attempts[0]<<"% ("<<Z_ctr[0]
    <<"/"<< measurement_attempts[0]<<")"<<endl;
    
    cout << endl << "<N>: " << (N_sum[0])/Z_ctr[0] << endl;

    cout << endl << "Elapsed time: " << duration << " seconds" << endl;

    return 0;
    
}
