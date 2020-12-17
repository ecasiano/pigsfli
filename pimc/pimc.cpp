//
//  pimc.cpp
//  pimc
//
//  Created by Emanuel Casiano-Diaz on 8/22/20.
//  Copyright © 2020 Emanuel Casiano-Díaz. All rights reserved.
//

#include "pimc.hpp"
#include <cxxopts.hpp>

// Main
int main(int argc, char** argv){


  /*------------------------- Command line arguments -----------------------------*/

    cxxopts::Options options("test", "A brief description");

    options.add_options()
        ("h,help", "Print usage")
        ("D","Dimension of hypercubic lattice",cxxopts::value<int>())
        ("L","Linear size of hypercube",cxxopts::value<int>())
        ("N","Total number of particles",cxxopts::value<int>())
        ("U","Interaction potential",cxxopts::value<double>())
        ("l","Linear size of hypercubic subregion",cxxopts::value<int>())
        ("sweeps","Number of Monte Carlo sweeps",cxxopts::value<unsigned long long int>())
        ("beta","Set length of imaginary time",cxxopts::value<double>())

        ("mu","Chemical potential",cxxopts::value<double>()->default_value("-3.0"))
        ("t","Tunneling parameter",cxxopts::value<double>()->default_value("1.0"))
        ("canonical", "set to false for grand canonical simulation",
            cxxopts::value<bool>()->default_value("true"))
        ("seed","Random seed value",cxxopts::value<int>()->default_value("0"))
        ("sweeps-pre","Number pre-equilibration sweeps",
            cxxopts::value<unsigned long long int>()->default_value("1000000"))
        ("bin-size","Number of measurements per bin",cxxopts::value<int>())
//        ("conventionals", "set to take conventional measurements (E,N,...)")
    ;

    auto result = options.parse(argc, argv);

  /*-----------------------------------------------------------------------------*/

    // Initialize a Mersenne Twister RNG for each replica
    int seed = result["seed"].as<int>();
    boost::random::mt19937 rng(seed);
    
    // Create a uniform distribution with support: [0.0,1.0)
    boost::random::uniform_real_distribution<double> rnum(0.0,1.0);

    // Create integer distribution with support: [0,16]
    boost::random::uniform_int_distribution<> updates(0,14);
    
    // Create integer distribution with support: [0,2]
    boost::random::uniform_int_distribution<> swap_updates(0,2);
    
    // Bose-Hubbard parameters
    int L,D,M,N;
    double t,U,mu;
    string boundary_condition;
    vector<int> initial_fock_state;
    
    // Simulation parameterss
    double eta,beta;
    bool canonical;
    unsigned long long int sweeps_pre,sweeps,sweep;
    int label; // random update label;
    
    // Adjacency matrix
    vector<vector<int> > adjacency_matrix;
    int total_nn;
    
    // Declare the data structure
    vector<vector<Kink> > paths;
        
    // Replicated trackers
    vector<int> num_kinks,head_idx,tail_idx,N_zero,N_beta,bin_ctr;
    vector<double> N_tracker;
    vector<vector<int> > last_kinks;
    vector<unsigned long long int> Z_ctr,measurement_attempts;
    
    // Replicated observables
    vector<double> N_sum,diagonal_energy,kinetic_energy;
    vector<vector<double> > tr_kinetic_energy,tr_diagonal_energy;
    
    // <SWAP> estimator settings and trackers
    int l_A; // subregion linear size
    int m_A; // subregion total size
    vector<int> sub_sites, swapped_sites;
    vector<int> swap_kinks;
    int num_swaps;
    vector<int> SWAP_histogram;
    vector<vector<int> > n_resolved_SWAP_histograms;
    
    // Measurement settings
    double measurement_center,measurement_plus_minus;
    int measurement_frequency,bin_size,writing_ctr;
    vector<double> measurement_centers;
    vector<int> fock_state_at_slice,fock_state_at_half_plus;
 
    // Declare data files
    vector<ofstream> kinetic_energy_file,diagonal_energy_file,total_energy_file,
    tr_kinetic_energy_file,tr_diagonal_energy_file;
    ofstream SWAP_histogram_file;
    vector<ofstream> n_resolved_SWAP_histogram_files;
    ofstream n_resolved_SWAP_histogram_file;

    // mu-calibration variables
    bool not_equilibrated;
    double mu_initial,N_hist_sum,P_N_peak,mu_right,mu_left,N_flats_mean;
    double Z_frac,N_mean_pre; // used only in eta-equilibration
    bool N_target_in_bins;
    vector<int> N_data,N_hist,N_bins;
    vector<double> P_N;
    int N_min,N_max,peak_idx,N_idx;
    unsigned long long int  dummy_counter,N_flats_samples;
    
    // SWAP
    int num_replicas;
    
    // Related to accessible entanglement
    int n_A; // total particles in subsystem
    
    
    // Attempt/Acceptance counters
    unsigned long long int insert_worm_attempts=0,insert_worm_accepts=0;
    unsigned long long int delete_worm_attempts=0,delete_worm_accepts=0;

    unsigned long long int insert_anti_attempts=0,insert_anti_accepts=0;
    unsigned long long int delete_anti_attempts=0,delete_anti_accepts=0;
    
    unsigned long long int insertZero_worm_attempts=0,insertZero_worm_accepts=0;
    unsigned long long int deleteZero_worm_attempts=0,deleteZero_worm_accepts=0;

    unsigned long long int insertZero_anti_attempts=0,insertZero_anti_accepts=0;
    unsigned long long int deleteZero_anti_attempts=0,deleteZero_anti_accepts=0;
    
    unsigned long long int insertBeta_worm_attempts=0,insertBeta_worm_accepts=0;
    unsigned long long int deleteBeta_worm_attempts=0,deleteBeta_worm_accepts=0;

    unsigned long long int insertBeta_anti_attempts=0,insertBeta_anti_accepts=0;
    unsigned long long int deleteBeta_anti_attempts=0,deleteBeta_anti_accepts=0;
    
    unsigned long long int advance_head_attempts=0,advance_head_accepts=0;
    unsigned long long int recede_head_attempts=0,recede_head_accepts=0;
    
    unsigned long long int advance_tail_attempts=0,advance_tail_accepts=0;
    unsigned long long int recede_tail_attempts=0,recede_tail_accepts=0;
    
    unsigned long long int ikbh_attempts=0,ikbh_accepts=0;
    unsigned long long int dkbh_attempts=0,dkbh_accepts=0;
    
    unsigned long long int ikah_attempts=0,ikah_accepts=0;
    unsigned long long int dkah_attempts=0,dkah_accepts=0;
    
    unsigned long long int ikbt_attempts=0,ikbt_accepts=0;
    unsigned long long int dkbt_attempts=0,dkbt_accepts=0;
    
    unsigned long long int ikat_attempts=0,ikat_accepts=0;
    unsigned long long int dkat_attempts=0,dkat_accepts=0;
    
    unsigned long long int insert_swap_kink_attempts=0,insert_swap_kink_accepts=0;
    unsigned long long int delete_swap_kink_attempts=0,delete_swap_kink_accepts=0;
    
    unsigned long long int swap_advance_head_attempts=0,swap_advance_head_accepts=0;
    unsigned long long int swap_recede_head_attempts=0,swap_recede_head_accepts=0;
        
/*------------------------- Initialize variables -----------------------------*/

    // SWAP
    num_replicas=2;
    
    // Bose-Hubbard parameters
    D=result["D"].as<int>();
    L=result["L"].as<int>();
    M=pow(L,D);
    N=result["N"].as<int>();
    t=result["t"].as<double>();
    U=result["U"].as<double>();
    mu=result["mu"].as<double>();
    boundary_condition="pbc";
   
    // Subsystem settings
    l_A = result["l"].as<int>(); // subsystem linear size
    if (l_A>=L){cout << "ERROR: Please choose l < L" << endl;exit(1);}
    m_A = pow(l_A,D);
    create_sub_sites(sub_sites,l_A,L,D,M);
    num_swaps=0;
    
    // Initialize Fock State
    initial_fock_state = random_boson_config(M,N,rng);
    
    // Simulation parameters
    eta=1/sqrt(M);
    beta=result["beta"].as<double>();
    canonical=result["canonical"].as<bool>();
    sweeps=result["sweeps"].as<unsigned long long int>();
    sweeps_pre=result["sweeps-pre"].as<unsigned long long int>();
    sweep=beta*M;
    if (sweep==0){sweep=M;} // in case beta<1.0
    
    // Adjacency matrix
    build_hypercube_adjacency_matrix(L,D,boundary_condition,adjacency_matrix);
    total_nn=0;
    for (int i=0;i<adjacency_matrix[0].size();i++){total_nn+=1;}
    
    // Replicated trackers
    for (int r=0;r<num_replicas;r++){
        num_kinks.push_back(M);
        N_tracker.push_back(N);
        head_idx.push_back(-1);
        tail_idx.push_back(-1);
        N_zero.push_back(N);
        N_beta.push_back(N);
        measurement_attempts.push_back(0);
        bin_ctr.push_back(0);
        
        // Initialize vector containing indices of last kinks at each site
        last_kinks.push_back(vector<int> (M,-1));
        for (int i=0;i<M;i++){last_kinks[r][i]=i;}
        
        // Worldlines data structure
        paths.push_back(create_paths(initial_fock_state,M,r));
        
        // Observables and other measurements
        N_sum.push_back(0);
        diagonal_energy.push_back(0);
        kinetic_energy.push_back(0);
        Z_ctr.push_back(0);
        measurement_attempts.push_back(0);
    }
    
    // just initializing
    for (int i=0; i<l_A; i++){
        swap_kinks.push_back(0);
    }

    // Measurement settings
    measurement_center=beta/2.0;
    measurement_plus_minus=0.10*beta;
    measurement_frequency=1;
    bin_size=result["bin-size"].as<int>();
    measurement_centers=get_measurement_centers(beta);
    for (int i=0;i<M;i++){
        fock_state_at_slice.push_back(0);
        fock_state_at_half_plus.push_back(0);
    }
    writing_ctr = 0;
    
    N_flats_mean=0.0;
    N_flats_samples=0;
        
/*------------------- Try drawing a pretty welcome message -------------------*/

for (int i=0; i<sub_sites.size(); i++){
    cout << sub_sites[i] << " ";
}
cout << endl;

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

    bool at_least_one_iteration = false;
    
    not_equilibrated=true;
    mu_initial=mu;
    dummy_counter=0;

    if (beta>=1.0){sweeps_pre*=(beta*M);}
    else {sweeps_pre*=M;}

    cout << "Stage (1/3): Determining mu and eta..." << endl << endl;

    // Iterate until particle distribution P(N) is peaked at target N
    while (true){
        
        if (!canonical){break;}

        // Restart data structure and trackers
        num_kinks.clear();
        N_tracker.clear();
        head_idx.clear();
        tail_idx.clear();
        N_zero.clear();
        N_beta.clear();
        last_kinks.clear();
        paths.clear();
        for (int r=0;r<num_replicas;r++){
            num_kinks.push_back(M);
            N_tracker.push_back(N);
            head_idx.push_back(-1);
            tail_idx.push_back(-1);
            N_zero.push_back(N);
            N_beta.push_back(N);
            
            last_kinks.push_back(vector<int> (M,-1));
            for (int i=0;i<M;i++){last_kinks[r][i]=i;}
            
            paths.push_back(create_paths(initial_fock_state,M,r));
        }
        
        N_data.clear();
        N_hist.clear();
        P_N.clear();
        N_bins.clear();
        N_hist_sum=0.0;
        N_min=-1;
        N_max=-1;
        N_mean_pre=0.0;
        
        N_flats_mean=0.0;
        N_flats_samples=0;
        
        Z_frac=0.0; // think about making this a vector too
        std::fill(measurement_attempts.begin(),measurement_attempts.end(),0);
        
        boost::random::uniform_int_distribution<> updates(0, 14);

        for (unsigned long long int m=0;m<sweeps_pre;m++){

              label = updates(rng);

              if (label==0){     // worm_insert
                  insert_worm(paths[0],num_kinks[0],head_idx[0],tail_idx[0],
                              M,N,U,mu,t,beta,eta,canonical,N_tracker[0],
                              N_zero[0],N_beta[0],last_kinks[0],
                              dummy_counter,dummy_counter,
                              dummy_counter,dummy_counter,rng);
              }
              else if (label==1){ // worm_delete
                  delete_worm(paths[0],num_kinks[0],head_idx[0],tail_idx[0],
                              M,N,U,mu,t,beta,eta,canonical,N_tracker[0],
                              N_zero[0],N_beta[0],last_kinks[0],
                              dummy_counter,dummy_counter,
                              dummy_counter,dummy_counter,rng);
              }
              else if (label==2){ // insertZero
                  insertZero(paths[0],num_kinks[0],head_idx[0],tail_idx[0],
                             M,N,U,mu,t,beta,eta,canonical,N_tracker[0],
                             N_zero[0],N_beta[0],last_kinks[0],
                             dummy_counter,dummy_counter,
                             dummy_counter,dummy_counter,rng);

              }
              else if (label==3){ // deleteZero
                  deleteZero(paths[0],num_kinks[0],head_idx[0],tail_idx[0],
                             M,N,U,mu,t,beta,eta,canonical,N_tracker[0],
                             N_zero[0],N_beta[0],last_kinks[0],
                             dummy_counter,dummy_counter,
                             dummy_counter,dummy_counter,rng);
              }
              else if (label==4){ // insertBeta
                  insertBeta(paths[0],num_kinks[0],head_idx[0],tail_idx[0],
                             M,N,U,mu,t,beta,eta,canonical,N_tracker[0],
                             N_zero[0],N_beta[0],last_kinks[0],
                             dummy_counter,dummy_counter,
                             dummy_counter,dummy_counter,rng);
              }
              else if (label==5){ // deleteBeta
                  deleteBeta(paths[0],num_kinks[0],head_idx[0],tail_idx[0],
                             M,N,U,mu,t,beta,eta,canonical,N_tracker[0],
                             N_zero[0],N_beta[0],last_kinks[0],
                             dummy_counter,dummy_counter,
                             dummy_counter,dummy_counter,rng);
              }
              else if (label==6){ // timeshift
                  timeshift(paths[0],num_kinks[0],head_idx[0],tail_idx[0],
                             M,N,U,mu,t,beta,eta,canonical,N_tracker[0],
                             N_zero[0],N_beta[0],last_kinks[0],
                             dummy_counter,dummy_counter,
                             dummy_counter,dummy_counter,
                             dummy_counter,dummy_counter,
                             dummy_counter,dummy_counter,rng);
              }
              else if (label==7){ // insert kink before head
                  insert_kink_before_head(paths[0],num_kinks[0],
                             head_idx[0],tail_idx[0],
                             M,N,U,mu,t,adjacency_matrix,total_nn,
                             beta,eta,canonical,N_tracker[0],
                             N_zero[0],N_beta[0],last_kinks[0],
                             dummy_counter,dummy_counter,rng);
              }
              else if (label==8){ // delete kink before head
                  delete_kink_before_head(paths[0],num_kinks[0],
                             head_idx[0],tail_idx[0],
                             M,N,U,mu,t,adjacency_matrix,total_nn,
                             beta,eta,canonical,N_tracker[0],
                             N_zero[0],N_beta[0],last_kinks[0],
                             dummy_counter,dummy_counter,rng);
              }
              else if (label==9){ // insert kink after head
                  insert_kink_after_head(paths[0],num_kinks[0],
                             head_idx[0],tail_idx[0],
                             M,N,U,mu,t,adjacency_matrix,total_nn,
                             beta,eta,canonical,N_tracker[0],
                             N_zero[0],N_beta[0],last_kinks[0],
                             dummy_counter,dummy_counter,rng);
              }
              else if (label==10){ // delete kink after head
                  delete_kink_after_head(paths[0],num_kinks[0],
                             head_idx[0],tail_idx[0],
                             M,N,U,mu,t,adjacency_matrix,total_nn,
                             beta,eta,canonical,N_tracker[0],
                             N_zero[0],N_beta[0],last_kinks[0],
                             dummy_counter,dummy_counter,rng);
                      }
              else if (label==11){ // insert kink before tail
                  insert_kink_before_tail(paths[0],num_kinks[0],
                             head_idx[0],tail_idx[0],
                             M,N,U,mu,t,adjacency_matrix,total_nn,
                             beta,eta,canonical,N_tracker[0],
                             N_zero[0],N_beta[0],last_kinks[0],
                             dummy_counter,dummy_counter,rng);
              }
              else if (label==12){ // delete kink before tail
                  delete_kink_before_tail(paths[0],num_kinks[0],
                             head_idx[0],tail_idx[0],
                             M,N,U,mu,t,adjacency_matrix,total_nn,
                             beta,eta,canonical,N_tracker[0],
                             N_zero[0],N_beta[0],last_kinks[0],
                             dummy_counter,dummy_counter,rng);
              }
              else if (label==13){ // insert kink after tail
                   insert_kink_after_tail(paths[0],num_kinks[0],
                              head_idx[0],tail_idx[0],
                              M,N,U,mu,t,adjacency_matrix,total_nn,
                              beta,eta,canonical,N_tracker[0],
                              N_zero[0],N_beta[0],last_kinks[0],
                              dummy_counter,dummy_counter,rng);
               }
               else if (label==14){ // delete kink after tail
                   delete_kink_after_tail(paths[0],num_kinks[0],
                              head_idx[0],tail_idx[0],
                              M,N,U,mu,t,adjacency_matrix,total_nn,
                              beta,eta,canonical,N_tracker[0],
                              N_zero[0],N_beta[0],last_kinks[0],
                              dummy_counter,dummy_counter,rng);
               }
              else{
                  // lol
              }

            // Measure the total number of particles
            if (m%(sweep*measurement_frequency)==0 && m>=0.25*sweeps_pre){
                measurement_attempts[0]+=1;
                if (head_idx[0]==-1 && tail_idx[0]==-1){
                    N_data.push_back(N_beta[0]);
                    Z_frac+=1.0;
                }
            }
            
            // Measure the number of flats
            N_flats_mean+=num_kinks[0];
            N_flats_samples+=1;
        }

        // Calculate diagonal fraction of Monte Carlo just exited
        Z_frac/=measurement_attempts[0];

        // If we did not collect data, decrease eta and try again.
        if (N_data.size()<5){eta*=0.5;continue;}

        // Find the minimum and maximum number of particles measured
        N_min=*min_element(N_data.begin(),N_data.end());
        N_max=*max_element(N_data.begin(),N_data.end());

        // Generate the support of the distribution & initialize the histogram
        N_target_in_bins=false;
        for (int i=N_min;i<=N_max;i++){
            N_bins.push_back(i);
            N_hist.push_back(0);
            P_N.push_back(0);
            if (i==N){N_target_in_bins=true;}
        }

        // Get the index of the target N in the support of the distribution
        N_idx = N-N_min;

        // Fill out the histogram
        for (int i=0;i<N_data.size();i++){
            N_hist[N_data[i]-N_min]+=1;
            N_hist_sum+=1.0;
        }

        // Build the normalized probability distribution P(N) & find its peak
        peak_idx=0;
        P_N_peak=P_N[peak_idx];
        for (int i=0;i<P_N.size();i++){
            P_N[i]=N_hist[i]/N_hist_sum;
            if (P_N[i]>P_N_peak){
                peak_idx=i;
                P_N_peak=P_N[peak_idx];
            }
        }

        // Print out current mu and draw the particle probability distribution
        cout << "mu: " << mu;
        cout << " eta: " << eta << " Z-frac: " << Z_frac*100 << "%" << endl;
        cout << "N     P(N)"<<endl;
        for (int i=0;i<N_bins.size();i++){
            cout << setw(6) << left << N_bins[i];
            for (int j=0;j<=static_cast<int>(100*P_N[i]);j++){
                cout<<"*";
            }
            cout<<endl;
            N_mean_pre+=N_bins[i]*P_N[i];
        }
        cout << "<N>: " << N_mean_pre << endl << endl;
        
        // Get average number of flats (relevant for eta-equilibration)
        N_flats_mean/=N_flats_samples;
        eta=1/sqrt(N_flats_mean);

        if (N_target_in_bins){
            // Stop the loop if the peak is at P(N)
            if (peak_idx==N_idx && abs(N_mean_pre-N)<0.33
                && at_least_one_iteration){break;}

            else{
                // Estimate mu via Eq. 15 in:https://arxiv.org/pdf/1312.6177.pdf
                if (std::count(N_bins.begin(),N_bins.end(),N-1) &&
                    std::count(N_bins.begin(),N_bins.end(),N+1)){
                    mu_right=mu-(1/beta)*log(P_N[N_idx+1]/P_N[N_idx]);
                    mu_left=mu-(1/beta)*log(P_N[N_idx]/P_N[N_idx-1]);
                    mu=0.5*(mu_left+mu_right);
                }
                else if (std::count(N_bins.begin(),N_bins.end(),N+1)){
                    mu_right=mu-(1/beta)*log(P_N[N_idx+1]/P_N[N_idx]);
                    mu=mu_right;
                }
                else if (std::count(N_bins.begin(),N_bins.end(),N-1)){
                    mu_left=mu-(1/beta)*log(P_N[N_idx]/P_N[N_idx-1]);
                    mu=mu_left;
                }
                else { // Peak is 100% at N... Yes. It can happen.
                    // We might've entered here b.c need at least 2 iterations
                }
            }
        }
        else{ // Target N not in P_N
            if (N_bins[peak_idx]>N){
                if (mu>1){mu*=0.5;}
                else if (mu<=-1){mu*=1.1;}
                else {mu-=1;}
            }
            else{
                if (mu>1){mu*=1.1;}
                else if (mu<=-1){mu*=0.5;}
                else {mu+=1;}
            }
        }
        at_least_one_iteration=true;
    }
    
/*---------------------------- Open files ------------------------------------*/
    
    // Declare conventional estimator files
    if (num_replicas<2){
        
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
        }
    }
    
    // Estimators in replicated configuration space
    else {
        string SWAP_histogram_name;
        vector<string> n_resolved_SWAP_histogram_names;
        
        if (canonical){ // name of file if canonical simulation
            SWAP_histogram_name=to_string(L)+"_"+to_string(N)+"_"+
            to_string(l_A)+"_"+to_string(D)+"D_"+
            to_string(U)+"_"+to_string(beta)+"_"+
            to_string(t)+"_"+to_string(sweeps)+"_"+
            to_string(seed)+"_"+
            "can_"+"SWAP.dat";
           
            // Create filenames of SWAP histograms for each n-Sector
            for (int i=0; i<=N; i++){
                n_resolved_SWAP_histogram_names.push_back(
                to_string(L)+"_"+to_string(N)+"_"+
                to_string(l_A)+"_"+to_string(D)+"D_"+
                to_string(U)+"_"+to_string(beta)+"_"+
                to_string(t)+"_"+to_string(sweeps)+"_"+
                to_string(seed)+"_"+"can_"+"SWAP_"+
                to_string(i)+"-sector.dat");
            }
        }
        else { // name of file if grand canonical simulation
            SWAP_histogram_name=to_string(L)+"_"+to_string(N)+"_"+
            to_string(l_A)+"_"+to_string(D)+"D_"+
            to_string(U)+"_"+to_string(beta)+"_"+
            to_string(t)+"_"+to_string(sweeps)+"_"+
            to_string(seed)+"_"+
            "grandcan_"+"SWAP.dat";
            
            // Create filenames of SWAP histograms for each n-Sector
            for (int i=0; i<=N; i++){
                n_resolved_SWAP_histogram_names.push_back(
                to_string(L)+"_"+to_string(N)+"_"+
                to_string(l_A)+"_"+to_string(D)+"D_"+
                to_string(U)+"_"+to_string(beta)+"_"+
                to_string(t)+"_"+to_string(sweeps)+"_"+
                to_string(seed)+"_"+"grandcan_"+"SWAP_"+
                to_string(i)+"-sector.dat");
            }
        }
            
        // Open SWAP histograms file
        SWAP_histogram_file.open(SWAP_histogram_name);
        for (int i=0; i<=N; i++){
            n_resolved_SWAP_histogram_file.open(
                                            n_resolved_SWAP_histogram_names[i]);
            n_resolved_SWAP_histogram_files.push_back(std::move(
                                               n_resolved_SWAP_histogram_file));
        }

        if( !SWAP_histogram_file ) { // file couldn't be opened
           cerr << "Error: SWAP histogram file could not be opened" << endl;
           exit(1);
        }

    } // End of replicated estimators else block
    
/*---------------------------- Monte Carlo -----------------------------------*/
    
    // Time main function execution
    auto start = high_resolution_clock::now();
    
    // Restart data structure and trackers
    num_kinks.clear();
    N_tracker.clear();
    head_idx.clear();
    tail_idx.clear();
    N_zero.clear();
    N_beta.clear();
    last_kinks.clear();
    paths.clear();
    for (int r=0;r<num_replicas;r++){
        num_kinks.push_back(M);
        N_tracker.push_back(N);
        head_idx.push_back(-1);
        tail_idx.push_back(-1);
        N_zero.push_back(N);
        N_beta.push_back(N);
        
        last_kinks.push_back(vector<int> (M,-1));
        for (int i=0;i<M;i++){last_kinks[r][i]=i;}
        
        paths.push_back(create_paths(initial_fock_state,M,r));
    }

    Z_frac=0.0;
    std::fill(measurement_attempts.begin(),measurement_attempts.end(),0);
    
    if (beta>=1.0){sweeps*=(beta*M);}
    else {sweeps*=M;}
    
    // Initialize vector estimators (conventional or SWAP)
    if (num_replicas<2) { // conventional vector estimators
        for (int r=0;r<num_replicas;r++){
            tr_kinetic_energy.push_back(vector<double>
                                        (measurement_centers.size(),0.0));
            tr_diagonal_energy.push_back(vector<double>
                                         (measurement_centers.size(),0.0));
        }
    }
    else { // SWAP vector estimators
        for (int i=0; i<=m_A; i++){
            SWAP_histogram.push_back(0); // just initializing
        }
        for (int i=0; i<=N; i++){
            n_resolved_SWAP_histograms.push_back(vector<int> (m_A+1,0));
        }
    }

    cout << "Stage (2/3): Equilibrating..." << endl << endl;
    
    for (unsigned long long int m=0; m < sweeps; m++){
    for (int r=0;r<num_replicas;r++){
        
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
    } // end of replica loop

        
        // SWAP Updates
        label = swap_updates(rng);

        if (label==0){ // insert_swap_kink
         insert_swap_kink(paths, num_kinks,
                         num_replicas, 0,
                         sub_sites, swapped_sites,
                         swap_kinks, num_swaps,
                         l_A, m_A,
                         head_idx,tail_idx,
                         M, N, U, mu, t,
                         adjacency_matrix, total_nn,
                         beta, eta, canonical, N_tracker,
                         N_zero, N_beta,
                         last_kinks,
                         insert_swap_kink_attempts,
                         insert_swap_kink_accepts,
                         rng);
         }
         else if (label==1) { // delete_swap_kink
             delete_swap_kink(paths, num_kinks,
                             num_replicas, 0,
                             sub_sites, swapped_sites,
                             swap_kinks, num_swaps,
                             l_A, m_A,
                             head_idx,tail_idx,
                             M, N, U, mu, t,
                             adjacency_matrix, total_nn,
                             beta, eta, canonical, N_tracker,
                             N_zero, N_beta,
                             last_kinks,
                             delete_swap_kink_attempts,
                             delete_swap_kink_accepts,
                             rng);
         }
         else if (label==2) {
             swap_timeshift_head(paths, num_kinks,
                             num_replicas, 0,
                             sub_sites, swapped_sites,
                             swap_kinks, num_swaps,
                             l_A, m_A,
                             head_idx,tail_idx,
                             M, N, U, mu, t,
                             adjacency_matrix, total_nn,
                             beta, eta, canonical, N_tracker,
                             N_zero, N_beta,
                             last_kinks,
                             swap_advance_head_attempts,
                             swap_advance_head_accepts,
                             swap_recede_head_attempts,
                             swap_recede_head_accepts,
                             rng);
         }
         else{
             // lol
         }
        
/*----------------------------- Measurements ---------------------------------*/
        
        if (m%(sweep*measurement_frequency)==0 && m>=sweeps*0.25){
 
            if (not_equilibrated){
                not_equilibrated=false;
                cout << "Stage (3/3): Main Monte Carlo loop..." << endl;
            }
            
            // Conventional measurements
            if (num_replicas<2){ // conventional measurements
                
                int r=0; // TEMPORARY (eventually might loop over 0,1)
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
            } // end of conventional measurements if statement
            
            // Non-conventional (SWAP) measurements
            if (num_replicas>1) {
                // add count to bin corresponding to number of swapped sites
                if (head_idx[0]==-1&&head_idx[1]==-1&&
                    tail_idx[0]==-1&&tail_idx[1]==-1){
                    if (N_zero[0]==N && N_beta[0]==N
                        && N_zero[1]==N && N_beta[1]==N){

                        // Add count to histogram of number of swapped sites
                        SWAP_histogram[num_swaps]+=1;

                        // Add count to swapped sites histogram of n-sector
                        n_A = 0;
                        get_fock_state(beta/2.0,M,fock_state_at_half_plus,paths[0]);
                        for (int i=0; i<m_A; i++){
                            n_A+=fock_state_at_half_plus[sub_sites[i]];
                        }
                        n_resolved_SWAP_histograms[n_A][num_swaps]+=1;
                        
                        // Add a count to measurements in the bin
                        writing_ctr+=1;
                    }
                }
            
                if (writing_ctr==bin_size){

                    // Save current histogram of swapped sites to file
                    for (int i=0; i<=m_A; i++){
                        SWAP_histogram_file<<fixed<<setprecision(17)<<
                        SWAP_histogram[i] << " ";
                    }
                    SWAP_histogram_file<<endl;

                    // Restart histogram
                    std::fill(SWAP_histogram.begin(),
                              SWAP_histogram.end(),0);
                    
                    // Save current n-resolved swapped sites histogram to file
                    for (int i=0; i<=N; i++){
                        for (int j=0; j<=m_A; j++){
                            n_resolved_SWAP_histogram_files[i]<<
                            fixed<<setprecision(17)<<
                            n_resolved_SWAP_histograms[i][j]<<" ";
                        }
                        n_resolved_SWAP_histogram_files[i]<<endl;
                        
                        // Restart histogram
                        std::fill(n_resolved_SWAP_histograms[i].begin(),
                                  n_resolved_SWAP_histograms[i].end(),0);
                        
                    }
                    


                    // Restart counter that tracks when to save to file
                    writing_ctr=0;
                }

            } // end of SWAP measurements if statement
        } // end of measurement after 25% equilibration if statement
    } // end of sweeps loop
    
    // Close data files
    if (num_replicas<2){
        for (int r=0;r<num_replicas;r++){
            kinetic_energy_file[r].close();
            diagonal_energy_file[r].close();
            tr_kinetic_energy_file[r].close();
            tr_diagonal_energy_file[r].close();
        }
    }
    else {
        SWAP_histogram_file.close();
        for (int i=0; i<=m_A; i++){
            n_resolved_SWAP_histogram_files[i].close();
        }
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
    
    cout<< endl <<"SWAP: "<<insert_swap_kink_accepts<<"/"<<
                               insert_swap_kink_attempts<<endl;
    cout <<"UNSWAP: "<<delete_swap_kink_accepts<<"/"<<
                               delete_swap_kink_attempts<<endl;
    
    cout<< endl <<"SWAP Advance Head: "<<swap_advance_head_accepts<<"/"<<
                               swap_advance_head_attempts<<endl;
    cout <<"SWAP Recede Head: "<<swap_recede_head_accepts<<"/"<<
                               swap_recede_head_attempts<<endl;
    
    auto end = high_resolution_clock::now();

    auto elapsed_time = duration_cast<nanoseconds>(end - start);
    double duration = elapsed_time.count() * 1e-9;
    
    cout << endl << "beta: " << beta << endl;
    cout << endl << "sweeps: " << sweeps/(beta*M) << endl;
    
    if (num_replicas<2){
    cout << "Z_ctr: " << Z_ctr[0] << endl;
    cout<<"Z_frac: "<<Z_ctr[0]*100.0/measurement_attempts[0]<<"% ("<<Z_ctr[0]
    <<"/"<< measurement_attempts[0]<<")"<<endl;
    
    cout << endl << "<N>: " << (N_sum[0])/Z_ctr[0] << endl;
    }

    cout << endl << "Elapsed time: " << duration << " seconds" << endl;

    return 0;
    
}
