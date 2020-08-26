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
    
    // Time main function execution
    auto start = high_resolution_clock::now();
    
    // Create a uniform distribution with support: [0.0,1.0)
    boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    
    // Bose-Hubbard parameters
    int L = 3, D = 2, N = pow(L,D);
    double t = 1.0, U = 5.0, mu = -2.6019;
    vector<int> alpha;
    int M = pow(L,D); // total sites
    string boundary_condition = "pbc";
    
    // Simulation parameters
    double eta = 0.03865, beta = 3.0;
    bool canonical = true;
    unsigned long long int sweeps=1000000,sweep=M*beta;
    
    // Adjacency matrix
    int total_nn = count_hypercube_nearest_neighbors(L,D,boundary_condition);
    vector<int> adjacency_matrix_rows (total_nn,0);
    vector<vector<int>> adjacency_matrix (M,adjacency_matrix_rows);
    build_hypercube_adjacency_matrix(L,D,boundary_condition,adjacency_matrix);
    
    // Trackers
    int num_kinks = M;
    double N_tracker = N;
    int head_idx = -1, tail_idx = -1;
    int N_zero=N, N_beta=N;
    vector<int> last_kinks (M,-1);
    
    // Attempt/Acceptance counters
    unsigned long long int  insert_worm_attempts=0, insert_worm_accepts=0;
    unsigned long long int  delete_worm_attempts=0, delete_worm_accepts=0;

    unsigned long long int  insert_anti_attempts=0, insert_anti_accepts=0;
    unsigned long long int  delete_anti_attempts=0, delete_anti_accepts=0;
    
    unsigned long long int  insertZero_worm_attempts=0, insertZero_worm_accepts=0;
    unsigned long long int  deleteZero_worm_attempts=0, deleteZero_worm_accepts=0;

    unsigned long long int  insertZero_anti_attempts=0, insertZero_anti_accepts=0;
    unsigned long long int  deleteZero_anti_attempts=0, deleteZero_anti_accepts=0;
    
    unsigned long long int  insertBeta_worm_attempts=0, insertBeta_worm_accepts=0;
    unsigned long long int  deleteBeta_worm_attempts=0, deleteBeta_worm_accepts=0;

    unsigned long long int  insertBeta_anti_attempts=0, insertBeta_anti_accepts=0;
    unsigned long long int  deleteBeta_anti_attempts=0, deleteBeta_anti_accepts=0;
    
    unsigned long long int  advance_head_attempts=0, advance_head_accepts=0;
    unsigned long long int  recede_head_attempts=0, recede_head_accepts=0;
    
    unsigned long long int  advance_tail_attempts=0, advance_tail_accepts=0;
    unsigned long long int  recede_tail_attempts=0, recede_tail_accepts=0;
    
    unsigned long long int  ikbh_attempts=0, ikbh_accepts=0;
    unsigned long long int  dkbh_attempts=0, dkbh_accepts=0;
    
    unsigned long long int  ikah_attempts=0, ikah_accepts=0;
    unsigned long long int  dkah_attempts=0, dkah_accepts=0;
    
    unsigned long long int  ikbt_attempts=0, ikbt_accepts=0;
    unsigned long long int  dkbt_attempts=0, dkbt_accepts=0;
    
    unsigned long long int  ikat_attempts=0, ikat_accepts=0;
    unsigned long long int  dkat_attempts=0, dkat_accepts=0;
    
    // Observables
    double N_sum=0,kinetic_energy=0,diagonal_energy=0;
    double measurement_center=beta/2,measurement_plus_minus=0.1;
    vector<int> fock_state_at_slice (M,0);
    
    // Non-observables
    unsigned long long int Z_ctr=0,measurement_attempts=0;
    
    // Generate a random fock state
    alpha = random_boson_config(M,N);

    // Generate the data structure (vector of kinks)
    vector<Kink> kinks_vector = create_kinks_vector(alpha, M);
    
    // Initialize array containing indices of last kinks at each site
    for (int i=0; i<M; i++){
        last_kinks[i] = i;
    }
    
    // Initialize vector that will store nearest neighbor indices
    vector<int> nn_sites (total_nn,0);
    
    cout << "Lattice PIGS started: " << endl << endl;
    
/*---------------------------- Open files ------------------------------------*/

    
/*---------------------------- Monte Carlo -----------------------------------*/

    boost::random::uniform_int_distribution<> updates(0, 14);
    int label;
    sweeps *= (beta*M);
    for (unsigned long long int m=0; m < sweeps; m++){
        label = updates(rng);
        
        if (label==0){     // worm_insert
            insert_worm(kinks_vector,num_kinks,head_idx,tail_idx,
                        M,N,U,mu,t,beta,eta,canonical,N_tracker,
                        N_zero, N_beta, last_kinks,
                        insert_worm_attempts,insert_worm_accepts,
                        insert_anti_attempts,insert_anti_accepts);
        }
        else if (label==1){ // worm_delete
            delete_worm(kinks_vector,num_kinks,head_idx,tail_idx,
                        M,N,U,mu,t,beta,eta,canonical,N_tracker,
                        N_zero, N_beta, last_kinks,
                        delete_worm_attempts,delete_worm_accepts,
                        delete_anti_attempts,delete_anti_accepts);
        }
        else if (label==2){ // insertZero
            insertZero(kinks_vector,num_kinks,head_idx,tail_idx,
                       M,N,U,mu,t,beta,eta,canonical,N_tracker,
                       N_zero, N_beta, last_kinks,
                       insertZero_worm_attempts,insertZero_worm_accepts,
                       insertZero_anti_attempts,insertZero_anti_accepts);
            
        }
        else if (label==3){ // deleteZero
            deleteZero(kinks_vector,num_kinks,head_idx,tail_idx,
                       M,N,U,mu,t,beta,eta,canonical,N_tracker,
                       N_zero, N_beta, last_kinks,
                       deleteZero_worm_attempts,deleteZero_worm_accepts,
                       deleteZero_anti_attempts,deleteZero_anti_accepts);
        }
        else if (label==4){ // insertBeta
            insertBeta(kinks_vector,num_kinks,head_idx,tail_idx,
                       M,N,U,mu,t,beta,eta,canonical,N_tracker,
                       N_zero, N_beta, last_kinks,
                       insertBeta_worm_attempts,insertBeta_worm_accepts,
                       insertBeta_anti_attempts,insertBeta_anti_accepts);
        }
        else if (label==5){ // deleteBeta
            deleteBeta(kinks_vector,num_kinks,head_idx,tail_idx,
                       M,N,U,mu,t,beta,eta,canonical,N_tracker,
                       N_zero, N_beta, last_kinks,
                       deleteBeta_worm_attempts,deleteBeta_worm_accepts,
                       deleteBeta_anti_attempts,deleteBeta_anti_accepts);
        }
        else if (label==6){ // timeshift
            timeshift(kinks_vector,num_kinks,head_idx,tail_idx,
                       M,N,U,mu,t,beta,eta,canonical,N_tracker,
                       N_zero, N_beta, last_kinks,
                       advance_head_attempts, advance_head_accepts,
                       recede_head_attempts, recede_head_accepts,
                       advance_tail_attempts, advance_tail_accepts,
                       recede_tail_attempts, recede_tail_accepts);
        }
        else if (label==7){ // insert kink before head
            insert_kink_before_head(kinks_vector,num_kinks,head_idx,tail_idx,
                       M,N,U,mu,t,adjacency_matrix,total_nn,
                       beta,eta,canonical,N_tracker,
                       N_zero, N_beta, last_kinks,
                       ikbh_attempts, ikbh_accepts);
        }
        else if (label==8){ // delete kink before head
            delete_kink_before_head(kinks_vector,num_kinks,head_idx,tail_idx,
                       M,N,U,mu,t,adjacency_matrix,total_nn,
                       beta,eta,canonical,N_tracker,
                       N_zero, N_beta, last_kinks,
                       dkbh_attempts, dkbh_accepts);
        }
        else if (label==9){ // insert kink after head
            insert_kink_after_head(kinks_vector,num_kinks,head_idx,tail_idx,
                       M,N,U,mu,t,adjacency_matrix,total_nn,
                       beta,eta,canonical,N_tracker,
                       N_zero, N_beta, last_kinks,
                       ikah_attempts, ikah_accepts);
        }
        else if (label==10){ // delete kink after head
            delete_kink_after_head(kinks_vector,num_kinks,head_idx,tail_idx,
                       M,N,U,mu,t,adjacency_matrix,total_nn,
                       beta,eta,canonical,N_tracker,
                       N_zero, N_beta, last_kinks,
                       dkah_attempts, dkah_accepts);
                }
        else if (label==11){ // insert kink before tail
            insert_kink_before_tail(kinks_vector,num_kinks,head_idx,tail_idx,
                       M,N,U,mu,t,adjacency_matrix,total_nn,
                       beta,eta,canonical,N_tracker,
                       N_zero, N_beta, last_kinks,
                       ikbt_attempts, ikbt_accepts);
        }
        else if (label==12){ // delete kink before tail
            delete_kink_before_tail(kinks_vector,num_kinks,head_idx,tail_idx,
                       M,N,U,mu,t,adjacency_matrix,total_nn,
                       beta,eta,canonical,N_tracker,
                       N_zero, N_beta, last_kinks,
                       dkbt_attempts, dkbt_accepts);
        }


        else if (label==13){ // insert kink after tail
             insert_kink_after_tail(kinks_vector,num_kinks,head_idx,tail_idx,
                        M,N,U,mu,t,adjacency_matrix,total_nn,
                        beta,eta,canonical,N_tracker,
                        N_zero, N_beta, last_kinks,
                        ikat_attempts, ikat_accepts);
         }
         else if (label==14){ // delete kink after tail
             delete_kink_after_tail(kinks_vector,num_kinks,head_idx,tail_idx,
                        M,N,U,mu,t,adjacency_matrix,total_nn,
                        beta,eta,canonical,N_tracker,
                        N_zero, N_beta, last_kinks,
                        dkat_attempts, dkat_accepts);
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

        
        if (m%sweep==0 && m>0.2*sweeps){
            measurement_attempts+=1;
            if (head_idx==-1 and tail_idx==-1 && N_beta==N){
                                
                // Get fock state at desired measurement center
                get_fock_state(measurement_center,M,fock_state_at_slice,
                               kinks_vector);
                
                // Measure N
                N_sum += N_tracker;
                Z_ctr += 1;
                
                // Measure <K>
                kinetic_energy += pimc_kinetic_energy(kinks_vector,num_kinks,
                                    measurement_center,measurement_plus_minus,
                                    M,t,beta);
                
//                // Measure <V>
                diagonal_energy += pimc_diagonal_energy(fock_state_at_slice,
                                                        M,U,mu);
//                V=0;
//                int ctr=0;
//                for (int site=0;site<M;site++){
//                    while kinks_vector
//                }
                
            }
        }
    }
        
/*--------------------------------- FIN --------------------------------------*/

    // Print out the head and tail indices
    cout << "head_idx: " << head_idx << endl;
    cout << "tail_idx: " << tail_idx << endl;

    // Print out the N_tracker
    cout << "N_tracker: " << N_tracker << endl;
    
    // Print out number of particles at path ends
    cout << "N_zero: " << N_zero << endl;
    cout << "N_beta: " << N_beta << endl;

    // Print out number of active kinks
    cout << "num_kinks: " << num_kinks << endl;

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
    
    cout << endl << "sweeps: " << sweeps/(M*beta) << endl;
    cout << "Z_ctr: " << Z_ctr << endl;
    cout << "Z_frac: " << Z_ctr*100.0/measurement_attempts << "% (" << Z_ctr
    << "/" << measurement_attempts << ")" << endl;
    
    cout << endl << "<N>: " << (N_sum)/Z_ctr << endl;
    cout << "<E>: " << (kinetic_energy+diagonal_energy)/Z_ctr+mu*N << endl;
    cout << "<K>: " << kinetic_energy/Z_ctr << endl;
    cout << "<V>: " << diagonal_energy/Z_ctr + mu*N << endl;

    cout << endl << "Elapsed time: " << duration << " seconds" << endl;

    cout << endl;
    cout << "Total neighbors: " << total_nn << endl << endl;
    
    return 0;
}

