//
//  pimc.hpp
//  pimc
//
//  Created by Emanuel Casiano-Diaz on 8/22/20.
//  Copyright © 2020 Emanuel Casiano-Díaz. All rights reserved.
//

#ifndef pimc_hpp
#define pimc_hpp

#include <stdio.h>
#include<iostream>
#include<vector>
// #include<boost/random.hpp>
#include<cmath>
#include<chrono>
#include<iomanip>  // for std::setprecision
#include <fstream>
#include <cstdlib> // for exit function
#include "RNG.h"
#include<sstream>
#include<string.h>

#include <boost/math/special_functions/lambert_w.hpp>
using boost::math::lambert_w0;
using boost::math::lambert_wm1;

using namespace std;
using namespace std::chrono;

/*--------------------------- Class Definitions ------------------------------*/

class Kink
{
    public:
    // Attribute declarations
    double tau;
    int n,src,dest,prev,next,src_replica,dest_replica;
    
    // Member function declarations (prototypes)
    Kink (double,int,int,int,int,int,int,int); // Kink constructor
    
    // Make "<<" a friend of the Kink class
    friend ostream& operator<<(ostream& os, const Kink& dt);
};

// Member function definitions
Kink::Kink (double time,int particles,int source_site,int destination_site,
            int prev_kink_idx,int next_kink_idx,
            int source_replica,int destination_replica){
    tau = time;
    n = particles;
    src = source_site;
    dest = destination_site;
    prev = prev_kink_idx;
    next = next_kink_idx;
    src_replica = source_replica;
    dest_replica = destination_replica;
    // bool swap_kink -> Add this
    // int replica_idx
}

// Overload "<<" operator
ostream& operator<<(ostream& os, const Kink& dt)
{
//    os << '<' << dt.tau << ',' << dt.n << ',' << dt.src << ','
//    << dt.dest << ',' << dt.prev << ',' << dt.next <<'>'
//    << dt.src_replica << dt.dest_replica;
    
        os << setprecision(17) << dt.tau << ' ' << dt.n << ' ' << dt.src << ' '
        << dt.dest << ' ' << dt.prev << ' ' << dt.next << ' '
        << dt.src_replica << ' ' << dt.dest_replica;

    return os;
}

/*-------------------------- Function Definitions ----------------------------*/

// function to convert decimal to binary
void decToBinary(int n)
{
    // Size of an integer is assumed to be 32 bits
    for (int i = 31; i >= 0; i--) {
        int k = n >> i;
        if (k & 1)
            cout << "1";
        else
            cout << "0";
    }
    return;
}

/*--------------------------------------------------------------------*/

// function to convert decimal to binary
int binaryToDecimal(vector<int> binary_word){
    int decimal=0;
    bool bit;
    
    for (size_t i=0;i<binary_word.size();i++){
        bit = binary_word[i];
        if (bit){
            decimal += pow(2,binary_word.size()-(i+1));
        }
    }

    return decimal;
}

/*--------------------------------------------------------------------*/

vector<int> random_boson_config(int M,int N,RNG &rng,bool restart){
    // Generates random Fock state of N bosons in M=L^D sites
    
    vector<int> alpha (M,0);
    unsigned long src;
    
    // Initialize the distribution object. Note: Support is fully closed [0,M-1]
    //boost::random::uniform_int_distribution<> sites(0, M-1);
    
    if (!restart){
        // Randomly sprinkle the N particles among sites
        for (int n=1; n<=N; n++){
            src = rng.randInt(M-1);
            alpha[src] += 1;
        }
    }
    else{ //restart
        for (int n=1; n<=N; n++){
            src = 1; // Throw all particles on first site for now.
            alpha[src] += 1;
        }
    }
    
    return alpha;
}

/*--------------------------------------------------------------------*/

vector<Kink> create_paths(vector<int> &fock_state,int M,int replica_idx){

    int num_empty_kinks;

    // Set the number of kinks to pre-allocate based on lattice size
    num_empty_kinks = M*1000;

    // Pre-allocate kinks. Recall: (tau,n,src,dest,prev,next)
    vector<Kink> paths(num_empty_kinks,Kink(-1.0,-1,-1,-1,-1,-1,-1,-1));

    // Initialize the first M=L^D kinks
    for (int site=0; site<M; site++){
        paths[site] = Kink(0.0,fock_state[site],site,site,-1,-1,
                           replica_idx,replica_idx);
    }
    return paths;
}

/*--------------------------------------------------------------------*/

ofstream save_paths(int D, int L, int N, int l_A,
                       double U, double t, double beta,
                       int bin_size, int bins_wanted, int seed,
                       string subgeometry, double mu, double eta,
                       int num_replicas, vector<int> num_kinks,
                       vector<vector<Kink> > paths,
                       vector<double> N_tracker,
                       unsigned long long int iteration_idx){
    
    // Saving last worldline configuration
    ofstream state_file;
    string state_name;
    
    // Name of system state file
    state_name=to_string(D)+"D_"+to_string(L)+
    "_"+to_string(N)+"_"+to_string(l_A)+"_"+
    to_string(U)+"_"+to_string(t)+"_"+
    to_string(beta)+"_"+to_string(bin_size)+"_"+
    "system-state_"+to_string(seed)+"_"+subgeometry+"_"+
    to_string(num_replicas)+".dat";
    
    state_file.open(state_name);
    
    // Find how many active kinks are in the replica with the most
    // active kinks
    int max_num_kinks=-1;
    for (int r=0; r<num_replicas; r++){
        if (num_kinks[r]>max_num_kinks){max_num_kinks=num_kinks[r];}
    }
    
    // First 8 attributes: first replica
    // Next 8 attributes: second replicas
    for (int k=0; k<max_num_kinks; k++){
        for (int r=0; r<num_replicas; r++){

            if (k<num_kinks[r]){
            state_file<<fixed<<setprecision(17)<<paths[r][k].tau<<" ";
            state_file<<fixed<<paths[r][k].n<<" ";
            state_file<<fixed<<paths[r][k].src<<" ";
            state_file<<fixed<<paths[r][k].dest<<" ";
            state_file<<fixed<<paths[r][k].prev<<" ";
            state_file<<fixed<<paths[r][k].next<<" ";
            state_file<<fixed<<paths[r][k].src_replica<<" ";
            state_file<<fixed<<paths[r][k].dest_replica<<" ";
            }
            else {
            state_file<<fixed<<setprecision(17)<<-1.0<<" ";
            state_file<<fixed<<-1<<" ";
            state_file<<fixed<<-1<<" ";
            state_file<<fixed<<-1<<" ";
            state_file<<fixed<<-1<<" ";
            state_file<<fixed<<-1<<" ";
            state_file<<fixed<<-1<<" ";
            state_file<<fixed<<-1<<" ";
            }
        }
        state_file<<fixed<<setprecision(17)<<mu<<" ";
        state_file<<fixed<<setprecision(17)<<eta<<" ";
        state_file<<fixed<<iteration_idx<<" ";
        
        state_file<<endl;
    }
    return state_file;
}

/*--------------------------------------------------------------------*/

vector<vector<Kink> > load_paths(int D, int L, int N, int l_A,
                                 double U, double t, double beta,
                                 int bin_size, int bins_wanted,
                                 int seed, string subgeometry,
                                 int num_replicas){

    int M = pow(L,D);
    string state_name;

    // Name of system state file
    state_name=to_string(D)+"D_"+to_string(L)+
    "_"+to_string(N)+"_"+to_string(l_A)+"_"+
    to_string(U)+"_"+to_string(t)+"_"+
    to_string(beta)+"_"+to_string(bin_size)+"_"+
    "system-state_"+to_string(seed)+"_"+subgeometry+"_"+
    to_string(num_replicas)+".dat";
    
    // File containing previous worldline configurations (paths)
    std::ifstream infile(state_name);
    
    if (!infile) {
            cout << "Unable to open file";
            exit(1); // terminate with error
        }

    // For each replica, initialize vector containinig all kinks
    vector<vector<Kink> > paths(num_replicas,vector<Kink> (M*1000,Kink(-1.0,-1,-1,-1,-1,-1,-1,-1)));
    
    if (num_replicas==2){
        // Assume 16 elements per line (two replicas max)
        int a1,a2,a3,a4,a5,a6,a7,a9,a10,a11,a12,a13,a14,a15;
        double a0,a8,a16,a17;
        unsigned long long int a18;
        
        int k=0; // kink idx counter
                
        while(infile >> a0 >> a1 >> a2 >> a3 >> a4 >> a5 >> a6 >> a7 >> a8 >> a9 >> a10 >> a11 >> a12 >> a13 >> a14 >> a15 >>
            a16 >> a17 >> a18){
                        
            // Fill paths of first replica
            paths[0][k].tau = a0;
            paths[0][k].n = a1;
            paths[0][k].src = a2;
            paths[0][k].dest = a3;
            paths[0][k].prev = a4;
            paths[0][k].next = a5;
            paths[0][k].src_replica = a6;
            paths[0][k].dest_replica = a7;
            
            // Fill paths of second replica
            paths[1][k].tau = a8;
            paths[1][k].n = a9;
            paths[1][k].src = a10;
            paths[1][k].dest = a11;
            paths[1][k].prev = a12;
            paths[1][k].next = a13;
            paths[1][k].src_replica = a14;
            paths[1][k].dest_replica = a15;
            
            k++;
        }
    }
    
    if (num_replicas==1){
    
        // Assume 16 elements per line (two replicas max)
        int a1,a2,a3,a4,a5,a6,a7;
        double a0,a8,a9;
        unsigned long long int a10;
        int k=0; // kink idx counter
        
        while(infile >> a0 >> a1 >> a2 >> a3 >> a4 >> a5 >> a6 >> a7 >>
              a8 >> a9 >> a10){
                
            // Fill paths of first replica
            paths[0][k].tau = a0;
            paths[0][k].n = a1;
            paths[0][k].src = a2;
            paths[0][k].dest = a3;
            paths[0][k].prev = a4;
            paths[0][k].next = a5;
            paths[0][k].src_replica = a6;
            paths[0][k].dest_replica = a7;
            
            k++;
        }
    }
    
    infile.close();
        
    return paths;
}

/*--------------------------------------------------------------------*/

vector<int> get_num_kinks(int D, int L, int N, int l_A,
                          double U, double t, double beta,
                          int bin_size, int bins_wanted,
                          int seed, string subgeometry,
                          int num_replicas){
        
    string state_name;

    // Name of system state file
    state_name=to_string(D)+"D_"+to_string(L)+
    "_"+to_string(N)+"_"+to_string(l_A)+"_"+
    to_string(U)+"_"+to_string(t)+"_"+
    to_string(beta)+"_"+to_string(bin_size)+"_"+
    "system-state_"+to_string(seed)+"_"+subgeometry+"_"+
    to_string(num_replicas)+".dat";
    
    // NOTE: For consistency, may rewrite function to get the number
    // of kinks from the path structure created with load_paths()
    std::ifstream infile(state_name);
    
    // Stores number of kinks in each replica
    vector<int> num_kinks(num_replicas,0);
    
    if (num_replicas==2){
        // Assume 16 elements per line (two replicas max)
        int a1,a2,a3,a4,a5,a6,a7,a9,a10,a11,a12,a13,a14,a15;
        double a0,a8,a16,a17;
        unsigned long long int a18;
        
        while(infile >> a0 >> a1 >> a2 >> a3 >> a4 >> a5 >> a6 >> a7 >> a8 >> a9 >> a10 >> a11 >> a12 >> a13 >> a14 >> a15 >>
            a16 >> a17 >> a18){
                
            // -1 elements indicate inactive kinks. Count only actives
            if (a1!=-1){num_kinks[0]+=1;}
            if (a9!=-1){num_kinks[1]+=1;}
        }
    }
    
    if (num_replicas==1){
        // Assume 16 elements per line (two replicas max)
        int a1,a2,a3,a4,a5,a6,a7;
        double a0,a8,a9;
        unsigned long long int a10;

        while(infile >> a0 >> a1 >> a2 >> a3 >> a4 >> a5 >> a6 >> a7 >>
              a8 >> a9 >> a10){
            // -1 elements indicate inactive kinks. Count only actives
            if (a1!=-1){num_kinks[0]+=1;}
        }
    }
    
    infile.close();
    
    return num_kinks;
}

/*--------------------------------------------------------------------*/

double get_mu(int D, int L, int N, int l_A,
              double U, double t, double beta,
              int bin_size, int bins_wanted,
              int seed, string subgeometry,
              int num_replicas){
        
    double mu;
    string state_name;
    
    // Name of system state file
    state_name=to_string(D)+"D_"+to_string(L)+
    "_"+to_string(N)+"_"+to_string(l_A)+"_"+
    to_string(U)+"_"+to_string(t)+"_"+
    to_string(beta)+"_"+to_string(bin_size)+"_"+
    "system-state_"+to_string(seed)+"_"+subgeometry+"_"+
    to_string(num_replicas)+".dat";
    
    // NOTE: For consistency, may rewrite function to get the number
    // of kinks from the path structure created with load_paths()
    std::ifstream infile(state_name);
    
    mu = -1;
    if (num_replicas==2){
        // Assume 16 elements per line (two replicas max)
        int a1,a2,a3,a4,a5,a6,a7,a9,a10,a11,a12,a13,a14,a15;
        double a0,a8,a16,a17;
        unsigned long long int a18;
        
        while(infile >> a0 >> a1 >> a2 >> a3 >> a4 >> a5 >> a6 >> a7 >> a8 >> a9 >> a10 >> a11 >> a12 >> a13 >> a14 >> a15 >>
            a16 >> a17 >> a18){
            mu = a16;
            break;
        }
    }
    
    if (num_replicas==1){
        // Assume 16 elements per line (two replicas max)
        int a1,a2,a3,a4,a5,a6,a7;
        double a0,a8,a9;
        unsigned long long int a10;
        
        while(infile >> a0 >> a1 >> a2 >> a3 >> a4 >> a5 >> a6 >> a7 >>
              a8 >> a9 >> a10){
            mu = a8;
            break;
        }
    }
    
    infile.close();
    
    return mu;
}

/*--------------------------------------------------------------------*/

unsigned long long int get_iteration_idx(int D, int L, int N, int l_A,
              double U, double t, double beta,
              int bin_size, int bins_wanted,
              int seed, string subgeometry,
              int num_replicas){
        
    unsigned long long int iteration_idx;
    string state_name;
    
    // Name of system state file
    state_name=to_string(D)+"D_"+to_string(L)+
    "_"+to_string(N)+"_"+to_string(l_A)+"_"+
    to_string(U)+"_"+to_string(t)+"_"+
    to_string(beta)+"_"+to_string(bin_size)+"_"+
    "system-state_"+to_string(seed)+"_"+subgeometry+"_"+
    to_string(num_replicas)+".dat";
        
    // NOTE: For consistency, may rewrite function to get the number
    // of kinks from the path structure created with load_paths()
    std::ifstream infile(state_name);
    
    iteration_idx = -1;
    if (num_replicas==2){
        // Assume 16 elements per line (two replicas max)
        int a1,a2,a3,a4,a5,a6,a7,a9,a10,a11,a12,a13,a14,a15;
        double a0,a8,a16,a17;
        unsigned long long int a18;
        
        while(infile >> a0 >> a1 >> a2 >> a3 >> a4 >> a5 >> a6 >> a7 >> a8 >> a9 >> a10 >> a11 >> a12 >> a13 >> a14 >> a15 >>
            a16 >> a17 >> a18){
            iteration_idx = a18;
            break;
        }
    }
    
    if (num_replicas==1){
        // Assume 16 elements per line (two replicas max)
        int a1,a2,a3,a4,a5,a6,a7;
        double a0,a8,a9;
        unsigned long long int a10;
        
        while(infile >> a0 >> a1 >> a2 >> a3 >> a4 >> a5 >> a6 >> a7 >>
              a8 >> a9 >> a10){
            iteration_idx = a10;
            break;
        }
    }
    
    infile.close();
    
    return iteration_idx;
}

/*--------------------------------------------------------------------*/

double get_eta(int D, int L, int N, int l_A,
              double U, double t, double beta,
              int bin_size, int bins_wanted,
              int seed, string subgeometry,
              int num_replicas){
        
    double eta;
    string state_name;
    
    // Name of system state file
    state_name=to_string(D)+"D_"+to_string(L)+
    "_"+to_string(N)+"_"+to_string(l_A)+"_"+
    to_string(U)+"_"+to_string(t)+"_"+
    to_string(beta)+"_"+to_string(bin_size)+"_"+
    "system-state_"+to_string(seed)+"_"+subgeometry+"_"+
    to_string(num_replicas)+".dat";
    
    // NOTE: For consistency, may rewrite function to get the number
    // of kinks from the path structure created with load_paths()
    std::ifstream infile(state_name);
    
    eta = -1;
    if (num_replicas==2){
        // Assume 16 elements per line (two replicas max)
        int a1,a2,a3,a4,a5,a6,a7,a9,a10,a11,a12,a13,a14,a15;
        double a0,a8,a16,a17;
        unsigned long long int a18;
        
        while(infile >> a0 >> a1 >> a2 >> a3 >> a4 >> a5 >> a6 >> a7 >> a8 >> a9 >> a10 >> a11 >> a12 >> a13 >> a14 >> a15 >>
            a16 >> a17 >> a18){
            eta = a17;
            break;
        }
    }
    
    if (num_replicas==1){
        // Assume 16 elements per line (two replicas max)
        int a1,a2,a3,a4,a5,a6,a7;
        double a0,a8,a9;
        unsigned long long int a10;
        
        while(infile >> a0 >> a1 >> a2 >> a3 >> a4 >> a5 >> a6 >> a7 >>
              a8 >> a9 >> a10){
            eta = a9;
            break;
        }
    }
    
    infile.close();
    
    return eta;
}

/*--------------------------------------------------------------------*/

vector<double> get_N_tracker(vector<vector<Kink> > paths,
                             int num_replicas,int M,double beta){

    vector<double> N_tracker (num_replicas,0.0);
    double l_path,dN;
    int current,next;

    for (int r=0; r<num_replicas; r++){
        for (int site=0; site<M; site++){
            current=site;
            next=paths[r][current].next;
            while (next!=-1){
                l_path = paths[r][next].tau - paths[r][current].tau;
                dN = paths[r][current].n * l_path/beta;

                N_tracker[r] += dN;

                current = next;
                next = paths[r][next].next;
            }
            l_path = beta - paths[r][current].tau;
            dN = paths[r][current].n * l_path/beta;

            N_tracker[r] += dN;
        }
        N_tracker[r] = round(N_tracker[r]);
                
    }
    return N_tracker;
}

/*--------------------------------------------------------------------*/

vector<int> get_head_idx(vector<vector<Kink> > paths,
                             int num_replicas,int M){

    vector<int> head_idx (num_replicas,-1);
    int current,next;
    
    
    for (int r=0; r<num_replicas; r++){
        for (int site=0; site<M; site++){
            current=site;
            next=paths[r][current].next;
            while (next!=-1){

                if ((paths[r][next].src==paths[r][next].dest)
                    && (paths[r][next].n-paths[r][current].n==-1
                    && paths[r][next].src_replica==
                       paths[r][next].dest_replica)){
                    head_idx[r] = next;
                }
                current = next;
                next = paths[r][next].next;
            }
        }
    }
    return head_idx;
}

/*--------------------------------------------------------------------*/

vector<int> get_tail_idx(vector<vector<Kink> > paths,
                             int num_replicas,int M){

    vector<int> tail_idx (num_replicas,-1);
    int current,next;
    
    for (int r=0; r<num_replicas; r++){
        for (int site=0; site<M; site++){
            current=site;
            next=paths[r][current].next;
            while (next!=-1){
                if ((paths[r][next].src==paths[r][next].dest)
                    && (paths[r][next].n-paths[r][current].n==1
                    && paths[r][next].src_replica==
                       paths[r][next].dest_replica)){
                    tail_idx[r] = next;
                }
                current = next;
                next = paths[r][next].next;
            }
        }
    }
    return tail_idx;
}

/*--------------------------------------------------------------------*/

vector<int> get_N_zero(vector<vector<Kink> > paths,
                             int num_replicas,int M){

    vector<int> N_zero(num_replicas,0);
    
    for (int r=0; r<num_replicas; r++){
        for (int site=0; site<M; site++){
            N_zero[r] += paths[r][site].n;
        }
    }
    return N_zero;
}

/*--------------------------------------------------------------------*/

vector<int> get_N_beta(vector<vector<Kink> > paths,
                             int num_replicas,int M){

    vector<int> N_beta(num_replicas,0);
    int current,next;
    
    for (int r=0; r<num_replicas; r++){
        for (int site=0; site<M; site++){
            current=site;
            next=paths[r][current].next;
            while (next!=-1){
                current = next;
                next = paths[r][next].next;
            }
            N_beta[r] += paths[r][current].n;
        }
    }
    return N_beta;
}

/*--------------------------------------------------------------------*/

vector<vector<int> > get_last_kinks(vector<vector<Kink> > paths,
                             int num_replicas,int M){

    vector<vector<int> > last_kinks(num_replicas,vector<int> (M,-1));
    int current,next;
    
    for (int r=0; r<num_replicas; r++){
        for (int site=0; site<M; site++){
            current=site;
            next=paths[r][current].next;
            while (next!=-1){
                current = next;
                next = paths[r][next].next;
            }
            last_kinks[r][site] = current;
        }
    }
    return last_kinks;
}


/*--------------------------------------------------------------------*/

int get_num_swaps(vector<vector<Kink> > paths,
                             int num_replicas,int M){

    int num_swaps,current,next;
    
    num_swaps = 0;
    
    for (int site=0; site<M; site++){
        current=site;
        next=paths[0][current].next;
        while (next!=-1){

            if (paths[0][next].src_replica!=
                paths[0][next].dest_replica){num_swaps+=1;}
            current = next;
            next = paths[0][next].next;
        }
    }
    
    return num_swaps;
}

/*--------------------------------------------------------------------*/

// Create function that calculates vector norm given array
double norm(vector<double> point){
    
    double squared_sum=0;
        
    for (size_t i=0; i<point.size(); i++){
        squared_sum += point[i]*point[i];
    }
    
    return sqrt(squared_sum);
}

/*--------------------------------------------------------------------*/

void build_hypercube_adjacency_matrix(int L,int D, string boundary_condition,
                                      vector<vector<int> > &adjacency_matrix){
    
    int M = pow(L,D);
    int site_left,site_right,site_up,site_down,site_high,site_low;
    int top_row_ctr,bottom_row_ctr;
    
    vector<int> rows (D*2,0);
    
    top_row_ctr=L;
    bottom_row_ctr=0;
    
    // Initialize adjacency matrix with placeholder zeroes
    for (int i=0; i<M; i++){adjacency_matrix.push_back(rows);}
    
    for (int site=0; site<M; site++){
        
        // Left neighbor
        if (site%L==0){site_left=site+(L-1);}
        else {site_left=site-1;}
        
        // Right neighbor
        if ((site+1)%L==0){site_right=site-(L-1);}
        else {site_right=site+1;}
        
        // Top neighbor
        if (site%(L*L)==0){top_row_ctr=L;}
        if (top_row_ctr>0){
            site_up = site+L*(L-1);
            top_row_ctr-=1;
        }
        else {site_up=site-L;}

        // Bottom neighbor
        if ((site-L*(L-1))%(L*L)==0){bottom_row_ctr=L;}
        if (bottom_row_ctr>0){
            site_down = site-L*(L-1);
            bottom_row_ctr-=1;
        }
        else {site_down=site+L;}
        
        // High neighbor (i.e, layer above in 3D; ran out of words, sorry)
        if (site-L*L<0){site_high=site+L*L*(L-1);}
        else {site_high=site-L*L;}
        
        // Low neighbor (i.e, layer below in 3D; ran out of words, sorry)
        if (site+L*L>L*L*L-1){site_low=site-L*L*(L-1);}
        else {site_low=site+L*L;}
        
        // Populate the adjacency matrix (in compact form)
        if (L<2){
            cout << "ERROR: CODE DOES NOT SUPPORT L<2 AT THE MOMENT" << endl;
        }
        else{
            adjacency_matrix[site][0]=site_left;
            adjacency_matrix[site][1]=site_right;
            if (D>1){
                adjacency_matrix[site][2]=site_up;
                adjacency_matrix[site][3]=site_down;
            }
            if (D>2){
                adjacency_matrix[site][4]=site_high;
                adjacency_matrix[site][5]=site_low;
            }
        }
    }
    return;
}
    
/*--------------------------------------------------------------------*/

void build_adjacency_matrix(int L,int D,string boundary_condition,
                            vector<vector<bool> >&adjacency_matrix){

    // Variable declarations
    int M = pow(L,D); // Number of lattice points
    int ctr,a1,a2,a3;
    double r_NN;
    vector<double> points_difference (D,0.0);

    // Initialize normalized basis vectors
    vector<double> a1_vec {1.0,0.0,0.0};
    vector<double> a2_vec {0.0,1.0,0.0};
    vector<double> a3_vec {0.0,0.0,1.0};
    
    // Initialize array that will store all the points
    vector<double> empty_point (D,0);
    vector<vector<double> > points (M,empty_point);
    
    // Norm of basis vectors
    a1 = norm(a1_vec);
    a2 = norm(a2_vec);
    a3 = norm(a3_vec);
        
    // Build the lattice vectors
    ctr = 0;
    if (D==1){
        for (int i1=0; i1<L; i1++){
            points[ctr][0] = i1*a1;
            ctr++;
        }
    }
    else if (D==2){
        for (int i1=0; i1<L; i1++){
            for (int i2=0; i2<L; i2++){
                points[ctr][0] = i1*a1;
                points[ctr][1] = i2*a2;
                ctr++;
            }
        }
    }
    else{ // D==3
        for (int i1=0; i1<L; i1++){
            for (int i2=0; i2<L; i2++){
                for (int i3=0; i3<L; i3++){
                    points[ctr][0] = i1*a1;
                    points[ctr][1] = i2*a2;
                    points[ctr][2] = i3*a3;
                    ctr++;
                }
            }
        }
    }
    
    // Set the nearest-neighbor distance
    r_NN = a1;
    
    // Set adjacency matrix elements by comparing inter-site distances
    for (int i=0; i<M; i++){
        for (int j=i+1; j<M; j++){
            for (int axis=0; axis<D; axis++){
                points_difference[axis] = points[i][axis]-points[j][axis];
            }
            if (boundary_condition=="pbc"){
                adjacency_matrix[i][j] = norm(points_difference) <= r_NN
                || norm(points_difference) == L-1;
            }
            else if (boundary_condition=="obc"){
                adjacency_matrix[i][j] = norm(points_difference) <= r_NN;
            }
            else{
                cout << "ERROR: Boundary condition not supported" << endl;
            }
            // Fill out the symmetric elements
            adjacency_matrix[j][i] = adjacency_matrix[i][j];
        }
    }
    
    return;
 }

/*--------------------------------------------------------------------*/

void create_sub_sites(vector<int> &sub_sites,int l_max,int L,int D,int M,
                      string geometry){
    // Hard coded to cluster of sites in 1D and to SQUARE region in 2D, for now.
    // Doesn't work for 3D yet.
    
    int ctr,next_sub_site,horizontal_direction,vertical_direction,
    horizontal_direction_old,y;
    size_t m_max;
    
    if (D==1 || L==2){ // cluster
        for (int i=0; i<l_max; i++){sub_sites.push_back(i);}
    }
    else if (D==2){
        
        if (geometry=="square"){
            m_max = pow(l_max,D); // maximum number of total subsystem sites
            ctr=0; // unused at the moment
            y=0;
            for (int l=0; l<l_max; l++){
                next_sub_site = l;
                sub_sites.push_back(next_sub_site);
                
                for (int j=1; j<=l; j++){
                    next_sub_site = l+j*L;
                    sub_sites.push_back(next_sub_site);
                    if (j==l){
                        for (int i=1; i<=l; i++){
                            next_sub_site -= 1;
                            sub_sites.push_back(next_sub_site);
                        }
                    }
                }
                y += L;
                ctr+=1; // unused at the moment
            }
        }
        
        else if (geometry=="strip"){
            m_max = L*l_max;
            next_sub_site = -1;
            horizontal_direction = +1;
            horizontal_direction_old = +1;
            vertical_direction = 0;
            ctr = 0;
            while (sub_sites.size()!=m_max){
                if (ctr==L){
                    vertical_direction = +L;
                    horizontal_direction = 0;
                    ctr=0;
                }
                else if (sub_sites.size()>2 && ctr==1){
                    vertical_direction = 0;
                    horizontal_direction = (-1)*horizontal_direction_old;
                    horizontal_direction_old = horizontal_direction;
                }
                else {
                    // nothing
                }
                next_sub_site += (horizontal_direction+vertical_direction);
                sub_sites.push_back(next_sub_site);
                ctr++;
            }
        }
        
        
    }
    
//    if (L>=2 && m_A>M){cout<<"ERROR: l_A needs to be smaller than L"<<endl; exit(1);}
//
//    if (D==1 || L==2){ // cluster
//        for (int i=0; i<l_A; i++){sub_sites.push_back(i);}
//    }
//    else if (D==2){

    
//    }
//    else if (D==3){
//        cout << "ERROR: create_sub_sites does not support 3D yet" << endl;
//        exit(1);
//    }
//    else{
//        // lol
//        cout << "ERROR: Please choose either  D=1,D=2, or D=3" << endl;
//        exit(1);
//    }
//    cout << endl;
    
    return;
}
/*--------------------------------------------------------------------*/

void insert_worm(vector<Kink> &paths, int &num_kinks, int &head_idx,
                 int &tail_idx, int M, int N, double U, double mu, double t,
                 double beta, double eta, bool canonical, double &N_tracker,
                 int &N_zero, int &N_beta, vector<int> &last_kinks,
                 unsigned long long int &insert_worm_attempts,
                 unsigned long long int &insert_worm_accepts,
                 unsigned long long int &insert_anti_attempts,
                 unsigned long long int &insert_anti_accepts,
                 RNG &rng){
    
    // Variable declarations
    int k,n,src,next,n_head,n_tail,src_replica;
    double tau,tau_h,tau_t,tau_prev,tau_next,tau_flat,l_path,dN,dV,p_iw,p_dw,R;
    bool is_worm;
    
    // Can only perform update if there are no worm ends
    if (head_idx != -1 || tail_idx != -1){return;}
        
    // Randomly sample a flat interval (or kink if you like)
    //boost::random::uniform_int_distribution<> flats(0, num_kinks-1);
    k = rng.randInt(num_kinks-1);
    
    // Extract the attributes of the kink at the bottom of the flat interval
    tau = paths[k].tau;
    n = paths[k].n;
    src = paths[k].src;
    next = paths[k].next;
    src_replica = paths[k].src_replica; //due to way replica indices are coded
    
    // Calculate the length of the flat interval
    tau_prev = tau;
    if (next != -1) // tau_next extractable iff sampled kink is not the last
        tau_next = paths[next].tau;
    else
        tau_next = beta;
    tau_flat = tau_next - tau_prev;
    
    // Randomly choose where to insert worm ends in the flat interval
    //boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    tau_h = tau_prev + tau_flat*rng.rand();
    tau_t = tau_prev + tau_flat*rng.rand();
    
    // Based on worm end time, determine worm type: antiworm or worm
    if (tau_h > tau_t){
        is_worm = true;
        insert_worm_attempts += 1; // Attempts counter
    }
    else{
        is_worm = false;
        insert_anti_attempts += 1;
    }
    
    // Determine the no. of particles after each worm end
    if (is_worm){
        n_tail = n + 1;
        n_head = n;
    }
    else{
        n_tail = n;
        n_head = n - 1;
    }
    
    // Reject update if illegal worm insertion is proposed
    if (n == 0 && !(is_worm)){insert_anti_attempts-=1;return;}
    if (tau_h == tau_prev || tau_t == tau_prev){return;}
    if (tau_h == tau_t){return;}
    
    // Determine length of modified path and particle change
    l_path = tau_h - tau_t;
    dN = l_path/beta;
    
    // Canonical simulations: Restrict updates to interval N:(N-1,N+1)
    if (canonical)
        if ((N_tracker+dN) < (N-1) || (N_tracker+dN) > (N+1)){return;}
    
    // Calculate the difference in diagonal energy dV = \epsilon_w - \epsilon
    dV = (U/2.0)*(n_tail*(n_tail-1)-n_head*(n_head-1)) - mu*(n_tail-n_head);
    
    // Build the Metropolis ratio (R)
    p_dw = 0.5;
    p_iw = 0.5;
    R = eta * eta * n_tail * exp(-dV*(tau_h-tau_t))* (p_dw/p_iw) *
    num_kinks * tau_flat * tau_flat;

    // Metropolis sampling
    if (rng.rand() < R){ // Accept

        // Activate the first two available kinks
        if (is_worm){
            paths[num_kinks]=Kink(tau_t,n_tail,src,src,k,num_kinks+1,
                                  src_replica,src_replica);
            paths[num_kinks+1]=Kink(tau_h,n_head,src,src,num_kinks,next,
                                    src_replica,src_replica);
            
            // Save indices of head & tail kinks
            head_idx = num_kinks+1;
            tail_idx = num_kinks;
            
            // Add to Acceptance counter
            insert_worm_accepts += 1;
        }
        else{ // Antiworm
            paths[num_kinks]=Kink(tau_h,n_head,src,src,k,num_kinks+1,
                                  src_replica,src_replica);
            paths[num_kinks+1]=Kink(tau_t,n_tail,src,src,num_kinks,next,
                                    src_replica,src_replica);
            
            // Save indices of head & tail kinks
            head_idx = num_kinks;
            tail_idx = num_kinks+1;
            
            // Add to Acceptance counter
            insert_anti_accepts += 1;
        }
        
        // "Connect" next of lower bound kink to nearest worm end
        paths[k].next = num_kinks;
        
        // "Connect" prev of next kink to nearest worm end
        if(next!=-1){paths[next].prev = num_kinks+1;}
        
        // Update trackers for: no of active kinks, total particles
        num_kinks += 2;
        N_tracker += dN;
        
        // If later worm end is last kink on site, update last kinks tracker vec
        if (next==-1){
            if (is_worm){last_kinks[src]=head_idx;}
            else {last_kinks[src]=tail_idx;}
        }
        
        return;
    }
    else // Reject
        return;
}

/*--------------------------------------------------------------------*/

void insert_worm_2(vector<Kink> &paths, int &num_kinks, int &head_idx,
                 int &tail_idx, int M, int N, double U, double mu, double t,
                 double beta, double eta, bool canonical, double &N_tracker,
                 int &N_zero, int &N_beta, vector<int> &last_kinks,
                 unsigned long long int &insert_worm_attempts,
                 unsigned long long int &insert_worm_accepts,
                 unsigned long long int &insert_anti_attempts,
                 unsigned long long int &insert_anti_accepts,
                 RNG &rng){
    
    // Variable declarations
    int k,n,src,next,n_head,n_tail,src_replica;
    double tau,tau_h,tau_t,tau_prev,tau_next,l_path,dN,dV,p_iw,p_dw,R;
    bool is_worm;
    double inv_e = boost::math::constants::exp_minus_one<double>(); // 1/e
    
    // Can only perform update if there are no worm ends
    if (head_idx != -1 || tail_idx != -1){return;}
        
    // Randomly sample a flat interval (or kink if you like)
    //boost::random::uniform_int_distribution<> flats(0, num_kinks-1);
    k = rng.randInt(num_kinks-1);
    
    // Extract the attributes of the kink at the bottom of the flat interval
    tau = paths[k].tau;
    n = paths[k].n;
    src = paths[k].src;
    next = paths[k].next;
    src_replica = paths[k].src_replica; //due to way replica indices are coded
   
    // Calculate the length of the flat interval
    tau_prev = tau;
    if (next != -1) // tau_next extractable iff sampled kink is not the last
        tau_next = paths[next].tau;
    else
        tau_next = beta;

    // Randomly choose to insert worm or antiworm
    if (rng.rand() < 0.5){
        is_worm = true;
        insert_worm_attempts += 1;
    }
    else{
        is_worm = false;
        insert_anti_attempts += 1;
    }
    // p_type = 1/2;
    
    // Determine the no. of particles after each worm end
    if (is_worm){
        n_tail = n + 1;
        n_head = n;
    }
    else{
        n_tail = n;
        n_head = n - 1;
    }
    // Reject antiworm insertion if no particles in the flat region
    if (n == 0 && !(is_worm)){insert_anti_attempts-=1;return;}

    // Calculate the difference in diagonal energy dV = \epsilon_w - \epsilon
    dV = (U/2.0)*(n_tail*(n_tail-1)-n_head*(n_head-1)) - mu*(n_tail-n_head);
    // To make acceptance ratio unity,antiworm needs to sample w/ dV=eps-eps_w
    if (!is_worm){dV *= -1;} // dV=eps-eps_w

    // debugging
    // if (!is_worm){return;}

    /* ==================== Truncated Exponential Sampling ================== */

    // sample tau_1
    double y,x,a,b,c,tau_1,tau_2,a_new;
    // long double arg;
    double arg;

    c = -dV;
    b = tau_next;
    a = tau_prev;

    x = rng.rand();

    // Compute normalization of truncated exponential dist.
    double Z_1, Z_2, Z;
    // if (dV == 0){cout << "NOOOOOO" << endl;}

    Z = (exp(c*(b-a)) + a*c - b*c - 1)/(c*c); // joint

    Z_1 = (1/c) * (exp(c*(b-a)) - 1) - (b - a); // marginalized dist.
    
    //
    y = Z_1*x - (1/c)*exp(c*(b-a)) - a;

    // this was the sampling without simplifying Ace^(cy) and no restriction
    // A = -(1/c)*exp(c*b);
    // arg = max(-1/exp(1), A*c*exp(c*y));
    // cout << c*(b+y) << " " << -exp(c*(b+y)) << endl;
    // if (abs(arg) < 1e-15){arg = 0.0;}

    // Determine LambertW branch & compute tau
    if (abs(c*(b+y)) < 600){
    arg = max(-inv_e, -exp(c*(b+y)));  // z in boost documentation

    }
    else{ // exponent c*y is too large. Better to ignore and carry on. 
        return;
    }


    if (c < 0){ // k = 0 branch
        tau_1 = (1/c)*lambert_w0(arg)-y;
    }
    else {      // k = -1 branch
        tau_1 = (1/c)*lambert_wm1(arg)-y;
    }

    // sample tau_2
    x = rng.rand();
    a_new = tau_1; // this is the new lower bound for the simple truncexpon sampling
    Z_2 = 1.0 - exp(-(-c)*(b-a_new)); // conditional dist.
    tau_2 = a_new - log(1.0-Z_2*x)  / (-c);

    if (is_worm){
        tau_t = tau_1;
        tau_h = tau_2;
    }
    else{
        tau_h = tau_1;
        tau_t = tau_2;
    }

    // cout << a << " "<< tau_1 << " "<< a_new << " "<< tau_2 << " "<< b << " " 
    // << c*y << " " << c+y << " " << arg << endl ;
    /* ====================================================================== */ 

    // Reject update if illegal worm insertion is proposed
    if (tau_h == tau_prev || tau_t == tau_prev){return;}
    if (tau_h == tau_t){return;}

    // Determine length of modified path and particle change
    l_path = tau_h - tau_t;
    dN = l_path/beta;
    
    // Canonical simulations: Restrict updates to interval N:(N-1,N+1)
    if (canonical)
        if ((N_tracker+dN) < (N-1) || (N_tracker+dN) > (N+1)){return;}
    
    // Build the Metropolis ratio (R)
    p_dw = 0.5;
    p_iw = 0.5;
    R = eta * eta * n_tail * Z * num_kinks * (p_dw/p_iw) * 2;
    // cout << R << endl;

    // Metropolis sampling
    if (rng.rand() < R){ // Accept

        // Activate the first two available kinks
        if (is_worm){
            paths[num_kinks]=Kink(tau_t,n_tail,src,src,k,num_kinks+1,
                                  src_replica,src_replica);
            paths[num_kinks+1]=Kink(tau_h,n_head,src,src,num_kinks,next,
                                    src_replica,src_replica);
            
            // Save indices of head & tail kinks
            head_idx = num_kinks+1;
            tail_idx = num_kinks;
            
            // Add to Acceptance counter
            insert_worm_accepts += 1;
        }
        else{ // Antiworm
            paths[num_kinks]=Kink(tau_h,n_head,src,src,k,num_kinks+1,
                                  src_replica,src_replica);
            paths[num_kinks+1]=Kink(tau_t,n_tail,src,src,num_kinks,next,
                                    src_replica,src_replica);
            
            // Save indices of head & tail kinks
            head_idx = num_kinks;
            tail_idx = num_kinks+1;
            
            // Add to Acceptance counter
            insert_anti_accepts += 1;
        }
        
        // "Connect" next of lower bound kink to nearest worm end
        paths[k].next = num_kinks;
        
        // "Connect" prev of next kink to nearest worm end
        if(next!=-1){paths[next].prev = num_kinks+1;}
        
        // Update trackers for: no of active kinks, total particles
        num_kinks += 2;
        N_tracker += dN;
        
        // If later worm end is last kink on site, update last kinks tracker vec
        if (next==-1){
            if (is_worm){last_kinks[src]=head_idx;}
            else {last_kinks[src]=tail_idx;}
        }
        
        return;
    }
    else // Reject
        return;
}

/*--------------------------------------------------------------------*/

void delete_worm(vector<Kink> &paths, int &num_kinks, int &head_idx,
                 int &tail_idx, int M, int N, double U, double mu, double t,
                 double beta, double eta, bool canonical, double &N_tracker,
                 int &N_zero, int &N_beta, vector<int> &last_kinks,
                 unsigned long long int &delete_worm_attempts,
                 unsigned long long int &delete_worm_accepts,
                 unsigned long long int &delete_anti_attempts,
                 unsigned long long int &delete_anti_accepts,
                 RNG &rng){
    
    // Variable declarations
    int src,prev,next,n_head,n_tail;
    int prev_h,next_h,prev_t,next_t,high_end,low_end;
    double tau_h,tau_t,tau_prev,tau_next,tau_flat,l_path,dN,dV,p_iw,p_dw,R;
    bool is_worm;
    
    // Can only propose worm deletion if both worm ends are present
    if (head_idx == -1 || tail_idx == -1){return;}
    
    // Can only delete worm if wormends are on same flat interval
    if (head_idx != paths[tail_idx].prev &&
        tail_idx != paths[head_idx].prev)
        return;
    
    // Extract worm end attributes
    tau_h = paths[head_idx].tau; // Head attributes
    n_head = paths[head_idx].n;
    src = paths[head_idx].src;
    prev_h = paths[head_idx].prev;
    next_h = paths[head_idx].next;
    
    tau_t = paths[tail_idx].tau; // Tail attributes
    n_tail = paths[tail_idx].n;
    src = paths[tail_idx].src;
    prev_t = paths[tail_idx].prev;
    next_t = paths[tail_idx].next;

    // Identify the type of worm
    if (tau_h > tau_t)
        is_worm = true;
    else
        is_worm = false; // antiworm
    
    // Identify lower and upper bound of flat interval where worm lives
    if(is_worm){
        tau_prev = paths[prev_t].tau;
        delete_worm_attempts += 1; // Attempts counter
        if(paths[head_idx].next == -1)
            tau_next = beta;
        else
            tau_next = paths[next_h].tau;
            }
    else{ // antiworm
        tau_prev = paths[prev_h].tau;
        delete_anti_attempts += 1; // Attempts counter
        if (paths[tail_idx].next == -1)
            tau_next = beta;
        else
            tau_next = paths[next_t].tau;
            }
    
    // Calculate the length of the flat interval
    tau_flat = tau_next - tau_prev;
    
    // Define upper,lower index variables independent of worm type
    if (is_worm){
        next = next_h;
        prev = prev_t;
        high_end = head_idx;
        low_end = tail_idx;
    }
    else{
        next = next_t;
        prev = prev_h;
        high_end = tail_idx;
        low_end = head_idx;
    }
    
    // Determine length of modified path and particle change
    l_path = tau_h-tau_t;
    dN = -l_path/beta;
    
    // Canonical simulations: Restrict updates to interval N:(N-1,N+1)
    if (canonical)
        if ((N_tracker+dN) < (N-1) || (N_tracker+dN) > (N+1)){return;}
        
    // Calculate the difference in diagonal energy dV = \epsilon_w - \epsilon
    dV = (U/2.0)*(n_tail*(n_tail-1)-n_head*(n_head-1)) - mu*(n_tail-n_head);
    
    // Build the Metropolis ratio (R)
    p_dw = 1.0;
    p_iw = 1.0;
    R = (eta*eta) * n_tail * exp(-dV*(tau_h-tau_t))* (p_dw/p_iw) *
    (num_kinks-2) * (tau_flat*tau_flat);
    R = 1.0/R;
    
    // Metropolis sampling
    //boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    if (rng.rand() < R){ // Accept
        
        // Add to Acceptance counter
        if (is_worm)
            delete_worm_accepts += 1;
        else
            delete_anti_accepts += 1;
        
        // Stage 1: Delete the higher worm end
        
        // num_kinks-1 will be swapped. Modify links to these.
        if (paths[num_kinks-1].next!=-1)
            paths[paths[num_kinks-1].next].prev = high_end;
        paths[paths[num_kinks-1].prev].next = high_end;
                
        swap(paths[high_end],paths[num_kinks-1]);
        
        // Upper,lower bounds of flat could've been swapped. Correct if so.
        if (prev==num_kinks-1){prev=high_end;}
        else if (next==num_kinks-1){next=high_end;}
        else if (low_end==num_kinks-1){low_end=high_end;}
        else {;}
        
        // The swapped kink could've been the last on its site.
        if (paths[high_end].next==-1){
            last_kinks[paths[high_end].src]=high_end;
        }
        
        // Connect upper,lower bounds to lower worm end
        if (next!=-1)
            paths[next].prev = low_end;
        paths[low_end].next = next;
        
        if (next==-1){last_kinks[src]=low_end;}

        // Stage 2: Delete the lower worm end
        
        // num_kinks-2 will be swapped. Modify links to these.
        if (paths[num_kinks-2].next!=-1)
            paths[paths[num_kinks-2].next].prev = low_end;
        paths[paths[num_kinks-2].prev].next = low_end;
        
        swap(paths[low_end],paths[num_kinks-2]);

        if (prev==num_kinks-2){prev=low_end;}
        else if (next==num_kinks-2){next=low_end;}
        else {;}

        if (paths[low_end].next==-1){
            last_kinks[paths[low_end].src]=low_end;
        }
        
        if (next!=-1)
            paths[next].prev = prev;
        paths[prev].next = next;
        
        if (next==-1){last_kinks[src]=prev;}

        // Deactivate the head,tail indices
        head_idx = -1;
        tail_idx = -1;
        
        // Update trackers for: num of active kinks, total particles
        num_kinks -= 2;
        N_tracker += dN;
        
        return;
    }
        
    else // Reject
            return;
}

/*--------------------------------------------------------------------*/

void delete_worm_2(vector<Kink> &paths, int &num_kinks, int &head_idx,
                 int &tail_idx, int M, int N, double U, double mu, double t,
                 double beta, double eta, bool canonical, double &N_tracker,
                 int &N_zero, int &N_beta, vector<int> &last_kinks,
                 unsigned long long int &delete_worm_attempts,
                 unsigned long long int &delete_worm_accepts,
                 unsigned long long int &delete_anti_attempts,
                 unsigned long long int &delete_anti_accepts,
                 RNG &rng){
    
    // Variable declarations
    int src,prev,next,n_head,n_tail;
    int prev_h,next_h,prev_t,next_t,high_end,low_end;
    double tau_h,tau_t,tau_prev,tau_next,l_path,dN,dV,p_iw,p_dw,R;
    bool is_worm;
    
    // Can only propose worm deletion if both worm ends are present
    if (head_idx == -1 || tail_idx == -1){return;}
    
    // Can only delete worm if wormends are on same flat interval
    if (head_idx != paths[tail_idx].prev &&
        tail_idx != paths[head_idx].prev)
        return;
    
    // Extract worm end attributes
    tau_h = paths[head_idx].tau; // Head attributes
    n_head = paths[head_idx].n;
    src = paths[head_idx].src;
    prev_h = paths[head_idx].prev;
    next_h = paths[head_idx].next;
    
    tau_t = paths[tail_idx].tau; // Tail attributes
    n_tail = paths[tail_idx].n;
    src = paths[tail_idx].src;
    prev_t = paths[tail_idx].prev;
    next_t = paths[tail_idx].next;

    // Identify the type of worm
    if (tau_h > tau_t)
        is_worm = true;
    else
        is_worm = false; // antiworm
    
    // Identify lower and upper bound of flat interval where worm lives
    if(is_worm){
        tau_prev = paths[prev_t].tau;
        delete_worm_attempts += 1; // Attempts counter
        if(paths[head_idx].next == -1)
            tau_next = beta;
        else
            tau_next = paths[next_h].tau;
            }
    else{ // antiworm
        tau_prev = paths[prev_h].tau;
        delete_anti_attempts += 1; // Attempts counter
        if (paths[tail_idx].next == -1)
            tau_next = beta;
        else
            tau_next = paths[next_t].tau;
            }
     
    // Define upper,lower index variables independent of worm type
    if (is_worm){
        next = next_h;
        prev = prev_t;
        high_end = head_idx;
        low_end = tail_idx;
    }
    else{
        next = next_t;
        prev = prev_h;
        high_end = tail_idx;
        low_end = head_idx;
    }
    
    // Determine length of modified path and particle change
    l_path = tau_h-tau_t;
    dN = -l_path/beta;
    
    // Canonical simulations: Restrict updates to interval N:(N-1,N+1)
    if (canonical)
        if ((N_tracker+dN) < (N-1) || (N_tracker+dN) > (N+1)){return;}
        
    // Calculate the difference in diagonal energy dV = \epsilon_w - \epsilon
    dV = (U/2.0)*(n_tail*(n_tail-1)-n_head*(n_head-1)) - mu*(n_tail-n_head);

    // In inverse update, antiworm could've been sampled w/ dV=\eps-\eps_w
    if (!is_worm){dV *= -1;} // dV=eps-eps_w

    // Compute normalization constant of joint distribution
    double Z,a,b,c;
    a = tau_prev;
    b = tau_next;
    c = -dV;
    Z = (exp(c*(b-a)) + a*c - b*c - 1)/(c*c);

    // In inverse update (insert), worm type must be randomly chosen
    // double p_type = 1/2;
    
    // Build the Metropolis ratio (R)
    p_dw = 1.0;
    p_iw = 1.0;
    R = eta * eta * n_tail * Z * (num_kinks-2) * (p_dw/p_iw) * 2;
    R = 1.0/R;
    
    // Metropolis sampling
    //boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    if (rng.rand() < R){ // Accept
        
        // Add to Acceptance counter
        if (is_worm)
            delete_worm_accepts += 1;
        else
            delete_anti_accepts += 1;
        
        // Stage 1: Delete the higher worm end
        
        // num_kinks-1 will be swapped. Modify links to these.
        if (paths[num_kinks-1].next!=-1)
            paths[paths[num_kinks-1].next].prev = high_end;
        paths[paths[num_kinks-1].prev].next = high_end;
                
        swap(paths[high_end],paths[num_kinks-1]);
        
        // Upper,lower bounds of flat could've been swapped. Correct if so.
        if (prev==num_kinks-1){prev=high_end;}
        else if (next==num_kinks-1){next=high_end;}
        else if (low_end==num_kinks-1){low_end=high_end;}
        else {;}
        
        // The swapped kink could've been the last on its site.
        if (paths[high_end].next==-1){
            last_kinks[paths[high_end].src]=high_end;
        }
        
        // Connect upper,lower bounds to lower worm end
        if (next!=-1)
            paths[next].prev = low_end;
        paths[low_end].next = next;
        
        if (next==-1){last_kinks[src]=low_end;}

        // Stage 2: Delete the lower worm end
        
        // num_kinks-2 will be swapped. Modify links to these.
        if (paths[num_kinks-2].next!=-1)
            paths[paths[num_kinks-2].next].prev = low_end;
        paths[paths[num_kinks-2].prev].next = low_end;
        
        swap(paths[low_end],paths[num_kinks-2]);

        if (prev==num_kinks-2){prev=low_end;}
        else if (next==num_kinks-2){next=low_end;}
        else {;}

        if (paths[low_end].next==-1){
            last_kinks[paths[low_end].src]=low_end;
        }
        
        if (next!=-1)
            paths[next].prev = prev;
        paths[prev].next = next;
        
        if (next==-1){last_kinks[src]=prev;}

        // Deactivate the head,tail indices
        head_idx = -1;
        tail_idx = -1;
        
        // Update trackers for: num of active kinks, total particles
        num_kinks -= 2;
        N_tracker += dN;
        
        return;
    }
        
    else // Reject
            return;
}

/*--------------------------------------------------------------------*/

void insertZero(vector<Kink> &paths, int &num_kinks, int &head_idx,
                int &tail_idx, int M, int N, double U, double mu, double t,
                double beta, double eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                unsigned long long int &insertZero_worm_attempts,
                unsigned long long int &insertZero_worm_accepts,
                unsigned long long int &insertZero_anti_attempts,
                unsigned long long int &insertZero_anti_accepts,
                RNG &rng){
    
    // Variable declarations
    int n,src,next,n_head,n_tail,i,dest_replica;
    double tau_flat,l_path,dN,dV,R,p_type,tau_new,p_wormend,C,W,
    p_dz,p_iz;
    bool is_worm;

    // Cannot insert if there's two worm ends present
    if (head_idx != -1 and tail_idx != -1){return;}
        
    // Randomly select site on which to insert worm/antiworm from tau=0
    //boost::random::uniform_int_distribution<> sites(0, M-1);
    i = rng.randInt(M-1);
    
    // Extract attributes of insertion flat
    n = paths[i].n;
    src = paths[i].src;
    next = paths[i].next;
    dest_replica = paths[i].dest_replica;
    
    // Determine the length of insertion flat interval
    if (next != -1)
        tau_flat = paths[next].tau;
    else
        tau_flat = beta;
    
    // Choose worm/antiworm insertion based on worm ends present
    //boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    if (head_idx==-1 and tail_idx==-1){ // no worm ends present
        if (n==0){ // can only insert worm, not antiworm
            is_worm = true;
            p_type = 1.0;
        }
        else{ // choose worm or antiworm insertion with equal probability
            if (rng.rand() < 0.5)
                is_worm = true;
            else
                is_worm = false;
            p_type = 0.5;
        }
    }
    else if (head_idx!=-1){ // only worm head present, can insert antiworm only
        if (n==0){
            insertZero_anti_attempts += 1;
            return; // cannot insert proposed antiworm, no particles present
        }
        else{
            is_worm = false;
            p_type = 1.0;
        }
    }
    else{ // only tail present, can insert worm only
        is_worm = true;
        p_type = 1.0;
    }
    
    // Add to worm/antiworm insertion attempt counters
    if (is_worm){insertZero_worm_attempts += 1;}
    else {insertZero_anti_attempts += 1;}
    
    // Randomly choose where to insert worm end on the flat interval
    tau_new = tau_flat*rng.rand();
    
    // Determine the no. of particles after each worm end
    if (is_worm){
        n_tail = n + 1;
        n_head = n;
    }
    else{
        n_tail = n;
        n_head = n - 1;
    }
    
    // Calculate the diagonal energy difference dV = \epsilon_w - \epsilon
    dV = (U/2.0)*(n_tail*(n_tail-1)-n_head*(n_head-1)) - mu*(n_tail-n_head);
    
    // deleteZero (reverse update) might've had to choose either head or tail
    if (head_idx==-1 and tail_idx==-1) // only one end after insertZero
        p_wormend = 1.0;
    else{                              // two worm ends after insertZero
        if (is_worm){
            if (paths[paths[tail_idx].prev].tau != 0)
                p_wormend = 1.0; //cannot choose tail.was not coming from tau=0.
            else
                p_wormend = 0.5; // deleteZero could choose either head or tail
        }
        else{ // if insert anti (i.e, a tail) the end present was a head
            if (paths[paths[head_idx].prev].tau != 0)
                p_wormend = 1.0;
            else
                p_wormend = 0.5;
        }
    }
    
    // Determine the length of the path to be modified
    l_path = tau_new;
    
    // Determine the total particle change based on worm type
    if (is_worm){
        dN = +1.0 * l_path/beta;
    }
    else{
        dN = -1.0 * l_path/beta;
    }
    
    // Canonical simulations: Restrict updates to interval N:(N-1,N+1)
    if (canonical)
        if ((N_tracker+dN) < (N-1) || (N_tracker+dN) > (N+1)){return;}
    

    // Build the weight ratio W'/W
     C = 1.0;
    if (is_worm){
//        C = sqrt(N_b+1)/sqrt(n+1);
        W = eta * sqrt(n_tail) * C * exp(-dV*tau_new);
    }
    else{
//        C = sqrt(n)/sqrt(N_b);
        W = eta * sqrt(n_tail) * C * exp(dV*tau_new);
    }
    
    // Build the Metropolis Ratio (R)
    p_dz = 0.5;
    p_iz = 0.5;
    R = W * (p_dz/p_iz) * M * p_wormend * tau_flat / p_type;

    // Metropolis sampling
    if (rng.rand() < R){ // Accept
        
        // Activate the first available kink
        if (is_worm){
            paths[num_kinks] = Kink(tau_new,n_head,src,src,src,next,
                                    dest_replica,dest_replica);
            
            // Save head index
            head_idx = num_kinks;
            
            // Update the number of particles in the initial kink at site i
            paths[i].n = n_tail;
            
            // Add to Acceptance counter
            insertZero_worm_accepts += 1;
            
            // Worm inserted, add one to tau=0 particle tracker
            N_zero += 1;
        }
        else{ // antiworm
            paths[num_kinks] = Kink(tau_new,n_tail,src,src,src,next,
                                    dest_replica,dest_replica);
            
            // Save head index
            tail_idx = num_kinks;
            
            // Update number of particles in initial kink of insertion site
            paths[i].n = n_head;
            
            // Add to Acceptance counter
            insertZero_anti_accepts += 1;
            
            // Antiworm inserted, subtract one to tau=0 particle tracker
            N_zero -= 1;
        }
        
        // "Connect" next of lower bound kink to the new worm end
        paths[src].next = num_kinks;
        
        // "Connect" prev of next kink to the new worm end
        if (next!=-1){paths[next].prev = num_kinks;}
        
        // Update trackers for: no. of active kinks, total particles
        num_kinks += 1;
        N_tracker += dN;
        
        // If new worm end is last kink on site, update last_kinks vector
        if (next==-1){
            if (is_worm){last_kinks[src]=head_idx;}
            else {last_kinks[src]=tail_idx;}
        }

        return;
    }
    else // Reject
        return;
}

/*--------------------------------------------------------------------*/

void insertZero_2(vector<Kink> &paths, int &num_kinks, int &head_idx,
                int &tail_idx, int M, int N, double U, double mu, double t,
                double beta, double eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                unsigned long long int &insertZero_worm_attempts,
                unsigned long long int &insertZero_worm_accepts,
                unsigned long long int &insertZero_anti_attempts,
                unsigned long long int &insertZero_anti_accepts,
                RNG &rng){
    
    // Variable declarations
    int n,src,next,n_head,n_tail,i,dest_replica;
    double tau_flat,l_path,dN,dV,R,p_type,tau_new,p_wormend,C,
    p_dz,p_iz;
    bool is_worm;

    // Cannot insert if there's two worm ends present
    if (head_idx != -1 and tail_idx != -1){return;}
        
    // Randomly select site on which to insert worm/antiworm from tau=0
    //boost::random::uniform_int_distribution<> sites(0, M-1);
    i = rng.randInt(M-1);
    
    // Extract attributes of insertion flat
    n = paths[i].n;
    src = paths[i].src;
    next = paths[i].next;
    dest_replica = paths[i].dest_replica;
    
    // Determine the length of insertion flat interval
    if (next != -1)
        tau_flat = paths[next].tau;
    else
        tau_flat = beta;
    
    // Choose worm/antiworm insertion based on worm ends present
    //boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    if (head_idx==-1 and tail_idx==-1){ // no worm ends present
        if (n==0){ // can only insert worm, not antiworm
            is_worm = true;
            p_type = 1.0;
        }
        else{ // choose worm or antiworm insertion with equal probability
            if (rng.rand() < 0.5)
                is_worm = true;
            else
                is_worm = false;
            p_type = 0.5;
        }
    }
    else if (head_idx!=-1){ // only worm head present, can insert antiworm only
        if (n==0){
            insertZero_anti_attempts += 1;
            return; // cannot insert proposed antiworm, no particles present
        }
        else{
            is_worm = false;
            p_type = 1.0;
        }
    }
    else{ // only tail present, can insert worm only
        is_worm = true;
        p_type = 1.0;
    }
    
    // Add to worm/antiworm insertion attempt counters
    if (is_worm){insertZero_worm_attempts += 1;}
    else {insertZero_anti_attempts += 1;}
    
    // Randomly choose where to insert worm end on the flat interval
    // tau_new = tau_flat*rng.rand();

    // Determine the no. of particles after each worm end
    if (is_worm){
        n_tail = n + 1;
        n_head = n;
    }
    else{
        n_tail = n;
        n_head = n - 1;
    }
    
    // Calculate the diagonal energy difference dV = \epsilon_w - \epsilon
    dV = (U/2.0)*(n_tail*(n_tail-1)-n_head*(n_head-1)) - mu*(n_tail-n_head);
    if (!is_worm){dV *= -1;}

    /* :::::::::::::::::::::::::: Truncated Sampling :::::::::::::::::::::::: */
    // Sample time on flat interval from truncated exponential for insertion
    double x,Z,a,b,c;

    a = 0.0;
    b = tau_flat;

    c = dV; 

    x = rng.rand();
    Z = 1.0 - exp(-c*(b-a)); //
    tau_new = a - log(1.0-Z*x)  / c;
    // if (!is_worm){cout << tau_new << endl;}
    /* :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */
    
    // deleteZero (reverse update) might've had to choose either head or tail
    if (head_idx==-1 and tail_idx==-1) // only one end after insertZero
        p_wormend = 1.0;
    else{                              // two worm ends after insertZero
        if (is_worm){
            if (paths[paths[tail_idx].prev].tau != 0)
                p_wormend = 1.0; //cannot choose tail.was not coming from tau=0.
            else
                p_wormend = 0.5; // deleteZero could choose either head or tail
        }
        else{ // if insert anti (i.e, a tail) the end present was a head
            if (paths[paths[head_idx].prev].tau != 0)
                p_wormend = 1.0;
            else
                p_wormend = 0.5;
        }
    }
    
    // Determine the length of the path to be modified
    l_path = tau_new;
    
    // Determine the total particle change based on worm type
    if (is_worm){
        dN = +1.0 * l_path/beta;
    }
    else{
        dN = -1.0 * l_path/beta;
    }
    
    // Canonical simulations: Restrict updates to interval N:(N-1,N+1)
    if (canonical)
        if ((N_tracker+dN) < (N-1) || (N_tracker+dN) > (N+1)){return;}
   

    // Build the weight ratio W'/W
    C = 1.0;
//     if (is_worm){
// //        C = sqrt(N_b+1)/sqrt(n+1);
//         W = eta * sqrt(n_tail) * C * exp(-dV*tau_new);
//     }
//     else{
// //        C = sqrt(n)/sqrt(N_b);
//         W = eta * sqrt(n_tail) * C * exp(dV*tau_new);
//     }
    
    // Build the Metropolis Ratio (R)
    p_dz = 0.5;
    p_iz = 0.5;
    R = eta * sqrt(n_tail) * C * (p_dz/p_iz) * M * p_wormend * (Z/dV) / p_type;
    // if (!is_worm){cout << R << " " << Z << " " << dV << endl;}

    // Metropolis sampling
    if (rng.rand() < R){ // Accept
        
        // Activate the first available kink
        if (is_worm){
            paths[num_kinks] = Kink(tau_new,n_head,src,src,src,next,
                                    dest_replica,dest_replica);
            
            // Save head index
            head_idx = num_kinks;
            
            // Update the number of particles in the initial kink at site i
            paths[i].n = n_tail;
            
            // Add to Acceptance counter
            insertZero_worm_accepts += 1;
            
            // Worm inserted, add one to tau=0 particle tracker
            N_zero += 1;
        }
        else{ // antiworm
            paths[num_kinks] = Kink(tau_new,n_tail,src,src,src,next,
                                    dest_replica,dest_replica);
            
            // Save head index
            tail_idx = num_kinks;
            
            // Update number of particles in initial kink of insertion site
            paths[i].n = n_head;
            
            // Add to Acceptance counter
            insertZero_anti_accepts += 1;
            
            // Antiworm inserted, subtract one to tau=0 particle tracker
            N_zero -= 1;
        }
        
        // "Connect" next of lower bound kink to the new worm end
        paths[src].next = num_kinks;
        
        // "Connect" prev of next kink to the new worm end
        if (next!=-1){paths[next].prev = num_kinks;}
        
        // Update trackers for: no. of active kinks, total particles
        num_kinks += 1;
        N_tracker += dN;
        
        // If new worm end is last kink on site, update last_kinks vector
        if (next==-1){
            if (is_worm){last_kinks[src]=head_idx;}
            else {last_kinks[src]=tail_idx;}
        }

        return;
    }
    else // Reject
        return;
}

/*--------------------------------------------------------------------*/

void deleteZero(vector<Kink> &paths, int &num_kinks, int &head_idx,
                int &tail_idx, int M, int N, double U, double mu, double t,
                double beta, double eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                unsigned long long int &deleteZero_worm_attempts,
                unsigned long long int &deleteZero_worm_accepts,
                unsigned long long int &deleteZero_anti_attempts,
                unsigned long long int &deleteZero_anti_accepts,
                RNG &rng){
    
    // Variable declarations
    int n,src,prev,next,n_head,n_tail,worm_end_idx;
    double tau,tau_next,tau_flat,l_path,dN,dV,R,p_type,p_wormend,C,W,p_dz,p_iz;
    bool delete_head;

    // Cannot delete if there are no worm ends present
    if (head_idx==-1 && tail_idx==-1){return;}

    // Cannot delete if there are no worm ends coming from tau=0
    if (head_idx!=-1 && tail_idx!=-1){
        if (paths[paths[head_idx].prev].tau != 0 and
            paths[paths[tail_idx].prev].tau != 0){return;}
    }
    else if (head_idx!=-1){ // only head present
        if (paths[paths[head_idx].prev].tau != 0){return;}
    }
    else{ // only tail present
        if (paths[paths[tail_idx].prev].tau != 0){return;}
    }

    // Decide which worm end to delete
    if (head_idx!=-1 && tail_idx!=-1){ // both wormends present
        if (paths[paths[head_idx].prev].tau == 0 &&
            paths[paths[tail_idx].prev].tau == 0){ //both near 0
            if (rng.rand() < 0.5)
                delete_head = true;
            else
                delete_head = false;
            p_wormend = 0.5;
        }
        else if (paths[paths[head_idx].prev].tau == 0){
            delete_head = true;
            p_wormend = 1.0;
        }
        else{ // only the tail is near zero
            delete_head = false;
            p_wormend = 1.0;
        }
    }
    else if (head_idx!=-1){ // only head present
        delete_head = true;
        p_wormend = 1.0;
    }
    else{ // only tail present
        delete_head = false;
        p_wormend = 1.0;
    }

    // Get index of worm end to be deleted
    if (delete_head)
        worm_end_idx = head_idx;
    else
        worm_end_idx = tail_idx;

    // Extract worm end attributes
    tau = paths[worm_end_idx].tau;
    n = paths[worm_end_idx].n;
    src = paths[worm_end_idx].src;
    prev = paths[worm_end_idx].prev;
    next = paths[worm_end_idx].next;

    // Calculate the length of the flat interval (excluding the wormend)
    if (next == -1)
        tau_next = beta;
    else
        tau_next = paths[next].tau;
    tau_flat = tau_next;

    // No. of particles before,after the worm end to be deleted
    if (delete_head) // delete worm
        n_tail = n+1;
    else //  delete antiworm
        n_tail = n;
    n_head = n_tail-1;
    
    // Worm insert (reverse update) probability of choosing worm or antiworm
    if (head_idx!=-1 and tail_idx!=-1) // worm end present before insertion
        p_type = 1.0;
    else{ // no worm ends present before insertion
        if (n==0)
            p_type = 1.0; // only worm can be inserted if no particles on flat
        else
            p_type = 0.5;
    }

    // Add to deleteZero PROPOSAL counters
    if (delete_head) // delete head (delete worm)
        deleteZero_worm_attempts += 1;
    else                   // delete tail (delete antiworm)
        deleteZero_anti_attempts += 1;
    
    // Determine the length of path to be modified
    l_path = tau;
    
    // Determine the total particle change based on worm type
    if (delete_head) // delete worm
        dN = -1.0 * l_path / beta;
    else // delete anti
        dN = +1.0 * l_path / beta;
    
    // Canonical simulations: Restrict updates to interval N:(N-1,N+1)
    if (canonical)
        if ((N_tracker+dN) < (N-1) || (N_tracker+dN) > (N+1)){return;}
        
    // Calculate diagonal energy difference
    dV = (U/2.0)*(n_tail*(n_tail-1)-n_head*(n_head-1)) - mu*(n_tail-n_head);
    

    // Build the weigh ratio W'/W
     C = 1.0;
    if (delete_head){ // delete worm
//        C = sqrt(N_b+1)/sqrt(n+1);
        W = eta * sqrt(n_tail) * C * exp(-dV*tau);
    }
    else{ // delete antiworm
//        C = sqrt(n)/sqrt(N_b);
        W = eta * sqrt(n_tail) * C * exp(dV*tau);
    }
    
    // Build the Metropolis Ratio  (R)
    p_dz = 0.5;
    p_iz = 0.5;
    R = W * (p_dz/p_iz) * M * p_wormend * tau_flat / p_type;
    R = 1.0/R;
    
    // Metropolis sampling
    if (rng.rand() < R){ // accept
        
        // Update the number of particles in the initial kink of worm end site
        paths[paths[worm_end_idx].prev].n = n;
        
        // num_kinks-1 (last available kink) will be swapped. Modify links to it
        if (paths[num_kinks-1].next!=-1) // avoids access with -1 index
            paths[paths[num_kinks-1].next].prev = worm_end_idx;
        paths[paths[num_kinks-1].prev].next = worm_end_idx;
        
        swap(paths[worm_end_idx],paths[num_kinks-1]);
        
        // Upper bound of flat could've been swapped. Correct if so.
        if (next==num_kinks-1){next=worm_end_idx;}
                
        // The other worm end could've been swapped. Reindex it if so.
        if (delete_head){
            if (tail_idx==num_kinks-1){tail_idx=worm_end_idx;}
        }
        else{
            if (head_idx==num_kinks-1){head_idx=worm_end_idx;}
        }
        
        // Whatever kink was swapped could've been the last on its site.
        if (paths[worm_end_idx].next==-1){
            last_kinks[paths[worm_end_idx].src]=worm_end_idx;
        }
        
        // If deleted worm end was last kink on site, update last kinks tracker
        if (next==-1){last_kinks[src]=prev;}
        
        // Reconnect the lower,upper bounds of the flat interval.
        paths[prev].next = next;
        if (next!=-1)
            paths[next].prev = prev;

        // Deactivate the worm end
        if (delete_head){
            head_idx = -1;
            
            // Add to Acceptance counter
            deleteZero_worm_accepts += 1;
            
            // Worm deleted, subtract one to tau=0 particle tracker
            N_zero -= 1;
        }
        else{
            tail_idx = -1;
            
            // Add to Acceptance counter
            deleteZero_anti_accepts += 1;
            
            // Antiworm deleted, add one to tau=0 particle tracker
            N_zero += 1;
        }
        
        // Update trackers for: num of active kinks, total particles
        num_kinks -= 1;
        N_tracker += dN;
        
        return;
    }
    else // reject
        return;
}
/*--------------------------------------------------------------------*/

void deleteZero_2(vector<Kink> &paths, int &num_kinks, int &head_idx,
                int &tail_idx, int M, int N, double U, double mu, double t,
                double beta, double eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                unsigned long long int &deleteZero_worm_attempts,
                unsigned long long int &deleteZero_worm_accepts,
                unsigned long long int &deleteZero_anti_attempts,
                unsigned long long int &deleteZero_anti_accepts,
                RNG &rng){
    
    // Variable declarations
    int n,src,prev,next,n_head,n_tail,worm_end_idx;
    double tau,tau_next,l_path,dN,dV,R,p_type,p_wormend,
    C,p_dz,p_iz,Z;
    bool delete_head;

    // Cannot delete if there are no worm ends present
    if (head_idx==-1 && tail_idx==-1){return;}

    // Cannot delete if there are no worm ends coming from tau=0
    if (head_idx!=-1 && tail_idx!=-1){
        if (paths[paths[head_idx].prev].tau != 0 and
            paths[paths[tail_idx].prev].tau != 0){return;}
    }
    else if (head_idx!=-1){ // only head present
        if (paths[paths[head_idx].prev].tau != 0){return;}
    }
    else{ // only tail present
        if (paths[paths[tail_idx].prev].tau != 0){return;}
    }

    // Decide which worm end to delete
    if (head_idx!=-1 && tail_idx!=-1){ // both wormends present
        if (paths[paths[head_idx].prev].tau == 0 &&
            paths[paths[tail_idx].prev].tau == 0){ //both near 0
            if (rng.rand() < 0.5)
                delete_head = true;
            else
                delete_head = false;
            p_wormend = 0.5;
        }
        else if (paths[paths[head_idx].prev].tau == 0){
            delete_head = true;
            p_wormend = 1.0;
        }
        else{ // only the tail is near zero
            delete_head = false;
            p_wormend = 1.0;
        }
    }
    else if (head_idx!=-1){ // only head present
        delete_head = true;
        p_wormend = 1.0;
    }
    else{ // only tail present
        delete_head = false;
        p_wormend = 1.0;
    }

    // Get index of worm end to be deleted
    if (delete_head)
        worm_end_idx = head_idx;
    else
        worm_end_idx = tail_idx;

    // Extract worm end attributes
    tau = paths[worm_end_idx].tau;
    n = paths[worm_end_idx].n;
    src = paths[worm_end_idx].src;
    prev = paths[worm_end_idx].prev;
    next = paths[worm_end_idx].next;

    // Calculate the length of the flat interval (excluding the wormend)
    if (next == -1)
        tau_next = beta;
    else
        tau_next = paths[next].tau;

    // No. of particles before,after the worm end to be deleted
    if (delete_head) // delete worm
        n_tail = n+1;
    else //  delete antiworm
        n_tail = n;
    n_head = n_tail-1;
    
    // Worm insert (reverse update) probability of choosing worm or antiworm
    if (head_idx!=-1 and tail_idx!=-1) // worm end present before insertion
        p_type = 1.0;
    else{ // no worm ends present before insertion
        if (n==0)
            p_type = 1.0; // only worm can be inserted if no particles on flat
        else
            p_type = 0.5;
    }

    // Add to deleteZero PROPOSAL counters
    if (delete_head) // delete head (delete worm)
        deleteZero_worm_attempts += 1;
    else                   // delete tail (delete antiworm)
        deleteZero_anti_attempts += 1;
    
    // Determine the length of path to be modified
    l_path = tau;
    
    // Determine the total particle change based on worm type
    if (delete_head) // delete worm
        dN = -1.0 * l_path / beta;
    else // delete anti
        dN = +1.0 * l_path / beta;
    
    // Canonical simulations: Restrict updates to interval N:(N-1,N+1)
    if (canonical)
        if ((N_tracker+dN) < (N-1) || (N_tracker+dN) > (N+1)){return;}
        
    // Calculate diagonal energy difference
    dV = (U/2.0)*(n_tail*(n_tail-1)-n_head*(n_head-1)) - mu*(n_tail-n_head);
    
  
//     // Build the weigh ratio W'/W
     C = 1.0;
//     if (delete_head){ // delete worm
// //        C = sqrt(N_b+1)/sqrt(n+1);
//         W = eta * sqrt(n_tail) * C * exp(-dV*tau);
//     }
//     else{ // delete antiworm
// //        C = sqrt(n)/sqrt(N_b);
//         W = eta * sqrt(n_tail) * C * exp(dV*tau);
//     }

    // Inverse move (insertZero) truncated sampling
    if (!delete_head){dV *= -1;}
    Z = 1.0 - exp(-dV*(tau_next));
    
    // Build the Metropolis Ratio  (R)
    p_dz = 0.5;
    p_iz = 0.5;
    // R = W * (p_dz/p_iz) * M * p_wormend * (Z/dV) / p_type;
    R = eta * sqrt(n_tail) * C * (p_dz/p_iz) * M * p_wormend * (Z/dV) / p_type;
    R = 1.0/R;
    
    // Metropolis sampling
    if (rng.rand() < R){ // accept
        
        // Update the number of particles in the initial kink of worm end site
        paths[paths[worm_end_idx].prev].n = n;
        
        // num_kinks-1 (last available kink) will be swapped. Modify links to it
        if (paths[num_kinks-1].next!=-1) // avoids access with -1 index
            paths[paths[num_kinks-1].next].prev = worm_end_idx;
        paths[paths[num_kinks-1].prev].next = worm_end_idx;
        
        swap(paths[worm_end_idx],paths[num_kinks-1]);
        
        // Upper bound of flat could've been swapped. Correct if so.
        if (next==num_kinks-1){next=worm_end_idx;}
                
        // The other worm end could've been swapped. Reindex it if so.
        if (delete_head){
            if (tail_idx==num_kinks-1){tail_idx=worm_end_idx;}
        }
        else{
            if (head_idx==num_kinks-1){head_idx=worm_end_idx;}
        }
        
        // Whatever kink was swapped could've been the last on its site.
        if (paths[worm_end_idx].next==-1){
            last_kinks[paths[worm_end_idx].src]=worm_end_idx;
        }
        
        // If deleted worm end was last kink on site, update last kinks tracker
        if (next==-1){last_kinks[src]=prev;}
        
        // Reconnect the lower,upper bounds of the flat interval.
        paths[prev].next = next;
        if (next!=-1)
            paths[next].prev = prev;

        // Deactivate the worm end
        if (delete_head){
            head_idx = -1;
            
            // Add to Acceptance counter
            deleteZero_worm_accepts += 1;
            
            // Worm deleted, subtract one to tau=0 particle tracker
            N_zero -= 1;
        }
        else{
            tail_idx = -1;
            
            // Add to Acceptance counter
            deleteZero_anti_accepts += 1;
            
            // Antiworm deleted, add one to tau=0 particle tracker
            N_zero += 1;
        }
        
        // Update trackers for: num of active kinks, total particles
        num_kinks -= 1;
        N_tracker += dN;
        
        return;
    }
    else // reject
        return;
}
/*--------------------------------------------------------------------*/

void insertBeta(vector<Kink> &paths, int &num_kinks, int &head_idx,
                int &tail_idx, int M, int N, double U, double mu, double t,
                double beta, double eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                unsigned long long int &insertBeta_worm_attempts,
                unsigned long long int &insertBeta_worm_accepts,
                unsigned long long int &insertBeta_anti_attempts,
                unsigned long long int &insertBeta_anti_accepts,
                RNG &rng){
    
    // Variable declarations
    int n,src,next,n_head,n_tail,i,src_replica;
    double tau_prev,tau_flat,
    l_path,dN,dV,R,p_type,tau_new,p_wormend,C,W,p_db,p_ib;
    bool is_worm;

    // Cannot insert if there's two worm ends present
    if (head_idx != -1 and tail_idx != -1){return;}
        
    // Randomly select site on which to insert worm/antiworm from tau=0
    //boost::random::uniform_int_distribution<> sites(0, M-1);
    i = rng.randInt(M-1);
    
    // Extract the flat interval where insertion is proposed & its attributes
    tau_prev = paths[last_kinks[i]].tau;
    n = paths[last_kinks[i]].n;
    src = paths[last_kinks[i]].src;
    next = paths[last_kinks[i]].next;
    src_replica = paths[last_kinks[i]].src_replica;
    
    // Determine the length of insertion flat interval
    tau_flat = beta - tau_prev;
    
    // Choose worm/antiworm insertion based on worm ends present
    if (head_idx==-1 and tail_idx==-1){ // no worm ends present
        if (n==0){ // can only insert worm, not antiworm
            is_worm = true;
            p_type = 1.0;
        }
        else{ // choose worm or antiworm insertion with equal probability
            if (rng.rand() < 0.5)
                is_worm = true;
            else
                is_worm = false;
            p_type = 0.5;
        }
    }
    else if (tail_idx!=-1){ // only worm tail present, can insert antiworm only
        if (n==0){
            insertBeta_anti_attempts += 1.0;
            return; // cannot insert proposed antiworm, no particles present
        }
        else{
            is_worm = false;
            p_type = 1.0;
        }
    }
    else{ // only head present, can insert worm only
        is_worm = true;
        p_type = 1.0;
    }
    
    // Add to worm/antiworm insertion attempt counters
    if (is_worm){insertBeta_worm_attempts += 1;}
    else {insertBeta_anti_attempts += 1;}
    
    // Randomly choose where to insert worm end on the flat interval
    tau_new = tau_prev + tau_flat*rng.rand();
    
    // Determine the no. of particles after each worm end
    if (is_worm){
        n_tail = n + 1;
        n_head = n;
    }
    else{
        n_tail = n;
        n_head = n - 1;
    }
    
    // Calculate the diagonal energy difference dV = \epsilon_w - \epsilon
    dV = (U/2.0)*(n_tail*(n_tail-1)-n_head*(n_head-1)) - mu*(n_tail-n_head);
    
    // deleteBeta (reverse update) might've had to choose either head or tail
    if (head_idx==-1 and tail_idx==-1) // only one end after insertZero
        p_wormend = 1.0;
    else{                              // two worm ends after insertZero
        if (is_worm){
            if (paths[head_idx].next != -1)
                p_wormend = 1.0; // cannot choose head.was not coming from beta.
            else
                p_wormend = 0.5; // deleteBeta could choose either head or tail
        }
        else{ // if insert anti (i.e, a head) the end present was a tail
            if (paths[tail_idx].next != -1)
                p_wormend = 1.0;
            else
                p_wormend = 0.5;
        }
    }
    
    // Determine the length of the path to be modified
    l_path = beta - tau_new;
    
    // Determine the total particle change based on worm type
    if (is_worm){
        dN = +1.0 * l_path/beta;
    }
    else{
        dN = -1.0 * l_path/beta;
    }
    
    // Canonical simulations: Restrict updates to interval N:(N-1,N+1)
    if (canonical)
        if ((N_tracker+dN) < (N-1) || (N_tracker+dN) > (N+1)){return;}
    
    
    // Build the weight ratio W'/W
     C = 1.0; // C_pre/C_post
    if (is_worm){
//        C = sqrt(N_b+1)/sqrt(n+1);
        W = eta * sqrt(n_tail) * C * exp(-dV*(beta-tau_new));
    }
    else{
//        C = sqrt(n)/sqrt(N_b);
        W = eta * sqrt(n_tail) * C * exp(-dV*(tau_new-beta));
    }
    
    // Build the Metropolis Ratio (R)
    p_db = 0.5;
    p_ib = 0.5;
    R = W * (p_db/p_ib) * M * p_wormend * tau_flat / p_type;

    // Metropolis sampling
    if (rng.rand() < R){ // Accept
        
        // Activate the first available kink
        if (is_worm){
            paths[num_kinks] = Kink (tau_new,n_tail,src,src,
                                            last_kinks[src],next,
                                            src_replica,src_replica);
            
            // Save tail index
            tail_idx = num_kinks;
            
            // Add to Acceptance counter
            insertBeta_worm_accepts += 1;
            
            // Worm inserted, add one to tau=beta particle tracker
            N_beta += 1;
        }
        else{ // antiworm
            paths[num_kinks] = Kink (tau_new,n_head,src,src,
                                            last_kinks[src],next,
                                            src_replica,src_replica);
            
            // Save head index
            head_idx = num_kinks;
            
            // Add to Acceptance counter
            insertBeta_anti_accepts += 1;
            
            // Antiworm inserted, subtract one to tau=beta particle tracker
            N_beta -= 1;
        }
        
        // "Connect" next of lower bound kink to the new worm end
        paths[last_kinks[src]].next = num_kinks;
        
        // Update trackers for: no. of active kinks, total particles
        num_kinks += 1;
        N_tracker += dN;
        
        // If new worm end is last kink on site, update last_kinks vector
        if (next==-1){
            if (is_worm){last_kinks[src]=tail_idx;}
            else {last_kinks[src]=head_idx;}
        }
        return;
    }
    else // Reject
        return;
}

/*--------------------------------------------------------------------*/

void insertBeta_2(vector<Kink> &paths, int &num_kinks, int &head_idx,
                int &tail_idx, int M, int N, double U, double mu, double t,
                double beta, double eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                unsigned long long int &insertBeta_worm_attempts,
                unsigned long long int &insertBeta_worm_accepts,
                unsigned long long int &insertBeta_anti_attempts,
                unsigned long long int &insertBeta_anti_accepts,
                RNG &rng){
    
    // Variable declarations
    int n,src,next,n_head,n_tail,i,src_replica;
    double tau_prev,
    l_path,dN,dV,R,p_type,tau_new,p_wormend,C,p_db,p_ib;
    bool is_worm; 

    // Cannot insert if there's two worm ends present
    if (head_idx != -1 and tail_idx != -1){return;}
        
    // Randomly select site on which to insert worm/antiworm from tau=0
    //boost::random::uniform_int_distribution<> sites(0, M-1);
    i = rng.randInt(M-1);
    
    // Extract the flat interval where insertion is proposed & its attributes
    tau_prev = paths[last_kinks[i]].tau;
    n = paths[last_kinks[i]].n;
    src = paths[last_kinks[i]].src;
    next = paths[last_kinks[i]].next;
    src_replica = paths[last_kinks[i]].src_replica;
     
    // Choose worm/antiworm insertion based on worm ends present
    if (head_idx==-1 and tail_idx==-1){ // no worm ends present
        if (n==0){ // can only insert worm, not antiworm
            is_worm = true;
            p_type = 1.0;
        }
        else{ // choose worm or antiworm insertion with equal probability
            if (rng.rand() < 0.5)
                is_worm = true;
            else
                is_worm = false;
            p_type = 0.5;
        }
    }
    else if (tail_idx!=-1){ // only worm tail present, can insert antiworm only
        if (n==0){
            insertBeta_anti_attempts += 1.0;
            return; // cannot insert proposed antiworm, no particles present
        }
        else{
            is_worm = false;
            p_type = 1.0;
        }
    }
    else{ // only head present, can insert worm only
        is_worm = true;
        p_type = 1.0;
    }
    
    // Add to worm/antiworm insertion attempt counters
    if (is_worm){insertBeta_worm_attempts += 1;}
    else {insertBeta_anti_attempts += 1;}
    
    // Randomly choose where to insert worm end on the flat interval
    // tau_new = tau_prev + tau_flat*rng.rand();
    
    // Determine the no. of particles after each worm end
    if (is_worm){
        n_tail = n + 1;
        n_head = n;
    }
    else{
        n_tail = n;
        n_head = n - 1;
    }
    
    // Calculate the diagonal energy difference dV = \epsilon_w - \epsilon
    dV = (U/2.0)*(n_tail*(n_tail-1)-n_head*(n_head-1)) - mu*(n_tail-n_head);
    if (!is_worm){dV *= -1;}

    /* :::::::::::::::::::::::::: Truncated Sampling :::::::::::::::::::::::: */
    // Sample time on flat interval from truncated exponential for insertion
    double x,Z,a,b,c;

    a = tau_prev;
    b = beta;

    c = dV; 

    x = rng.rand();
    Z = 1.0 - exp(-c*(b-a)); //
    tau_new = b + log(1.0-Z*x)  / c;
    // if (!is_worm){cout << tau_new << endl;}
    /* :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */
    
    // deleteBeta (reverse update) might've had to choose either head or tail
    if (head_idx==-1 and tail_idx==-1) // only one end after insertZero
        p_wormend = 1.0;
    else{                              // two worm ends after insertZero
        if (is_worm){
            if (paths[head_idx].next != -1)
                p_wormend = 1.0; // cannot choose head.was not coming from beta.
            else
                p_wormend = 0.5; // deleteBeta could choose either head or tail
        }
        else{ // if insert anti (i.e, a head) the end present was a tail
            if (paths[tail_idx].next != -1)
                p_wormend = 1.0;
            else
                p_wormend = 0.5;
        }
    }
    
    // Determine the length of the path to be modified
    l_path = beta - tau_new;
    
    // Determine the total particle change based on worm type
    if (is_worm){
        dN = +1.0 * l_path/beta;
    }
    else{
        dN = -1.0 * l_path/beta;
    }
    
    // Canonical simulations: Restrict updates to interval N:(N-1,N+1)
    if (canonical)
        if ((N_tracker+dN) < (N-1) || (N_tracker+dN) > (N+1)){return;}
    
    
    // Build the weight ratio W'/W
     C = 1.0; // C_pre/C_post
    if (is_worm){
//        C = sqrt(N_b+1)/sqrt(n+1);
        // W = eta * sqrt(n_tail) * C * exp(-dV*(beta-tau_new));
    }
    else{
//        C = sqrt(n)/sqrt(N_b);
        // W = eta * sqrt(n_tail) * C * exp(dV*(beta-tau_new));
    }
    
    // Build the Metropolis Ratio (R)
    p_db = 0.5;
    p_ib = 0.5;
    R = eta * sqrt(n_tail) * C * (p_db/p_ib) * M * p_wormend * (Z/dV) / p_type;

    // Metropolis sampling
    if (rng.rand() < R){ // Accept
        
        // Activate the first available kink
        if (is_worm){
            paths[num_kinks] = Kink (tau_new,n_tail,src,src,
                                            last_kinks[src],next,
                                            src_replica,src_replica);
            
            // Save tail index
            tail_idx = num_kinks;
            
            // Add to Acceptance counter
            insertBeta_worm_accepts += 1;
            
            // Worm inserted, add one to tau=beta particle tracker
            N_beta += 1;
        }
        else{ // antiworm
            paths[num_kinks] = Kink (tau_new,n_head,src,src,
                                            last_kinks[src],next,
                                            src_replica,src_replica);
            
            // Save head index
            head_idx = num_kinks;
            
            // Add to Acceptance counter
            insertBeta_anti_accepts += 1;
            
            // Antiworm inserted, subtract one to tau=beta particle tracker
            N_beta -= 1;
        }
        
        // "Connect" next of lower bound kink to the new worm end
        paths[last_kinks[src]].next = num_kinks;
        
        // Update trackers for: no. of active kinks, total particles
        num_kinks += 1;
        N_tracker += dN;
        
        // If new worm end is last kink on site, update last_kinks vector
        if (next==-1){
            if (is_worm){last_kinks[src]=tail_idx;}
            else {last_kinks[src]=head_idx;}
        }
        return;
    }
    else // Reject
        return;
}

/*--------------------------------------------------------------------*/

void deleteBeta(vector<Kink> &paths, int &num_kinks, int &head_idx,
                int &tail_idx, int M, int N, double U, double mu, double t,
                double beta, double eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                unsigned long long int &deleteBeta_worm_attempts,
                unsigned long long int &deleteBeta_worm_accepts,
                unsigned long long int &deleteBeta_anti_attempts,
                unsigned long long int &deleteBeta_anti_accepts,
                RNG &rng){
    
    // Variable declarations
    int n,src,prev,next,n_head,n_tail,worm_end_idx;
    double tau,tau_prev,tau_flat,l_path,dN,dV,R,p_type,p_wormend,C,W,p_db,p_ib;
    bool delete_head;

    // Cannot delete if there are no worm ends present
    if (head_idx==-1 && tail_idx==-1){return;}

    // Cannot delete if there are no worm ends coming from tau=beta
    if (head_idx!=-1 && tail_idx!=-1){
        if (paths[head_idx].next != -1 and
            paths[tail_idx].next != -1){return;}
    }
    else if (head_idx!=-1){ // only head present
        if (paths[head_idx].next != -1){return;}
    }
    else{ // only tail present
        if (paths[tail_idx].next != -1){return;}
    }
    
    // Decide which worm end to delete
    //boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    if (head_idx!=-1 && tail_idx!=-1){ // both wormends present
        if (paths[head_idx].next == -1 &&
            paths[tail_idx].next == -1){ // both last
            if (rng.rand() < 0.5)
                delete_head = true;  // delete antiworm
            else
                delete_head = false; // delete worm
            p_wormend = 0.5;
        }
        else if (paths[head_idx].next == -1){
            delete_head = true;
            p_wormend = 1.0;
        }
        else{ // only the tail is near zero
            delete_head = false;
            p_wormend = 1.0;
        }
    }
    else if (head_idx!=-1){ // only head present
        delete_head = true;
        p_wormend = 1.0;
    }
    else{ // only tail present
        delete_head = false;
        p_wormend = 1.0;
    }

    // Get index of worm end to be deleted
    if (delete_head)
        worm_end_idx = head_idx;
    else
        worm_end_idx = tail_idx;

    // Extract worm end attributes
    tau = paths[worm_end_idx].tau;
    n = paths[worm_end_idx].n;
    src = paths[worm_end_idx].src;
    prev = paths[worm_end_idx].prev;
    next = paths[worm_end_idx].next;

    // Calculate the length of the flat interval (excluding the worm end)
    tau_prev = paths[prev].tau;
    tau_flat = beta - tau_prev;

    // No. of particles before,after the worm end to be deleted
    if (delete_head) // delete antiworm
        n_tail = n+1;
    else //  delete worm
        n_tail = n;
    n_head = n_tail-1;
    
    // Worm insert (reverse update) probability of choosing worm or antiworm
    if (head_idx!=-1 and tail_idx!=-1) // worm end present before insertion
        p_type = 1.0;
    else{ // no worm ends present before insertion
        if (paths[paths[worm_end_idx].prev].n==0)
            p_type = 1.0; // only worm can be inserted if no particles on flat
        else
            p_type = 0.5;
    }

    // Add to deleteBeta PROPOSAL counters
    if (delete_head) // delete antiworm
        deleteBeta_anti_attempts += 1;
    else                   // delete worm
        deleteBeta_worm_attempts += 1;
    
    // Determine the length of path to be modified
    l_path = beta - tau;
    
    // Determine the total particle change based on worm type
    if (delete_head) // delete antiworm
        dN = +1.0 * l_path / beta;
    else // delete worm
        dN = -1.0 * l_path / beta;
    
    // Canonical simulations: Restrict updates to interval N:(N-1,N+1)
    if (canonical)
        if ((N_tracker+dN) < (N-1) || (N_tracker+dN) > (N+1)){return;}
    
    // Calculate diagonal energy difference
    dV = (U/2.0)*(n_tail*(n_tail-1)-n_head*(n_head-1)) - mu*(n_tail-n_head);
    
 
    // Build the weigh ratio W'/W
    C = 1.0;
    if (!delete_head){ // delete worm
//        C = sqrt(N_b+1)/sqrt(n);
        W = eta * sqrt(n_tail) * C * exp(-dV*(beta-tau));
    }
    else{ // delete antiworm
//        C = sqrt(n+1)/sqrt(N_b);
        W = eta * sqrt(n_tail) * C * exp(-dV*(tau-beta));
    }
    
    // Build the Metropolis Ratio  (R)
    p_db = 0.5;
    p_ib = 0.5;
    R = W * (p_db/p_ib) * M * p_wormend * tau_flat / p_type;
    R = 1.0/R;
    
    // Metropolis sampling
    if (rng.rand() < R){ // accept
                
        // num_kinks-1 (last available kink) will be swapped. Modify links to it
        if (paths[num_kinks-1].next!=-1) // avoids access with -1 index
            paths[paths[num_kinks-1].next].prev = worm_end_idx;
        paths[paths[num_kinks-1].prev].next = worm_end_idx;
        
        swap(paths[worm_end_idx],paths[num_kinks-1]);
        
        // Lower bound of flat could've been swapped. Correct if so.
        if (prev==num_kinks-1){prev=worm_end_idx;}
        
        // The other worm end could've been swapped. Reindex it if so.
        if (delete_head){
            if (tail_idx==num_kinks-1){tail_idx=worm_end_idx;}
        }
        else{
            if (head_idx==num_kinks-1){head_idx=worm_end_idx;}
        }
        
        // Whatever kink was swapped could've been the last on its site.
        if (paths[worm_end_idx].next==-1){
            last_kinks[paths[worm_end_idx].src]=worm_end_idx;
        }
        
        // Deleted worm end was last kink on site, update last kinks tracker
        last_kinks[src]=prev;
        
        // Reconnect the lower,upper bounds of the flat interval.
        paths[prev].next = next;

        // Deactivate the worm end
        if (delete_head){
            head_idx = -1;
            
            // Add to Acceptance counter
            deleteBeta_anti_accepts += 1;
            
            // Antiworm deleted, add one to tau=beta particle tracker
            N_beta += 1;
        }
        else{
            tail_idx = -1;
            
            // Add to Acceptance counter
            deleteBeta_worm_accepts += 1;
            
            // Worm deleted, subtracts one to tau=beta particle tracker
            N_beta -= 1;
        }
        
        // Update trackers for: num of active kinks, total particles
        num_kinks -= 1;
        N_tracker += dN;
        
        return;
    }
    else // reject
        return;
}
/*--------------------------------------------------------------------*/

void deleteBeta_2(vector<Kink> &paths, int &num_kinks, int &head_idx,
                int &tail_idx, int M, int N, double U, double mu, double t,
                double beta, double eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                unsigned long long int &deleteBeta_worm_attempts,
                unsigned long long int &deleteBeta_worm_accepts,
                unsigned long long int &deleteBeta_anti_attempts,
                unsigned long long int &deleteBeta_anti_accepts,
                RNG &rng){
    
    // Variable declarations
    int n,src,prev,next,n_head,n_tail,worm_end_idx;
    double tau,tau_prev,l_path,dN,dV,R,p_type,p_wormend,C,p_db,p_ib;
    bool delete_head;

    // Cannot delete if there are no worm ends present
    if (head_idx==-1 && tail_idx==-1){return;}

    // Cannot delete if there are no worm ends coming from tau=beta
    if (head_idx!=-1 && tail_idx!=-1){
        if (paths[head_idx].next != -1 and
            paths[tail_idx].next != -1){return;}
    }
    else if (head_idx!=-1){ // only head present
        if (paths[head_idx].next != -1){return;}
    }
    else{ // only tail present
        if (paths[tail_idx].next != -1){return;}
    }
    
    // Decide which worm end to delete
    //boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    if (head_idx!=-1 && tail_idx!=-1){ // both wormends present
        if (paths[head_idx].next == -1 &&
            paths[tail_idx].next == -1){ // both last
            if (rng.rand() < 0.5)
                delete_head = true;  // delete antiworm
            else
                delete_head = false; // delete worm
            p_wormend = 0.5;
        }
        else if (paths[head_idx].next == -1){
            delete_head = true;
            p_wormend = 1.0;
        }
        else{ // only the tail is near zero
            delete_head = false;
            p_wormend = 1.0;
        }
    }
    else if (head_idx!=-1){ // only head present
        delete_head = true;
        p_wormend = 1.0;
    }
    else{ // only tail present
        delete_head = false;
        p_wormend = 1.0;
    }

    // Get index of worm end to be deleted
    if (delete_head)
        worm_end_idx = head_idx;
    else
        worm_end_idx = tail_idx;

    // Extract worm end attributes
    tau = paths[worm_end_idx].tau;
    n = paths[worm_end_idx].n;
    src = paths[worm_end_idx].src;
    prev = paths[worm_end_idx].prev;
    next = paths[worm_end_idx].next;

    // Calculate the length of the flat interval (excluding the worm end)
    tau_prev = paths[prev].tau;

    // No. of particles before,after the worm end to be deleted
    if (delete_head) // delete antiworm
        n_tail = n+1;
    else //  delete worm
        n_tail = n;
    n_head = n_tail-1;
    
    // Worm insert (reverse update) probability of choosing worm or antiworm
    if (head_idx!=-1 and tail_idx!=-1) // worm end present before insertion
        p_type = 1.0;
    else{ // no worm ends present before insertion
        if (paths[paths[worm_end_idx].prev].n==0)
            p_type = 1.0; // only worm can be inserted if no particles on flat
        else
            p_type = 0.5;
    }

    // Add to deleteBeta PROPOSAL counters
    if (delete_head) // delete antiworm
        deleteBeta_anti_attempts += 1;
    else                   // delete worm
        deleteBeta_worm_attempts += 1;
    
    // Determine the length of path to be modified
    l_path = beta - tau;
    
    // Determine the total particle change based on worm type
    if (delete_head) // delete antiworm
        dN = +1.0 * l_path / beta;
    else // delete worm
        dN = -1.0 * l_path / beta;
    
    // Canonical simulations: Restrict updates to interval N:(N-1,N+1)
    if (canonical)
        if ((N_tracker+dN) < (N-1) || (N_tracker+dN) > (N+1)){return;}
    
    // Calculate diagonal energy difference
    dV = (U/2.0)*(n_tail*(n_tail-1)-n_head*(n_head-1)) - mu*(n_tail-n_head);
    if (delete_head){dV *= -1;}
     
    // Build the weigh ratio W'/W
    C = 1.0;
    if (!delete_head){ // delete worm
//        C = sqrt(N_b+1)/sqrt(n);
        // W = eta * sqrt(n_tail) * C * exp(-dV*(beta-tau));
    }
    else{ // delete antiworm
//        C = sqrt(n+1)/sqrt(N_b);
        // W = eta * sqrt(n_tail) * C * exp(-dV*(tau-beta));
    }
    
    // inverse move (insertBeta) truncated exponential sampling
    double Z,a,b,c;
    a = tau_prev;
    b = beta;
    c = dV; 
    Z = 1.0 - exp(-c*(b-a));

    // Build the Metropolis Ratio  (R)
    p_db = 0.5;
    p_ib = 0.5;
    R = eta * sqrt(n_tail) * C * (p_db/p_ib) * M * p_wormend * (Z/dV) / p_type;
    R = 1.0/R;
    
    // Metropolis sampling
    if (rng.rand() < R){ // accept
                
        // num_kinks-1 (last available kink) will be swapped. Modify links to it
        if (paths[num_kinks-1].next!=-1) // avoids access with -1 index
            paths[paths[num_kinks-1].next].prev = worm_end_idx;
        paths[paths[num_kinks-1].prev].next = worm_end_idx;
        
        swap(paths[worm_end_idx],paths[num_kinks-1]);
        
        // Lower bound of flat could've been swapped. Correct if so.
        if (prev==num_kinks-1){prev=worm_end_idx;}
        
        // The other worm end could've been swapped. Reindex it if so.
        if (delete_head){
            if (tail_idx==num_kinks-1){tail_idx=worm_end_idx;}
        }
        else{
            if (head_idx==num_kinks-1){head_idx=worm_end_idx;}
        }
        
        // Whatever kink was swapped could've been the last on its site.
        if (paths[worm_end_idx].next==-1){
            last_kinks[paths[worm_end_idx].src]=worm_end_idx;
        }
        
        // Deleted worm end was last kink on site, update last kinks tracker
        last_kinks[src]=prev;
        
        // Reconnect the lower,upper bounds of the flat interval.
        paths[prev].next = next;

        // Deactivate the worm end
        if (delete_head){
            head_idx = -1;
            
            // Add to Acceptance counter
            deleteBeta_anti_accepts += 1;
            
            // Antiworm deleted, add one to tau=beta particle tracker
            N_beta += 1;
        }
        else{
            tail_idx = -1;
            
            // Add to Acceptance counter
            deleteBeta_worm_accepts += 1;
            
            // Worm deleted, subtracts one to tau=beta particle tracker
            N_beta -= 1;
        }
        
        // Update trackers for: num of active kinks, total particles
        num_kinks -= 1;
        N_tracker += dN;
        
        return;
    }
    else // reject
        return;
}
/*--------------------------------------------------------------------*/

void timeshift_uniform(vector<Kink> &paths, int &num_kinks,int &head_idx,
                int &tail_idx, int M, int N, double U, double mu, double t,
                double beta, double eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                unsigned long long int &advance_head_attempts,
                unsigned long long int &advance_head_accepts,
                unsigned long long int &recede_head_attempts,
                unsigned long long int &recede_head_accepts,
                unsigned long long int &advance_tail_attempts,
                unsigned long long int &advance_tail_accepts,
                unsigned long long int &recede_tail_attempts,
                unsigned long long int &recede_tail_accepts,
                RNG &rng){
    
    // Variable declarations
    int n,prev,next,worm_end_idx;
    double tau,tau_prev,tau_next,tau_flat,l_path,dN,dV,R,tau_new,W;
    bool shift_head;
    
    // Reject update if there are is no worm end present
    if (head_idx==-1 && tail_idx==-1){return;}

    // Choose which worm end to move
    //boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    if (head_idx!=-1 && tail_idx!=-1){ // both worm ends present

        // Randomly choose to shift HEAD or TAIL
        if (rng.rand() < 0.5)
            shift_head = true;
        else
            shift_head = false;
        }
    else if (head_idx!=-1){ // only head present
        shift_head = true;
    }
    else{ // only tail present
        shift_head = false;
    }
    
    // Save the kink index of the end that will be shifted
    if (shift_head){worm_end_idx=head_idx;}
    else {worm_end_idx=tail_idx;}
    
    // Extract worm end attributes
    tau = paths[worm_end_idx].tau;
    n = paths[worm_end_idx].n;
    prev = paths[worm_end_idx].prev;
    next = paths[worm_end_idx].next;
    
    // Measure diagonal energy difference dV
    if (shift_head)
        dV=U*n-mu;
    else
        dV=U*(n-1)-mu;
    
    // Determine the lower and upper bounds of the worm end to be timeshifted
    if (next==-1)
        tau_next = beta;
    else
        tau_next = paths[next].tau;
    tau_prev = paths[prev].tau;
    
    // Calculate length of flat interval
    tau_flat = tau_next - tau_prev;

    // Sample the new time of the worm end from uniform distribution
    tau_new = tau_prev + tau_flat*rng.rand();
    
    // Add to PROPOSAL counter
    if (shift_head){
        if (tau_new > tau){advance_head_attempts+=1;}
        else{recede_head_attempts+=1;}
    }
    else{ // shift tail
        if (tau_new > tau){advance_tail_attempts+=1;}
        else{recede_tail_attempts+=1;}
    }
    
    // Determine the length of path to be modified
    l_path = tau_new - tau;
    
    // Determine the total particle change based on wormend to be shifted
    if (shift_head)
        dN = +1.0 * l_path/beta;
    else // shift tail
        dN = -1.0 * l_path/beta;
    
    // Canonical simulations: Restrict updates to interval N:(N-1,N+1)
    if (canonical)
        if ((N_tracker+dN) < (N-1) || (N_tracker+dN) > (N+1)){return;}
    
    // Build the Metropolis condition (R)
    if (shift_head){
        W = exp(-dV*(tau_new-tau));
    }
    else{
        W = exp(-dV*(tau-tau_new));
    }
    R = W; // recall that timeshift is it's own inverse

    // Metropolis sampling
    if (rng.rand() < R){
        
        // Add to ACCEPTANCE counter
        if (shift_head){
            if (tau_new > tau){advance_head_accepts+=1;}
            else{recede_head_accepts+=1;}
        }
        else{ // shift tail
            if (tau_new > tau){advance_tail_accepts+=1;}
            else{recede_tail_accepts+=1;}
        }
        
        // Modify the worm end time
        paths[worm_end_idx].tau = tau_new;
        
        // Modify total particle number tracker
        N_tracker += dN;
        
        return;
    }
    else // Reject
        return;
}

/*--------------------------------------------------------------------*/

void timeshift(vector<Kink> &paths, int &num_kinks, int &head_idx,
                int &tail_idx, int M, int N, double U, double mu, double t,
                double beta, double eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                unsigned long long int &advance_head_attempts,
                unsigned long long int &advance_head_accepts,
                unsigned long long int &recede_head_attempts,
                unsigned long long int &recede_head_accepts,
                unsigned long long int &advance_tail_attempts,
                unsigned long long int &advance_tail_accepts,
                unsigned long long int &recede_tail_attempts,
                unsigned long long int &recede_tail_accepts,
                RNG &rng){
    
    // Variable declarations
    int n,prev,next,worm_end_idx;
    double tau,tau_prev,tau_next,l_path,dN,dV,R,tau_new,Z;
    bool shift_head;
    
    // Reject update if there is no worm end present
    if (head_idx==-1 && tail_idx==-1){return;}

    // Choose which worm end to move
    //boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    if (head_idx!=-1 && tail_idx!=-1){ // both worm ends present

        // Randomly choose to shift HEAD or TAIL
        if (rng.rand() < 0.5)
            shift_head = true;
        else
            shift_head = false;
        }
    else if (head_idx!=-1){ // only head present
        shift_head = true;
    }
    else{ // only tail present
        shift_head = false;
    }
    
    // Save the kink index of the end that will be shifted
    if (shift_head){worm_end_idx=head_idx;}
    else {worm_end_idx=tail_idx;}
    
    // Extract worm end attributes
    tau = paths[worm_end_idx].tau;
    n = paths[worm_end_idx].n;
    prev = paths[worm_end_idx].prev;
    next = paths[worm_end_idx].next;
    
    // Diagonal energy difference in simplified form
    dV=U*(n-!shift_head)-mu;
    
    // To make acceptance ratio unity,shift tail needs to sample w/ dV=eps-eps_w
    if (!shift_head){dV *= -1;} // dV=eps-eps_w
        
    // Determine the lower and upper bounds of the worm end to be timeshifted
    if (next==-1)
        tau_next = beta;
    else
        tau_next = paths[next].tau;
    tau_prev = paths[prev].tau;
    
    // Sample the new time of the worm end from truncated exponential dist.
    /*:::::::::::::::::::: Truncated Exponential RVS :::::::::::::::::::::::::*/
    Z = 1.0 - exp(-dV*(tau_next-tau_prev));
    tau_new = tau_prev - log(1.0-Z*rng.rand())  / dV;
    /*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
    
//    cout << tau << " " << tau_new << endl;

    
    // Add to PROPOSAL counter
    if (shift_head){
        if (tau_new > tau){advance_head_attempts+=1;}
        else{recede_head_attempts+=1;}
    }
    else{ // shift tail
        if (tau_new > tau){advance_tail_attempts+=1;}
        else{recede_tail_attempts+=1;}
    }
    
    // Determine the length of path to be modified
    l_path = tau_new - tau;
    
    // Determine the total particle change based on wormend to be shifted
    if (shift_head)
        dN = +1.0 * l_path/beta;
    else // shift tail
        dN = -1.0 * l_path/beta;
    
    // Canonical simulations: Restrict updates to interval N:(N-1,N+1)
    if (canonical)
        if ((N_tracker+dN) < (N-1) || (N_tracker+dN) > (N+1)){return;}
    
    // Build the Metropolis condition (R)
    R = 1.0; // Sampling worm end time from truncated exponential makes R unity.

    // Metropolis sampling
    if (rng.rand() < R){
        
        // Add to ACCEPTANCE counter
        if (shift_head){
            if (tau_new > tau){advance_head_accepts+=1;}
            else{recede_head_accepts+=1;}
        }
        else{ // shift tail
            if (tau_new > tau){advance_tail_accepts+=1;}
            else{recede_tail_accepts+=1;}
        }
        
        // Modify the worm end time
        paths[worm_end_idx].tau = tau_new;
        
        // Modify total particle number tracker
        N_tracker += dN;
        
        return;
    }
    else // Reject
        return;
}

/*--------------------------------------------------------------------*/

void insert_kink_before_head(vector<Kink> &paths, int &num_kinks,
                int &head_idx,int &tail_idx,
                int M, int N, double U, double mu, double t,
                vector<vector<int> > &adjacency_matrix, int total_nn,
                double beta, double eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                unsigned long long int &ikbh_attempts,
                unsigned long long int &ikbh_accepts,
                RNG &rng){
    
//    if (t==0.0){return;}
    // Variable declarations
    int prev,i,j,n_i,n_wi,n_j,n_wj,prev_i,prev_j,next_i,next_j,
    src_replica,dest_replica;
    double tau,tau_h,p_site,W,R,p_dkbh,p_ikbh,tau_prev_i,tau_prev_j,
    tau_kink,tau_min,dV_i,dV_j;
        
    // Update only possible if worm head present
    if (head_idx==-1){return;}
    
    // Need at least two sites to perform a spaceshift
    if (M<2){return;}
        
    // Add to proposal counter
    ikbh_attempts += 1;
    
    // Extract the worm head site and replica
    i = paths[head_idx].src;
    src_replica = paths[head_idx].src_replica;
    dest_replica = paths[head_idx].dest_replica;
    
    // Randomly choose a nearest neighbor site
    //boost::random::uniform_int_distribution<> random_nn(0, total_nn-1);
    j = adjacency_matrix[i][rng.randInt(total_nn-1)];
    p_site = 1.0/total_nn;
    
    // Retrieve the time of the worm head
    tau_h = paths[head_idx].tau;
    
    // Determine index of lower/upper kinks of flat where head is (site i)
    prev_i = paths[head_idx].prev;
    next_i = paths[head_idx].next;
    
    // Determine index of lower/upper kinks of flat where head jumps to (site j)
    tau = 0.0;            // tau_prev_j candidate
    prev = j;           // prev_j candidate
    prev_j = j;         // this avoids "variable maybe not initialized" warning
    while (tau<tau_h){
        // Set the lower bound index
        prev_j = prev;
        
        // Update lower bound index and tau candidates for next iteration
        prev = paths[prev].next;
        if (prev==-1){break;}
        tau = paths[prev].tau;
    }
    next_j=prev;
    
    // Determine upper,lower bound times on both sites (upper time not needed)
    tau_prev_i = paths[prev_i].tau;
    tau_prev_j = paths[prev_j].tau;
    
    // Determine lowest time at which kink could've been inserted
    if (tau_prev_i>tau_prev_j){tau_min=tau_prev_i;}
    else {tau_min=tau_prev_j;}
    
    // Randomly choose the time of the kink
    //boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    tau_kink = tau_min + rng.rand()*(tau_h-tau_min);
    if (tau_kink == tau_min){return;}
    
    // Extract no. of particles in the flats adjacent to the new kink
    n_wi = paths[prev_i].n;
    n_i = n_wi-1;
    n_j = paths[prev_j].n;
    n_wj = n_j+1;                   // "w": segment with the extra particle
        
    // Calculate the diagonal energy difference on both sites
    dV_i = (U/2.0)*(n_wi*(n_wi-1)-n_i*(n_i-1)) - mu*(n_wi-n_i);
    dV_j = (U/2.0)*(n_wj*(n_wj-1)-n_j*(n_j-1)) - mu*(n_wj-n_j);
    
    // Calculate the weight ratio W'/W
    W = t * n_wj * exp((dV_i-dV_j)*(tau_h-tau_kink));
    
    // Build the Metropolis ratio (R)
    p_dkbh = 0.5;
    p_ikbh = 0.5;
    R = W * (p_dkbh/p_ikbh) * (tau_h-tau_min)/p_site;
    
    // Metropolis Sampling
    if (rng.rand() < R){ // Accept
        
        // Add to acceptance counter
        ikbh_accepts += 1;
                
        // Change kink that stored head information to a regular kink
        paths[head_idx].tau = tau_kink;
        paths[head_idx].n = n_i;
        paths[head_idx].src = i;  // more like site
        paths[head_idx].dest = j; // more like connecting site
        paths[head_idx].prev = prev_i;
        paths[head_idx].next = next_i;
        
        // Create the kinks on the destination site
        paths[num_kinks]=Kink(tau_kink,n_wj,j,i,prev_j,num_kinks+1,
                              src_replica,dest_replica);
        paths[num_kinks+1]=Kink(tau_h,n_j,j,j,num_kinks,next_j,
                                src_replica,dest_replica);
        
        // Set new worm head index
        head_idx = num_kinks+1;
                
        // "Connect" next of lower bound kink to new kink
        paths[prev_j].next = num_kinks;
        
        // "Connect" prev of next kink to worm head
        if(next_j!=-1){paths[next_j].prev = head_idx;}
                
        // Update number of kinks tracker
        num_kinks += 2;
        
        // If worm head is last kink on site j, update last kinks tracker vector
        if (next_j==-1){last_kinks[j]=head_idx;}
        
        return;
        
    }
    else // Reject
        return;
    }

/*--------------------------------------------------------------------*/

void insert_kink_before_head_2(vector<Kink> &paths, int &num_kinks,
                int &head_idx,int &tail_idx,
                int M, int N, double U, double mu, double t,
                vector<vector<int> > &adjacency_matrix, int total_nn,
                double beta, double eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                unsigned long long int &ikbh_attempts,
                unsigned long long int &ikbh_accepts,
                RNG &rng){
    
//    if (t==0.0){return;}
    // Variable declarations
    int prev,i,j,n_i,n_wi,n_j,n_wj,prev_i,prev_j,next_i,next_j,
    src_replica,dest_replica;
    double tau,tau_h,p_site,R,p_dkbh,p_ikbh,tau_prev_i,tau_prev_j,
    tau_kink,tau_min,dV_i,dV_j,dV;
        
    // Update only possible if worm head present
    if (head_idx==-1){return;}
    
    // Need at least two sites to perform a spaceshift
    if (M<2){return;}
        
    // Add to proposal counter
    ikbh_attempts += 1;
    
    // Extract the worm head site and replica
    i = paths[head_idx].src;
    src_replica = paths[head_idx].src_replica;
    dest_replica = paths[head_idx].dest_replica;
    
    // Randomly choose a nearest neighbor site
    //boost::random::uniform_int_distribution<> random_nn(0, total_nn-1);
    j = adjacency_matrix[i][rng.randInt(total_nn-1)];
    p_site = 1.0/total_nn;
    
    // Retrieve the time of the worm head
    tau_h = paths[head_idx].tau;
    
    // Determine index of lower/upper kinks of flat where head is (site i)
    prev_i = paths[head_idx].prev;
    next_i = paths[head_idx].next;
    
    // Determine index of lower/upper kinks of flat where head jumps to (site j)
    tau = 0.0;            // tau_prev_j candidate
    prev = j;           // prev_j candidate
    prev_j = j;         // this avoids "variable maybe not initialized" warning
    while (tau<tau_h){
        // Set the lower bound index
        prev_j = prev;
        
        // Update lower bound index and tau candidates for next iteration
        prev = paths[prev].next;
        if (prev==-1){break;}
        tau = paths[prev].tau;
    }
    next_j=prev;
    
    // Determine upper,lower bound times on both sites (upper time not needed)
    tau_prev_i = paths[prev_i].tau;
    tau_prev_j = paths[prev_j].tau;
    
    // Determine lowest time at which kink could've been inserted
    if (tau_prev_i>tau_prev_j){tau_min=tau_prev_i;}
    else {tau_min=tau_prev_j;}
    
    // Randomly choose the time of the kink
    //boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    
    // Extract no. of particles in the flats adjacent to the new kink
    n_wi = paths[prev_i].n;
    n_i = n_wi-1;
    n_j = paths[prev_j].n;
    n_wj = n_j+1;                   // "w": segment with the extra particle
        
    // Calculate the diagonal energy difference on both sites
    dV_i = (U/2.0)*(n_wi*(n_wi-1)-n_i*(n_i-1)) - mu*(n_wi-n_i);
    dV_j = (U/2.0)*(n_wj*(n_wj-1)-n_j*(n_j-1)) - mu*(n_wj-n_j);
    dV = dV_j - dV_i;
    if (dV == 0){dV = 0.0;}

    /* :::::::::::::::::::::::::: Truncated Sampling :::::::::::::::::::::::: */
    // Sample time on flat interval from truncated exponential for insertion
    double x,Z,a,b,c;

    a = tau_min;
    b = tau_h;

    c = dV; 

    x = rng.rand();
    Z = 1.0 - exp(-c*(b-a));
    if (dV != 0)
        tau_kink = b + log(1.0-Z*x)  / c;
    else // dV == 0
        tau_kink = b + x*(a-b); // L'hopitale was used
    // if (!is_worm){cout << tau_new << endl;}
    /* :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */
    
    // Calculate the weight ratio W'/W
    // W = t * n_wj * exp((dV_i-dV_j)*(tau_h-tau_kink));
    
    // Build the Metropolis ratio (R)
    p_dkbh = 0.5;
    p_ikbh = 0.5;
    // R = W * (p_dkbh/p_ikbh) * (tau_h-tau_min)/p_site;
    if (dV != 0)
        R = t * n_wj * (p_dkbh/p_ikbh) * (Z/dV) / p_site;
    else
        R = t * n_wj * (p_dkbh/p_ikbh) * (b-a) / p_site; // L'hopitale

    // cout << a << " " << b << " " << c << " " << dV << " " << Z << " " << R << " " << tau_kink << endl;
    
    // Metropolis Sampling
    if (rng.rand() < R){ // Accept
        
        // Add to acceptance counter
        ikbh_accepts += 1;
                
        // Change kink that stored head information to a regular kink
        paths[head_idx].tau = tau_kink;
        paths[head_idx].n = n_i;
        paths[head_idx].src = i;  // more like site
        paths[head_idx].dest = j; // more like connecting site
        paths[head_idx].prev = prev_i;
        paths[head_idx].next = next_i;
        
        // Create the kinks on the destination site
        paths[num_kinks]=Kink(tau_kink,n_wj,j,i,prev_j,num_kinks+1,
                              src_replica,dest_replica);
        paths[num_kinks+1]=Kink(tau_h,n_j,j,j,num_kinks,next_j,
                                src_replica,dest_replica);
        
        // Set new worm head index
        head_idx = num_kinks+1;
                
        // "Connect" next of lower bound kink to new kink
        paths[prev_j].next = num_kinks;
        
        // "Connect" prev of next kink to worm head
        if(next_j!=-1){paths[next_j].prev = head_idx;}
                
        // Update number of kinks tracker
        num_kinks += 2;
        
        // If worm head is last kink on site j, update last kinks tracker vector
        if (next_j==-1){last_kinks[j]=head_idx;}
        
        return;
        
    }
    else // Reject
        return;
    }

/*--------------------------------------------------------------------*/

void delete_kink_before_head(vector<Kink> &paths, int &num_kinks,
                int &head_idx,int &tail_idx,
                int M, int N, double U, double mu, double t,
                vector<vector<int> > &adjacency_matrix, int total_nn,
                double beta, double eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                unsigned long long int &dkbh_attempts,
                unsigned long long int &dkbh_accepts,
                RNG &rng){

    // Variable declarations
    int prev,i,j,n_i,n_wi,n_j,n_wj,prev_i,prev_j,next_i,next_j,
    kink_idx_i,kink_idx_j,src_replica,dest_replica;
    double tau,tau_h,p_site,W,R,p_dkbh,p_ikbh,tau_prev_i,tau_prev_j,
    tau_kink,tau_min,dV_i,dV_j,tau_next_i;

    // Update only possible if worm head present
    if (head_idx==-1){return;}

    // There has to be a regular kink before the worm head
    if (paths[paths[head_idx].prev].src
    ==paths[paths[head_idx].prev].dest){return;}

    // Need at least two sites to perform a spaceshift
    if (M<2){return;}

    // Indices of: upper bound kink, kink before head, lower bound kink ; site j
    next_j = paths[head_idx].next;
    kink_idx_j = paths[head_idx].prev;
    prev_j = paths[kink_idx_j].prev;

    // Times of: worm head, kink before head, lower bound kink; site j
    tau_h = paths[head_idx].tau;
    tau_kink = paths[kink_idx_j].tau;
    tau_prev_j = paths[prev_j].tau;

    // Only kinks in which the particle hops from i TO j can be deleted
    if (paths[kink_idx_j].n-paths[prev_j].n<0){return;}

    // Retrieve worm head site (j) and connecting site (i)
    j = paths[kink_idx_j].src;
    i = paths[kink_idx_j].dest;
    src_replica = paths[kink_idx_j].src_replica;
    dest_replica = paths[kink_idx_j].dest_replica;
    
    // Determine index of lower/upper bounds of flat where kink connects to (i)
    tau = 0.0;            // tau_prev_i candidate
    prev = i;           // prev_i candidate
    prev_i = i;         // this avoids "variable maybe not initialized" warning
    while (tau<tau_kink){
        // Set the lower bound index
        prev_i = prev;

        // Update lower bound index and tau candidates for next iteration
        prev = paths[prev].next;
        if (prev==-1){break;}
        tau = paths[prev].tau;
    }
    kink_idx_i = prev;
    next_i=paths[kink_idx_i].next;

    // Retrieve time of lower,upper bounds on connecting site (i)
    tau_prev_i = paths[prev_i].tau;
    if (next_i!=-1){tau_next_i = paths[next_i].tau;}
    else{tau_next_i=beta;}

    // Deletion cannot interfere w/ kinks on other site
    if (tau_h >= tau_next_i){return;}

    // Add to proposal counter
    dkbh_attempts += 1;

    // Determine lowest time at which kink could've been inserted
    if (tau_prev_i>tau_prev_j){tau_min=tau_prev_i;}
    else {tau_min=tau_prev_j;}

    // Probability of inverse move (ikbh) of choosing site where worm end is
    p_site = 1.0/total_nn;

    // Extract no. of particles in the flats adjacent to the new kink
    n_wi = paths[prev_i].n;
    n_i = n_wi-1;
    n_j = paths[prev_j].n;
    n_wj = n_j+1;                   // "w": segment with the extra particle

    // Calculate the diagonal energy difference on both sites
    dV_i = (U/2.0)*(n_wi*(n_wi-1)-n_i*(n_i-1)) - mu*(n_wi-n_i);
    dV_j = (U/2.0)*(n_wj*(n_wj-1)-n_j*(n_j-1)) - mu*(n_wj-n_j);

    // Calculate the weight ratio W'/W
    W = t * n_wj * exp((dV_i-dV_j)*(tau_h-tau_kink));

    // Build the Metropolis ratio (R)
    p_dkbh = 0.5;
    p_ikbh = 0.5;
    R = W * (p_dkbh/p_ikbh) * (tau_h-tau_min)/p_site;
    R = 1.0/R;

    // Metropolis Sampling
    //boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    if (rng.rand() < R){ // Accept

        // Add to acceptance counter
        dkbh_accepts += 1;

        // Stage 1: Delete kink on i
        if (paths[num_kinks-1].next!=-1)
            paths[paths[num_kinks-1].next].prev = kink_idx_i;
        paths[paths[num_kinks-1].prev].next = kink_idx_i;

        swap(paths[kink_idx_i],paths[num_kinks-1]);

        // Important kinks might've been at end of paths vector
        if (prev_i==num_kinks-1){prev_i=kink_idx_i;}
        else if (next_i==num_kinks-1){next_i=kink_idx_i;}
        else if (prev_j==num_kinks-1){prev_j=kink_idx_i;}
        else if (next_j==num_kinks-1){next_j=kink_idx_i;}
        else if (kink_idx_j==num_kinks-1){kink_idx_j=kink_idx_i;}
        else if (head_idx==num_kinks-1){head_idx=kink_idx_i;}
        else {;}

        // I don't remember why I left this statement out of the block above :(
        if (tail_idx==num_kinks-1){tail_idx=kink_idx_i;}

        // The kink sent to where deleted kink was might be last on it's site
        if (paths[kink_idx_i].next==-1){
            last_kinks[paths[kink_idx_i].src]=kink_idx_i;
        }
        
        // Reconnect upper and lower bounds of the flat
        if (next_i!=-1)
            paths[next_i].prev = prev_i;
        paths[prev_i].next = next_i;

        if (next_i==-1){last_kinks[i]=prev_i;}

        // Stage 2: Delete worm head on j
        if (paths[num_kinks-2].next!=-1)
            paths[paths[num_kinks-2].next].prev = head_idx;
        paths[paths[num_kinks-2].prev].next = head_idx;

        swap(paths[head_idx],paths[num_kinks-2]);

        if (prev_i==num_kinks-2){prev_i=head_idx;}
        else if (next_i==num_kinks-2){next_i=head_idx;}
        else if (prev_j==num_kinks-2){prev_j=head_idx;}
        else if (next_j==num_kinks-2){next_j=head_idx;}
        else if (kink_idx_j==num_kinks-2){kink_idx_j=head_idx;}
        else {;}

        if (tail_idx==num_kinks-2){tail_idx=head_idx;}

        if (paths[head_idx].next==-1){
            last_kinks[paths[head_idx].src]=head_idx;
        }

        if (next_j!=-1)
            paths[next_j].prev = kink_idx_j;
        paths[kink_idx_j].next = next_j;

        if (next_j==-1){last_kinks[j]=kink_idx_j;}

        // Stage 3: Delete kink on j
        if (paths[num_kinks-3].next!=-1)
            paths[paths[num_kinks-3].next].prev = kink_idx_j;
        paths[paths[num_kinks-3].prev].next = kink_idx_j;

        swap(paths[kink_idx_j],paths[num_kinks-3]);

        if (prev_i==num_kinks-3){prev_i=kink_idx_j;}
        else if (next_i==num_kinks-3){next_i=kink_idx_j;}
        else if (prev_j==num_kinks-3){prev_j=kink_idx_j;}
        else if (next_j==num_kinks-3){next_j=kink_idx_j;}
        else {;}

        if (tail_idx==num_kinks-3){tail_idx=kink_idx_j;}

        if (paths[kink_idx_j].next==-1){
            last_kinks[paths[kink_idx_j].src]=kink_idx_j;
        }

        if (next_j!=-1)
            paths[next_j].prev = prev_j;
        paths[prev_j].next = next_j;

        if (next_j==-1){last_kinks[j]=prev_j;}

        // Stage 4: Insert worm head on i
        paths[num_kinks-3]=Kink(tau_h,n_i,i,i,prev_i,next_i,
                                src_replica,dest_replica);

        head_idx = num_kinks-3;

        paths[prev_i].next = head_idx;
        if(next_i!=-1){paths[next_i].prev = head_idx;}

        if (next_i==-1){last_kinks[i]=head_idx;}

        // Update number of kinks tracker
        num_kinks -= 2;

        return;

    }
    else // Reject
        return;
    }

/*--------------------------------------------------------------------*/

void delete_kink_before_head_2(vector<Kink> &paths, int &num_kinks,
                int &head_idx,int &tail_idx,
                int M, int N, double U, double mu, double t,
                vector<vector<int> > &adjacency_matrix, int total_nn,
                double beta, double eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                unsigned long long int &dkbh_attempts,
                unsigned long long int &dkbh_accepts,
                RNG &rng){

    // Variable declarations
    int prev,i,j,n_i,n_wi,n_j,n_wj,prev_i,prev_j,next_i,next_j,
    kink_idx_i,kink_idx_j,src_replica,dest_replica;
    double tau,tau_h,p_site,R,p_dkbh,p_ikbh,tau_prev_i,tau_prev_j,
    tau_kink,tau_min,dV_i,dV_j,tau_next_i,dV;

    // Update only possible if worm head present
    if (head_idx==-1){return;}

    // There has to be a regular kink before the worm head
    if (paths[paths[head_idx].prev].src
    ==paths[paths[head_idx].prev].dest){return;}

    // Need at least two sites to perform a spaceshift
    if (M<2){return;}

    // Indices of: upper bound kink, kink before head, lower bound kink ; site j
    next_j = paths[head_idx].next;
    kink_idx_j = paths[head_idx].prev;
    prev_j = paths[kink_idx_j].prev;

    // Times of: worm head, kink before head, lower bound kink; site j
    tau_h = paths[head_idx].tau;
    tau_kink = paths[kink_idx_j].tau;
    tau_prev_j = paths[prev_j].tau;

    // Only kinks in which the particle hops from i TO j can be deleted
    if (paths[kink_idx_j].n-paths[prev_j].n<0){return;}

    // Retrieve worm head site (j) and connecting site (i)
    j = paths[kink_idx_j].src;
    i = paths[kink_idx_j].dest;
    src_replica = paths[kink_idx_j].src_replica;
    dest_replica = paths[kink_idx_j].dest_replica;
    
    // Determine index of lower/upper bounds of flat where kink connects to (i)
    tau = 0.0;            // tau_prev_i candidate
    prev = i;           // prev_i candidate
    prev_i = i;         // this avoids "variable maybe not initialized" warning
    while (tau<tau_kink){
        // Set the lower bound index
        prev_i = prev;

        // Update lower bound index and tau candidates for next iteration
        prev = paths[prev].next;
        if (prev==-1){break;}
        tau = paths[prev].tau;
    }
    kink_idx_i = prev;
    next_i=paths[kink_idx_i].next;

    // Retrieve time of lower,upper bounds on connecting site (i)
    tau_prev_i = paths[prev_i].tau;
    if (next_i!=-1){tau_next_i = paths[next_i].tau;}
    else{tau_next_i=beta;}

    // Deletion cannot interfere w/ kinks on other site
    if (tau_h >= tau_next_i){return;}

    // Add to proposal counter
    dkbh_attempts += 1;

    // Determine lowest time at which kink could've been inserted
    if (tau_prev_i>tau_prev_j){tau_min=tau_prev_i;}
    else {tau_min=tau_prev_j;}

    // Probability of inverse move (ikbh) of choosing site where worm end is
    p_site = 1.0/total_nn;

    // Extract no. of particles in the flats adjacent to the new kink
    n_wi = paths[prev_i].n;
    n_i = n_wi-1;
    n_j = paths[prev_j].n;
    n_wj = n_j+1;                   // "w": segment with the extra particle

    // Calculate the diagonal energy difference on both sites
    dV_i = (U/2.0)*(n_wi*(n_wi-1)-n_i*(n_i-1)) - mu*(n_wi-n_i);
    dV_j = (U/2.0)*(n_wj*(n_wj-1)-n_j*(n_j-1)) - mu*(n_wj-n_j);
    dV = dV_j - dV_i;
    // if (dV == 0){cout << "NOOOO" << endl;}

    // Calculate the weight ratio W'/W
    // W = t * n_wj * exp((dV_i-dV_j)*(tau_h-tau_kink));

    // inverse move (insert kink before head) truncated sampling
    double Z,a,b,c;
    a = tau_min;
    b = tau_h;
    c = dV;
    Z = 1.0 - exp(-c*(b-a)); //

    // Build the Metropolis ratio (R)
    p_dkbh = 0.5;
    p_ikbh = 0.5;
    if (dV != 0)
        R = t * n_wj * (p_dkbh/p_ikbh) * (Z/dV) / p_site;
    else
        R = t * n_wj * (p_dkbh/p_ikbh) * (b-a) / p_site;
    // R = W * (p_dkbh/p_ikbh) * (tau_h-tau_min)/p_site;
    R = 1.0/R;

    // Metropolis Sampling
    //boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    if (rng.rand() < R){ // Accept

        // Add to acceptance counter
        dkbh_accepts += 1;

        // Stage 1: Delete kink on i
        if (paths[num_kinks-1].next!=-1)
            paths[paths[num_kinks-1].next].prev = kink_idx_i;
        paths[paths[num_kinks-1].prev].next = kink_idx_i;

        swap(paths[kink_idx_i],paths[num_kinks-1]);

        // Important kinks might've been at end of paths vector
        if (prev_i==num_kinks-1){prev_i=kink_idx_i;}
        else if (next_i==num_kinks-1){next_i=kink_idx_i;}
        else if (prev_j==num_kinks-1){prev_j=kink_idx_i;}
        else if (next_j==num_kinks-1){next_j=kink_idx_i;}
        else if (kink_idx_j==num_kinks-1){kink_idx_j=kink_idx_i;}
        else if (head_idx==num_kinks-1){head_idx=kink_idx_i;}
        else {;}

        // I don't remember why I left this statement out of the block above :(
        if (tail_idx==num_kinks-1){tail_idx=kink_idx_i;}

        // The kink sent to where deleted kink was might be last on it's site
        if (paths[kink_idx_i].next==-1){
            last_kinks[paths[kink_idx_i].src]=kink_idx_i;
        }
        
        // Reconnect upper and lower bounds of the flat
        if (next_i!=-1)
            paths[next_i].prev = prev_i;
        paths[prev_i].next = next_i;

        if (next_i==-1){last_kinks[i]=prev_i;}

        // Stage 2: Delete worm head on j
        if (paths[num_kinks-2].next!=-1)
            paths[paths[num_kinks-2].next].prev = head_idx;
        paths[paths[num_kinks-2].prev].next = head_idx;

        swap(paths[head_idx],paths[num_kinks-2]);

        if (prev_i==num_kinks-2){prev_i=head_idx;}
        else if (next_i==num_kinks-2){next_i=head_idx;}
        else if (prev_j==num_kinks-2){prev_j=head_idx;}
        else if (next_j==num_kinks-2){next_j=head_idx;}
        else if (kink_idx_j==num_kinks-2){kink_idx_j=head_idx;}
        else {;}

        if (tail_idx==num_kinks-2){tail_idx=head_idx;}

        if (paths[head_idx].next==-1){
            last_kinks[paths[head_idx].src]=head_idx;
        }

        if (next_j!=-1)
            paths[next_j].prev = kink_idx_j;
        paths[kink_idx_j].next = next_j;

        if (next_j==-1){last_kinks[j]=kink_idx_j;}

        // Stage 3: Delete kink on j
        if (paths[num_kinks-3].next!=-1)
            paths[paths[num_kinks-3].next].prev = kink_idx_j;
        paths[paths[num_kinks-3].prev].next = kink_idx_j;

        swap(paths[kink_idx_j],paths[num_kinks-3]);

        if (prev_i==num_kinks-3){prev_i=kink_idx_j;}
        else if (next_i==num_kinks-3){next_i=kink_idx_j;}
        else if (prev_j==num_kinks-3){prev_j=kink_idx_j;}
        else if (next_j==num_kinks-3){next_j=kink_idx_j;}
        else {;}

        if (tail_idx==num_kinks-3){tail_idx=kink_idx_j;}

        if (paths[kink_idx_j].next==-1){
            last_kinks[paths[kink_idx_j].src]=kink_idx_j;
        }

        if (next_j!=-1)
            paths[next_j].prev = prev_j;
        paths[prev_j].next = next_j;

        if (next_j==-1){last_kinks[j]=prev_j;}

        // Stage 4: Insert worm head on i
        paths[num_kinks-3]=Kink(tau_h,n_i,i,i,prev_i,next_i,
                                src_replica,dest_replica);

        head_idx = num_kinks-3;

        paths[prev_i].next = head_idx;
        if(next_i!=-1){paths[next_i].prev = head_idx;}

        if (next_i==-1){last_kinks[i]=head_idx;}

        // Update number of kinks tracker
        num_kinks -= 2;

        return;

    }
    else // Reject
        return;
    }

/*--------------------------------------------------------------------*/

void insert_kink_after_head(vector<Kink> &paths, int &num_kinks,
                int &head_idx,int &tail_idx,
                int M, int N, double U, double mu, double t,
                vector<vector<int> > &adjacency_matrix, int total_nn,
                double beta, double eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                unsigned long long int &ikah_attempts,
                unsigned long long int &ikah_accepts,
                RNG &rng){
    
    // Variable declarations
    int prev,i,j,n_i,n_wi,n_j,n_wj,prev_i,prev_j,next_i,next_j,
    src_replica,dest_replica;
    double tau,tau_h,p_site,W,R,p_dkah,p_ikah,
    tau_kink,tau_max,dV_i,dV_j,tau_next_i,tau_next_j;
    
    // Update only possible if worm head present
    if (head_idx==-1){return;}
    
    // Need at least two sites to perform a spaceshift
    if (M<2){return;}
    
    // Add to proposal counter
    ikah_attempts += 1;
    
    // Extract the worm head site
    i = paths[head_idx].src;
    src_replica = paths[head_idx].src_replica;
    dest_replica = paths[head_idx].dest_replica;

    // Randomly choose a nearest neighbor site
    //boost::random::uniform_int_distribution<> random_nn(0, total_nn-1);
    j = adjacency_matrix[i][rng.randInt(total_nn-1)];
    p_site = 1.0/total_nn;
    
    // Retrieve the time of the worm head
    tau_h = paths[head_idx].tau;
    
    // Determine index of lower/upper kinks of flat where head is (site i)
    prev_i = paths[head_idx].prev;
    next_i = paths[head_idx].next;
    
    // Determine index of lower/upper kinks of flat where head jumps to (site j)
    tau = 0.0;            // tau_prev_j candidate
    prev = j;           // prev_j candidate
    prev_j = j;         // this avoids "variable maybe not initialized" warning
    while (tau<tau_h){

        // Set the lower bound index
        prev_j = prev;
        
        // Update lower bound index and tau candidates for next iteration
        prev = paths[prev].next;
        if (prev==-1){break;}
        tau = paths[prev].tau;
    }
    next_j=prev;
    
    // Determine upper,lower bound times on both sites
    if (next_i!=-1)
        tau_next_i = paths[next_i].tau;
    else
        tau_next_i = beta;
    if (next_j!=-1)
        tau_next_j = paths[next_j].tau;
    else
        tau_next_j = beta;
    
    // Determine highest time at which kink could've been inserted
    if (tau_next_i<tau_next_j){tau_max=tau_next_i;}
    else {tau_max=tau_next_j;}
    
    // Randomly choose the time of the kink
    //boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    tau_kink = tau_h + rng.rand()*(tau_max-tau_h);
    if (tau_kink==tau_h){return;}
    
     // Extract no. of particles in the flats adjacent to the new kink
     n_wi = paths[prev_i].n;
     n_i = n_wi-1;
     n_wj = paths[prev_j].n;
     n_j = n_wj-1;                   // "w": segment with the extra particle
    
    // Update not possible if no particles on destinaton site (j)
    if (n_wj==0){return;}

    // Calculate the diagonal energy difference on both sites
    dV_i = (U/2.0)*(n_wi*(n_wi-1)-n_i*(n_i-1)) - mu*(n_wi-n_i);
    dV_j = (U/2.0)*(n_wj*(n_wj-1)-n_j*(n_j-1)) - mu*(n_wj-n_j);
    
    // Calculate the weight ratio W'/W
    W = t * n_wj * exp((-dV_i+dV_j)*(tau_kink-tau_h));
    
    // Build the Metropolis ratio (R)
    p_dkah = 0.5;
    p_ikah = 0.5;
    R = W * (p_dkah/p_ikah) * (tau_max-tau_h)/p_site;
    
    // Metropolis Sampling
    if (rng.rand() < R){ // Accept
        
        // Add to acceptance counter
        ikah_accepts += 1;
                
        // Change kink that stored head information to a regular kink
        paths[head_idx].tau = tau_kink;
        paths[head_idx].n = n_i;
        paths[head_idx].src = i;
        paths[head_idx].dest = j;
        paths[head_idx].prev = prev_i;
        paths[head_idx].next = next_i;
        
        // Create the kinks on the destination site
        paths[num_kinks]=Kink(tau_h,n_j,j,j,prev_j,num_kinks+1,
                              src_replica,dest_replica);
        paths[num_kinks+1]=Kink(tau_kink,n_wj,j,i,num_kinks,next_j,
                                src_replica,dest_replica);
        
        // Set new worm head index
        head_idx = num_kinks;
                
        // "Connect" next of lower bound kink to worm head
        paths[prev_j].next = head_idx;
        
        // "Connect" prev of next kink to new kink
        if(next_j!=-1){paths[next_j].prev = num_kinks+1;}
        
        // If new kink is last kink on site j, update last kinks tracker vector
        if (next_j==-1){last_kinks[j]=num_kinks+1;}
        
        // Update number of kinks tracker
        num_kinks += 2;

        return;
        }
        else // Reject
            return;
        }

/*--------------------------------------------------------------------*/

void insert_kink_after_head_2(vector<Kink> &paths, int &num_kinks,
                int &head_idx,int &tail_idx,
                int M, int N, double U, double mu, double t,
                vector<vector<int> > &adjacency_matrix, int total_nn,
                double beta, double eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                unsigned long long int &ikah_attempts,
                unsigned long long int &ikah_accepts,
                RNG &rng){
    
    // Variable declarations
    int prev,i,j,n_i,n_wi,n_j,n_wj,prev_i,prev_j,next_i,next_j,
    src_replica,dest_replica;
    double tau,tau_h,p_site,R,p_dkah,p_ikah,
    tau_kink,tau_max,dV_i,dV_j,tau_next_i,tau_next_j,dV;
    
    // Update only possible if worm head present
    if (head_idx==-1){return;}
    
    // Need at least two sites to perform a spaceshift
    if (M<2){return;}
    
    // Add to proposal counter
    ikah_attempts += 1;
    
    // Extract the worm head site
    i = paths[head_idx].src;
    src_replica = paths[head_idx].src_replica;
    dest_replica = paths[head_idx].dest_replica;

    // Randomly choose a nearest neighbor site
    //boost::random::uniform_int_distribution<> random_nn(0, total_nn-1);
    j = adjacency_matrix[i][rng.randInt(total_nn-1)];
    p_site = 1.0/total_nn;
    
    // Retrieve the time of the worm head
    tau_h = paths[head_idx].tau;
    
    // Determine index of lower/upper kinks of flat where head is (site i)
    prev_i = paths[head_idx].prev;
    next_i = paths[head_idx].next;
    
    // Determine index of lower/upper kinks of flat where head jumps to (site j)
    tau = 0.0;            // tau_prev_j candidate
    prev = j;           // prev_j candidate
    prev_j = j;         // this avoids "variable maybe not initialized" warning
    while (tau<tau_h){

        // Set the lower bound index
        prev_j = prev;
        
        // Update lower bound index and tau candidates for next iteration
        prev = paths[prev].next;
        if (prev==-1){break;}
        tau = paths[prev].tau;
    }
    next_j=prev;
    
    // Determine upper,lower bound times on both sites
    if (next_i!=-1)
        tau_next_i = paths[next_i].tau;
    else
        tau_next_i = beta;
    if (next_j!=-1)
        tau_next_j = paths[next_j].tau;
    else
        tau_next_j = beta;
    
    // Determine highest time at which kink could've been inserted
    if (tau_next_i<tau_next_j){tau_max=tau_next_i;}
    else {tau_max=tau_next_j;}
    
    // Randomly choose the time of the kink
    //boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    // tau_kink = tau_h + rng.rand()*(tau_max-tau_h);
    // if (tau_kink==tau_h){return;}
    
     // Extract no. of particles in the flats adjacent to the new kink
     n_wi = paths[prev_i].n;
     n_i = n_wi-1;
     n_wj = paths[prev_j].n;
     n_j = n_wj-1;                   // "w": segment with the extra particle
    
    // Update not possible if no particles on destinaton site (j)
    if (n_wj==0){return;}

    // Calculate the diagonal energy difference on both sites
    dV_i = (U/2.0)*(n_wi*(n_wi-1)-n_i*(n_i-1)) - mu*(n_wi-n_i);
    dV_j = (U/2.0)*(n_wj*(n_wj-1)-n_j*(n_j-1)) - mu*(n_wj-n_j);
    dV = dV_i - dV_j;
    if (dV == 0){dV = 0.0;}

    /* :::::::::::::::::::::::::: Truncated Sampling :::::::::::::::::::::::: */
    // Sample time on flat interval from truncated exponential for insertion
    double x,Z,a,b,c;

    a = tau_h;
    b = tau_max;

    c = dV; 

    x = rng.rand();
    Z = 1.0 - exp(-c*(b-a));
    if (dV != 0)
        tau_kink = a - log(1.0-Z*x)  / c;
    else // dV == 0
        tau_kink = a - x*(a-b); // L'Hopitale
    // if (!is_worm){cout << tau_new << endl;}
    /* :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */
    
    // Calculate the weight ratio W'/W
    // W = t * n_wj * exp((-dV_i+dV_j)*(tau_kink-tau_h));
    
    // Build the Metropolis ratio (R)
    p_dkah = 0.5;
    p_ikah = 0.5;
    if (dV != 0)
        R = t * n_wj * (p_dkah/p_ikah) * (Z/dV) / p_site;
    else
        R = t * n_wj * (p_dkah/p_ikah) * (b-a) / p_site; // L'Hopitale


    // cout << a << " " << b << " " << c << " " << dV << " " << Z << " " << R << " " << tau_kink << endl;
    
    // Metropolis Sampling
    if (rng.rand() < R){ // Accept
        
        // Add to acceptance counter
        ikah_accepts += 1;
                
        // Change kink that stored head information to a regular kink
        paths[head_idx].tau = tau_kink;
        paths[head_idx].n = n_i;
        paths[head_idx].src = i;
        paths[head_idx].dest = j;
        paths[head_idx].prev = prev_i;
        paths[head_idx].next = next_i;
        
        // Create the kinks on the destination site
        paths[num_kinks]=Kink(tau_h,n_j,j,j,prev_j,num_kinks+1,
                              src_replica,dest_replica);
        paths[num_kinks+1]=Kink(tau_kink,n_wj,j,i,num_kinks,next_j,
                                src_replica,dest_replica);
        
        // Set new worm head index
        head_idx = num_kinks;
                
        // "Connect" next of lower bound kink to worm head
        paths[prev_j].next = head_idx;
        
        // "Connect" prev of next kink to new kink
        if(next_j!=-1){paths[next_j].prev = num_kinks+1;}
        
        // If new kink is last kink on site j, update last kinks tracker vector
        if (next_j==-1){last_kinks[j]=num_kinks+1;}
        
        // Update number of kinks tracker
        num_kinks += 2;

        return;
        }
        else // Reject
            return;
        }

/*--------------------------------------------------------------------*/

void delete_kink_after_head(vector<Kink> &paths, int &num_kinks,
                int &head_idx,int &tail_idx,
                int M, int N, double U, double mu, double t,
                vector<vector<int> > &adjacency_matrix, int total_nn,
                double beta, double eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                unsigned long long int &dkah_attempts,
                unsigned long long int &dkah_accepts,
                RNG &rng){
    
    // Variable declarations
    int prev,i,j,n_i,n_wi,n_j,n_wj,prev_i,prev_j,next_i,next_j,
    kink_idx_i,kink_idx_j,src_replica,dest_replica;
    double tau,tau_h,p_site,W,R,p_dkah,p_ikah,tau_prev_i,
    tau_kink,tau_max,dV_i,dV_j,tau_next_i,tau_next_j;
    
    // Update only possible if worm head present
    if (head_idx==-1){return;}
    
    // There has to be a regular kink after the worm head
    if (paths[head_idx].next==tail_idx ||
        paths[head_idx].next==-1){return;}

    // Need at least two sites to perform a spaceshift
    if (M<2){return;}

    // Indices of: upper bound kink, kink before head, lower bound kink ; site j
    kink_idx_j = paths[head_idx].next;
    next_j = paths[kink_idx_j].next;
    prev_j = paths[head_idx].prev;

    // Times of: worm head, kink before head, lower bound kink; site j
    if (next_j!=-1)
        tau_next_j = paths[next_j].tau;
    else
        tau_next_j = beta;
    tau_kink = paths[kink_idx_j].tau;
    tau_h = paths[head_idx].tau;
    
    // Only kinks in which the particle hops from i TO j can be deleted
    if (paths[kink_idx_j].n-paths[head_idx].n<0){return;}
    
    // Retrieve worm head site (j) and connecting site (i)
    j = paths[head_idx].src;
    i = paths[kink_idx_j].dest;
    src_replica = paths[kink_idx_j].src_replica;
    dest_replica = paths[kink_idx_j].dest_replica;

    // Determine index of lower/upper bounds of flat where kink connects to (i)
    tau = 0.0;            // tau_prev_i candidate
    prev = i;           // prev_i candidate
    prev_i = i;         // this avoids "variable maybe not initialized" warning
    while (tau<tau_kink){
        // Set the lower bound index
        prev_i = prev;

        // Update lower bound index and tau candidates for next iteration
        prev = paths[prev].next;
        if (prev==-1){break;}
        tau = paths[prev].tau;
    }
    kink_idx_i = prev;
    next_i=paths[kink_idx_i].next;

    // Retrieve time of lower,upper bounds on connecting site (i)
    tau_prev_i = paths[prev_i].tau;
    if (next_i!=-1){tau_next_i = paths[next_i].tau;}
    else{tau_next_i=beta;}
    
    // Deletion cannot interfere w/ kinks on other site
    if (tau_h <= tau_prev_i){return;}

    // Add to proposal counter
    dkah_attempts += 1;

    // Determine highest time at which kink could've been inserted
    if (tau_next_i<tau_next_j){tau_max=tau_next_i;}
    else {tau_max=tau_next_j;}

    // Probability of inverse move (ikah) choosing site where worm end is
    p_site = 1.0/total_nn;

    // Extract no. of particles in the flats adjacent to the new kink
    n_wi = paths[prev_i].n;
    n_i = n_wi-1;
    n_wj = paths[prev_j].n;
    n_j = n_wj-1;                   // "w": segment with the extra particle

    // Calculate the diagonal energy difference on both sites
    dV_i = (U/2.0)*(n_wi*(n_wi-1)-n_i*(n_i-1)) - mu*(n_wi-n_i);
    dV_j = (U/2.0)*(n_wj*(n_wj-1)-n_j*(n_j-1)) - mu*(n_wj-n_j);

    // Calculate the weight ratio W'/W
    W = t * n_wj * exp((-dV_i+dV_j)*(tau_kink-tau_h));

    // Build the Metropolis ratio (R)
    p_dkah = 0.5;
    p_ikah = 0.5;
    R = W * (p_dkah/p_ikah) * (tau_max-tau_h)/p_site;
    R = 1.0/R;

    // Metropolis Sampling
    //boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    if (rng.rand() < R){ // Accept

        // Add to acceptance counter
        dkah_accepts += 1;
        
        // Stage 1: Delete kink on i
        if (paths[num_kinks-1].next!=-1)
            paths[paths[num_kinks-1].next].prev = kink_idx_i;
        paths[paths[num_kinks-1].prev].next = kink_idx_i;
        
        swap(paths[kink_idx_i],paths[num_kinks-1]);
        
        if (prev_i==num_kinks-1){prev_i=kink_idx_i;}
        else if (next_i==num_kinks-1){next_i=kink_idx_i;}
        else if (prev_j==num_kinks-1){prev_j=kink_idx_i;}
        else if (next_j==num_kinks-1){next_j=kink_idx_i;}
        else if (kink_idx_j==num_kinks-1){kink_idx_j=kink_idx_i;}
        else if (head_idx==num_kinks-1){head_idx=kink_idx_i;}
        else {;}
        
        if (tail_idx==num_kinks-1){tail_idx=kink_idx_i;}
        
        if (paths[kink_idx_i].next==-1){
            last_kinks[paths[kink_idx_i].src]=kink_idx_i;
        }
        
        if (next_i!=-1)
            paths[next_i].prev = prev_i;
        paths[prev_i].next = next_i;
        
        if (next_i==-1){last_kinks[i]=prev_i;}

        // Stage 2: Delete kink on j
        if (paths[num_kinks-2].next!=-1)
            paths[paths[num_kinks-2].next].prev = kink_idx_j;
        paths[paths[num_kinks-2].prev].next = kink_idx_j;
        
        swap(paths[kink_idx_j],paths[num_kinks-2]);
        
        if (prev_i==num_kinks-2){prev_i=kink_idx_j;}
        else if (next_i==num_kinks-2){next_i=kink_idx_j;}
        else if (prev_j==num_kinks-2){prev_j=kink_idx_j;}
        else if (next_j==num_kinks-2){next_j=kink_idx_j;}
        else if (head_idx==num_kinks-2){head_idx=kink_idx_j;}
        else {;}
        
        if (tail_idx==num_kinks-2){tail_idx=kink_idx_j;}
        
        if (paths[kink_idx_j].next==-1){
            last_kinks[paths[kink_idx_j].src]=kink_idx_j;
        }
        
        if (next_j!=-1)
            paths[next_j].prev = head_idx;
        paths[head_idx].next = next_j;
        
        if (next_j==-1){last_kinks[j]=head_idx;}
        
        // Stage 3: Delete worm head on j
        if (paths[num_kinks-3].next!=-1)
            paths[paths[num_kinks-3].next].prev = head_idx;
        paths[paths[num_kinks-3].prev].next = head_idx;
        
        swap(paths[head_idx],paths[num_kinks-3]);
        
        if (prev_i==num_kinks-3){prev_i=head_idx;}
        else if (next_i==num_kinks-3){next_i=head_idx;}
        else if (prev_j==num_kinks-3){prev_j=head_idx;}
        else if (next_j==num_kinks-3){next_j=head_idx;}
        else {;}
        
        if (tail_idx==num_kinks-3){tail_idx=head_idx;}
        
        if (paths[head_idx].next==-1){
            last_kinks[paths[head_idx].src]=head_idx;
        }
        
        if (next_j!=-1)
            paths[next_j].prev = prev_j;
        paths[prev_j].next = next_j;
        
        if (next_j==-1){last_kinks[j]=prev_j;}
        
        // Stage 4: Insert worm head on i
        paths[num_kinks-3]=Kink(tau_h,n_i,i,i,prev_i,next_i,
                                src_replica,dest_replica);
        
        head_idx = num_kinks-3;
        
        paths[prev_i].next = head_idx;
        if(next_i!=-1){paths[next_i].prev = head_idx;}
        
        if (next_i==-1){last_kinks[i]=head_idx;}
        
        // Update number of kinks tracker
        num_kinks -= 2;

        return;

    }
    else // Reject
        return;
    }

/*--------------------------------------------------------------------*/

void delete_kink_after_head_2(vector<Kink> &paths, int &num_kinks,
                int &head_idx,int &tail_idx,
                int M, int N, double U, double mu, double t,
                vector<vector<int> > &adjacency_matrix, int total_nn,
                double beta, double eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                unsigned long long int &dkah_attempts,
                unsigned long long int &dkah_accepts,
                RNG &rng){
    
    // Variable declarations
    int prev,i,j,n_i,n_wi,n_j,n_wj,prev_i,prev_j,next_i,next_j,
    kink_idx_i,kink_idx_j,src_replica,dest_replica;
    double tau,tau_h,p_site,R,p_dkah,p_ikah,tau_prev_i,
    tau_kink,tau_max,dV_i,dV_j,tau_next_i,tau_next_j,dV;
    
    // Update only possible if worm head present
    if (head_idx==-1){return;}
    
    // There has to be a regular kink after the worm head
    if (paths[head_idx].next==tail_idx ||
        paths[head_idx].next==-1){return;}

    // Need at least two sites to perform a spaceshift
    if (M<2){return;}

    // Indices of: upper bound kink, kink before head, lower bound kink ; site j
    kink_idx_j = paths[head_idx].next;
    next_j = paths[kink_idx_j].next;
    prev_j = paths[head_idx].prev;

    // Times of: worm head, kink before head, lower bound kink; site j
    if (next_j!=-1)
        tau_next_j = paths[next_j].tau;
    else
        tau_next_j = beta;
    tau_kink = paths[kink_idx_j].tau;
    tau_h = paths[head_idx].tau;
    
    // Only kinks in which the particle hops from i TO j can be deleted
    if (paths[kink_idx_j].n-paths[head_idx].n<0){return;}
    
    // Retrieve worm head site (j) and connecting site (i)
    j = paths[head_idx].src;
    i = paths[kink_idx_j].dest;
    src_replica = paths[kink_idx_j].src_replica;
    dest_replica = paths[kink_idx_j].dest_replica;

    // Determine index of lower/upper bounds of flat where kink connects to (i)
    tau = 0.0;            // tau_prev_i candidate
    prev = i;           // prev_i candidate
    prev_i = i;         // this avoids "variable maybe not initialized" warning
    while (tau<tau_kink){
        // Set the lower bound index
        prev_i = prev;

        // Update lower bound index and tau candidates for next iteration
        prev = paths[prev].next;
        if (prev==-1){break;}
        tau = paths[prev].tau;
    }
    kink_idx_i = prev;
    next_i=paths[kink_idx_i].next;

    // Retrieve time of lower,upper bounds on connecting site (i)
    tau_prev_i = paths[prev_i].tau;
    if (next_i!=-1){tau_next_i = paths[next_i].tau;}
    else{tau_next_i=beta;}
    
    // Deletion cannot interfere w/ kinks on other site
    if (tau_h <= tau_prev_i){return;}

    // Add to proposal counter
    dkah_attempts += 1;

    // Determine highest time at which kink could've been inserted
    if (tau_next_i<tau_next_j){tau_max=tau_next_i;}
    else {tau_max=tau_next_j;}

    // Probability of inverse move (ikah) choosing site where worm end is
    p_site = 1.0/total_nn;

    // Extract no. of particles in the flats adjacent to the new kink
    n_wi = paths[prev_i].n;
    n_i = n_wi-1;
    n_wj = paths[prev_j].n;
    n_j = n_wj-1;                   // "w": segment with the extra particle

    // Calculate the diagonal energy difference on both sites
    dV_i = (U/2.0)*(n_wi*(n_wi-1)-n_i*(n_i-1)) - mu*(n_wi-n_i);
    dV_j = (U/2.0)*(n_wj*(n_wj-1)-n_j*(n_j-1)) - mu*(n_wj-n_j);
    dV = dV_i - dV_j;
    // if (dV == 0){cout << "NOOOO" << endl;}

    // Calculate the weight ratio W'/W
    // W = t * n_wj * exp((-dV_i+dV_j)*(tau_kink-tau_h));

    // inverse move (insert kink after head) tuncated sampling
    double Z,a,b,c;
    a = tau_h;
    b = tau_max;
    c = dV;
    Z = 1.0 - exp(-c*(b-a)); 

    // Build the Metropolis ratio (R)
    p_dkah = 0.5;
    p_ikah = 0.5;
    if (dV != 0)
        R = t * n_wj * (p_dkah/p_ikah) * (Z/dV) / p_site;
    else
        R = t * n_wj * (p_dkah/p_ikah) * (b-a) / p_site;
    // R = W * (p_dkah/p_ikah) * (tau_max-tau_h)/p_site;
    R = 1.0/R;

    // Metropolis Sampling
    //boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    if (rng.rand() < R){ // Accept

        // Add to acceptance counter
        dkah_accepts += 1;
        
        // Stage 1: Delete kink on i
        if (paths[num_kinks-1].next!=-1)
            paths[paths[num_kinks-1].next].prev = kink_idx_i;
        paths[paths[num_kinks-1].prev].next = kink_idx_i;
        
        swap(paths[kink_idx_i],paths[num_kinks-1]);
        
        if (prev_i==num_kinks-1){prev_i=kink_idx_i;}
        else if (next_i==num_kinks-1){next_i=kink_idx_i;}
        else if (prev_j==num_kinks-1){prev_j=kink_idx_i;}
        else if (next_j==num_kinks-1){next_j=kink_idx_i;}
        else if (kink_idx_j==num_kinks-1){kink_idx_j=kink_idx_i;}
        else if (head_idx==num_kinks-1){head_idx=kink_idx_i;}
        else {;}
        
        if (tail_idx==num_kinks-1){tail_idx=kink_idx_i;}
        
        if (paths[kink_idx_i].next==-1){
            last_kinks[paths[kink_idx_i].src]=kink_idx_i;
        }
        
        if (next_i!=-1)
            paths[next_i].prev = prev_i;
        paths[prev_i].next = next_i;
        
        if (next_i==-1){last_kinks[i]=prev_i;}

        // Stage 2: Delete kink on j
        if (paths[num_kinks-2].next!=-1)
            paths[paths[num_kinks-2].next].prev = kink_idx_j;
        paths[paths[num_kinks-2].prev].next = kink_idx_j;
        
        swap(paths[kink_idx_j],paths[num_kinks-2]);
        
        if (prev_i==num_kinks-2){prev_i=kink_idx_j;}
        else if (next_i==num_kinks-2){next_i=kink_idx_j;}
        else if (prev_j==num_kinks-2){prev_j=kink_idx_j;}
        else if (next_j==num_kinks-2){next_j=kink_idx_j;}
        else if (head_idx==num_kinks-2){head_idx=kink_idx_j;}
        else {;}
        
        if (tail_idx==num_kinks-2){tail_idx=kink_idx_j;}
        
        if (paths[kink_idx_j].next==-1){
            last_kinks[paths[kink_idx_j].src]=kink_idx_j;
        }
        
        if (next_j!=-1)
            paths[next_j].prev = head_idx;
        paths[head_idx].next = next_j;
        
        if (next_j==-1){last_kinks[j]=head_idx;}
        
        // Stage 3: Delete worm head on j
        if (paths[num_kinks-3].next!=-1)
            paths[paths[num_kinks-3].next].prev = head_idx;
        paths[paths[num_kinks-3].prev].next = head_idx;
        
        swap(paths[head_idx],paths[num_kinks-3]);
        
        if (prev_i==num_kinks-3){prev_i=head_idx;}
        else if (next_i==num_kinks-3){next_i=head_idx;}
        else if (prev_j==num_kinks-3){prev_j=head_idx;}
        else if (next_j==num_kinks-3){next_j=head_idx;}
        else {;}
        
        if (tail_idx==num_kinks-3){tail_idx=head_idx;}
        
        if (paths[head_idx].next==-1){
            last_kinks[paths[head_idx].src]=head_idx;
        }
        
        if (next_j!=-1)
            paths[next_j].prev = prev_j;
        paths[prev_j].next = next_j;
        
        if (next_j==-1){last_kinks[j]=prev_j;}
        
        // Stage 4: Insert worm head on i
        paths[num_kinks-3]=Kink(tau_h,n_i,i,i,prev_i,next_i,
                                src_replica,dest_replica);
        
        head_idx = num_kinks-3;
        
        paths[prev_i].next = head_idx;
        if(next_i!=-1){paths[next_i].prev = head_idx;}
        
        if (next_i==-1){last_kinks[i]=head_idx;}
        
        // Update number of kinks tracker
        num_kinks -= 2;

        return;

    }
    else // Reject
        return;
    }

/*--------------------------------------------------------------------*/

void insert_kink_before_tail(vector<Kink> &paths, int &num_kinks,
                int &head_idx,int &tail_idx,
                int M, int N, double U, double mu, double t,
                vector<vector<int> > &adjacency_matrix, int total_nn,
                double beta, double eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                unsigned long long int &ikbt_attempts,
                unsigned long long int &ikbt_accepts,
                RNG &rng){
    
    // Variable declarations
    int prev,i,j,n_i,n_wi,n_j,n_wj,prev_i,prev_j,next_i,next_j,
    src_replica,dest_replica;
    double tau,tau_t,p_site,W,R,p_dkbt,p_ikbt,tau_prev_i,tau_prev_j,
    tau_kink,tau_min,dV_i,dV_j;
    
    // Update only possible if worm tail present
    if (tail_idx==-1){return;}
    
    // Need at least two sites to perform a spaceshift
    if (M<2){return;}
    
    // Extract the worm tail site
    i = paths[tail_idx].src;
    src_replica = paths[tail_idx].src_replica;
    dest_replica = paths[tail_idx].dest_replica;
    
    // Randomly choose a nearest neighbor site
    //boost::random::uniform_int_distribution<> random_nn(0, total_nn-1);
    j = adjacency_matrix[i][rng.randInt(total_nn-1)];
    p_site = 1.0/total_nn;
    
    // Retrieve the time of the worm tail
    tau_t = paths[tail_idx].tau;
    
    // Determine index of lower/upper kinks of flat where tail is (site i)
    prev_i = paths[tail_idx].prev;
    next_i = paths[tail_idx].next;
    
    // Determine index of lower/upper kinks of flat where tail jumps to (site j)
    tau = 0.0;            // tau_prev_j candidate
    prev = j;           // prev_j candidate
    prev_j = j;         // this avoids "variable maybe not initialized" warning
    while (tau<tau_t){
        // Set the lower bound index
        prev_j = prev;
        
        // Update lower bound index and tau candidates for next iteration
        prev = paths[prev].next;
        if (prev==-1){break;}
        tau = paths[prev].tau;
    }
    next_j=prev;
    
    // Determine upper,lower bound times on both sites
    tau_prev_i = paths[prev_i].tau;
    tau_prev_j = paths[prev_j].tau;
    
    // Determine lowest time at which kink could've been inserted
    if (tau_prev_i>tau_prev_j){tau_min=tau_prev_i;}
    else {tau_min=tau_prev_j;}
    
    // Randomly choose the time of the kink
    //boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    tau_kink = tau_min + rng.rand()*(tau_t-tau_min);
    if (tau_kink == tau_min){return;}
    
     // Extract no. of particles in the flats adjacent to the new kink
     n_i = paths[prev_i].n;
     n_wi = n_i+1;
     n_wj = paths[prev_j].n;
     n_j = n_wj-1;                   // "w": segment with the extra particle
    
    // Update not possible if no particles on destinaton site (j)
    if (n_wj == 0){return;}
    
    // Add to proposal counter
    ikbt_attempts += 1;
    
    // Calculate the diagonal energy difference on both sites
    dV_i = (U/2.0)*(n_wi*(n_wi-1)-n_i*(n_i-1)) - mu*(n_wi-n_i);
    dV_j = (U/2.0)*(n_wj*(n_wj-1)-n_j*(n_j-1)) - mu*(n_wj-n_j);
    
    // Calculate the weight ratio W'/W
    W = t * n_wj * exp((-dV_i+dV_j)*(tau_t-tau_kink));

    // Build the Metropolis ratio (R)
    p_dkbt = 0.5;
    p_ikbt = 0.5;
    R = W * (p_dkbt/p_ikbt) * (tau_t-tau_min)/p_site;
    
    // Metropolis Sampling
    if (rng.rand() < R){ // Accept
        
        // Add to acceptance counter
        ikbt_accepts += 1;
                
        // Change kink that stored tail information to a regular kink
        paths[tail_idx].tau = tau_kink;
        paths[tail_idx].n = n_wi;
        paths[tail_idx].src = i;
        paths[tail_idx].dest = j;
        paths[tail_idx].prev = prev_i;
        paths[tail_idx].next = next_i;
        
        // Create the kinks on the destination site
        paths[num_kinks]=Kink(tau_kink,n_j,j,i,prev_j,num_kinks+1,
                              src_replica,dest_replica);
        paths[num_kinks+1]=Kink(tau_t,n_wj,j,j,num_kinks,next_j,
                              src_replica,dest_replica);
        
        // Set new worm tail index
        tail_idx = num_kinks+1;
                
        // "Connect" next of lower bound kink to new kink
        paths[prev_j].next = num_kinks;
        
        // "Connect" prev of next kink to worm tail
        if(next_j!=-1){paths[next_j].prev = tail_idx;}
                
        // Update number of kinks tracker
        num_kinks += 2;
        
        // If worm tail is last kink on site j, update last kinks tracker vector
        if (next_j==-1){last_kinks[j]=tail_idx;}
        
        return;
            
        }
        else // Reject
            return;
        }

/*--------------------------------------------------------------------*/

void insert_kink_before_tail_2(vector<Kink> &paths, int &num_kinks,
                int &head_idx,int &tail_idx,
                int M, int N, double U, double mu, double t,
                vector<vector<int> > &adjacency_matrix, int total_nn,
                double beta, double eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                unsigned long long int &ikbt_attempts,
                unsigned long long int &ikbt_accepts,
                RNG &rng){
    
    // Variable declarations
    int prev,i,j,n_i,n_wi,n_j,n_wj,prev_i,prev_j,next_i,next_j,
    src_replica,dest_replica;
    double tau,tau_t,p_site,R,p_dkbt,p_ikbt,tau_prev_i,tau_prev_j,
    tau_kink,tau_min,dV_i,dV_j,dV;
    
    // Update only possible if worm tail present
    if (tail_idx==-1){return;}
    
    // Need at least two sites to perform a spaceshift
    if (M<2){return;}
    
    // Extract the worm tail site
    i = paths[tail_idx].src;
    src_replica = paths[tail_idx].src_replica;
    dest_replica = paths[tail_idx].dest_replica;
    
    // Randomly choose a nearest neighbor site
    //boost::random::uniform_int_distribution<> random_nn(0, total_nn-1);
    j = adjacency_matrix[i][rng.randInt(total_nn-1)];
    p_site = 1.0/total_nn;
    
    // Retrieve the time of the worm tail
    tau_t = paths[tail_idx].tau;
    
    // Determine index of lower/upper kinks of flat where tail is (site i)
    prev_i = paths[tail_idx].prev;
    next_i = paths[tail_idx].next;
    
    // Determine index of lower/upper kinks of flat where tail jumps to (site j)
    tau = 0.0;            // tau_prev_j candidate
    prev = j;           // prev_j candidate
    prev_j = j;         // this avoids "variable maybe not initialized" warning
    while (tau<tau_t){
        // Set the lower bound index
        prev_j = prev;
        
        // Update lower bound index and tau candidates for next iteration
        prev = paths[prev].next;
        if (prev==-1){break;}
        tau = paths[prev].tau;
    }
    next_j=prev;
    
    // Determine upper,lower bound times on both sites
    tau_prev_i = paths[prev_i].tau;
    tau_prev_j = paths[prev_j].tau;
    
    // Determine lowest time at which kink could've been inserted
    if (tau_prev_i>tau_prev_j){tau_min=tau_prev_i;}
    else {tau_min=tau_prev_j;}
    
    // Randomly choose the time of the kink
    //boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    // tau_kink = tau_min + rng.rand()*(tau_t-tau_min);
    // if (tau_kink == tau_min){return;}
    
     // Extract no. of particles in the flats adjacent to the new kink
     n_i = paths[prev_i].n;
     n_wi = n_i+1;
     n_wj = paths[prev_j].n;
     n_j = n_wj-1;                   // "w": segment with the extra particle
    
    // Update not possible if no particles on destinaton site (j)
    if (n_wj == 0){return;}
    
    // Add to proposal counter
    ikbt_attempts += 1;
    
    // Calculate the diagonal energy difference on both sites
    dV_i = (U/2.0)*(n_wi*(n_wi-1)-n_i*(n_i-1)) - mu*(n_wi-n_i);
    dV_j = (U/2.0)*(n_wj*(n_wj-1)-n_j*(n_j-1)) - mu*(n_wj-n_j);
    dV = dV_i - dV_j;
    if (dV == 0){dV = 0.0;}

    /* :::::::::::::::::::::::::: Truncated Sampling :::::::::::::::::::::::: */
    // Sample time on flat interval from truncated exponential for insertion
    double x,Z,a,b,c;

    a = tau_min;
    b = tau_t;

    c = dV; 

    x = rng.rand();
    Z = 1.0 - exp(-c*(b-a)); //
    if (dV != 0)
        tau_kink = b + log(1.0-Z*x)  / c;
    else // dV == 0
        tau_kink = b + x*(a-b); // L'hopitale was used
    /* :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */
    
    // Calculate the weight ratio W'/W
    // W = t * n_wj * exp((-dV_i+dV_j)*(tau_t-tau_kink));

    // Build the Metropolis ratio (R)
    p_dkbt = 0.5;
    p_ikbt = 0.5;
    if (dV != 0)
        R = t * n_wj * (p_dkbt/p_ikbt) * (Z/dV) /p_site;
    else
        R = t * n_wj * (p_dkbt/p_ikbt) * (b-a) /p_site;
    // R = W * (p_dkbt/p_ikbt) * (tau_t-tau_min)/p_site;

    // cout << a << " " << b << " " << c << " " << dV << " " << Z << " " << R << " " << tau_kink << endl;

    
    // Metropolis Sampling
    if (rng.rand() < R){ // Accept
        
        // Add to acceptance counter
        ikbt_accepts += 1;
                
        // Change kink that stored tail information to a regular kink
        paths[tail_idx].tau = tau_kink;
        paths[tail_idx].n = n_wi;
        paths[tail_idx].src = i;
        paths[tail_idx].dest = j;
        paths[tail_idx].prev = prev_i;
        paths[tail_idx].next = next_i;
        
        // Create the kinks on the destination site
        paths[num_kinks]=Kink(tau_kink,n_j,j,i,prev_j,num_kinks+1,
                              src_replica,dest_replica);
        paths[num_kinks+1]=Kink(tau_t,n_wj,j,j,num_kinks,next_j,
                              src_replica,dest_replica);
        
        // Set new worm tail index
        tail_idx = num_kinks+1;
                
        // "Connect" next of lower bound kink to new kink
        paths[prev_j].next = num_kinks;
        
        // "Connect" prev of next kink to worm tail
        if(next_j!=-1){paths[next_j].prev = tail_idx;}
                
        // Update number of kinks tracker
        num_kinks += 2;
        
        // If worm tail is last kink on site j, update last kinks tracker vector
        if (next_j==-1){last_kinks[j]=tail_idx;}
        
        return;
            
        }
        else // Reject
            return;
        }

/*--------------------------------------------------------------------*/

void delete_kink_before_tail(vector<Kink> &paths, int &num_kinks,
                int &head_idx,int &tail_idx,
                int M, int N, double U, double mu, double t,
                vector<vector<int> > &adjacency_matrix, int total_nn,
                double beta, double eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                unsigned long long int &dkbt_attempts,
                unsigned long long int &dkbt_accepts,
                RNG &rng){
    
    // Variable declarations
    int prev,i,j,n_i,n_wi,n_j,n_wj,prev_i,prev_j,next_i,next_j,
    kink_idx_i,kink_idx_j,src_replica,dest_replica;
    double tau,tau_t,p_site,W,R,p_dkbt,p_ikbt,tau_prev_i,tau_prev_j,
    tau_kink,tau_min,dV_i,dV_j,tau_next_i;
    
    // Update only possible if worm tail present
    if (tail_idx==-1){return;}
    
    // There has to be a regular kink after the worm tail
    if (paths[tail_idx].prev==head_idx ||
        paths[paths[tail_idx].prev].tau==0){return;}

    // Need at least two sites to perform a spaceshift
    if (M<2){return;}
    
    // Indices of: upper bound kink, kink before tail, lower bound kink ; site j
    next_j = paths[tail_idx].next;
    kink_idx_j = paths[tail_idx].prev;
    prev_j = paths[kink_idx_j].prev;

    // Times of: worm tail, kink before tail, lower bound kink; site j
    tau_t = paths[tail_idx].tau;
    tau_kink = paths[kink_idx_j].tau;
    tau_prev_j = paths[prev_j].tau;
    
    // Only kinks in which the particle hops from j TO i can be deleted
    if (paths[kink_idx_j].n-paths[prev_j].n>0){return;}
    
    // Retrieve worm tail site (j) and connecting site (i)
    j = paths[kink_idx_j].src;
    i = paths[kink_idx_j].dest;
    src_replica = paths[kink_idx_j].src_replica;
    dest_replica = paths[kink_idx_j].dest_replica;
    
    // Determine index of lower/upper bounds of flat where kink connects to (i)
    tau = 0.0;            // tau_prev_i candidate
    prev = i;           // prev_i candidate
    prev_i = i;         // this avoids "variable maybe not initialized" warning
    while (tau<tau_kink){
        // Set the lower bound index
        prev_i = prev;

        // Update lower bound index and tau candidates for next iteration
        prev = paths[prev].next;
        if (prev==-1){break;}
        tau = paths[prev].tau;
    }
    kink_idx_i = prev;
    next_i=paths[kink_idx_i].next;

    // Retrieve time of lower,upper bounds on connecting site (i)
    tau_prev_i = paths[prev_i].tau;
    if (next_i==-1){tau_next_i=beta;}
    else {tau_next_i=paths[next_i].tau;};
    
    // Deletion cannot interfere w/ kinks on other site
    if (tau_t>=tau_next_i){return;}

    // Add to proposal counter
    dkbt_attempts += 1;

    // Determine lowest time at which kink could've been inserted
    if (tau_prev_i>tau_prev_j){tau_min=tau_prev_i;}
    else {tau_min=tau_prev_j;}
    
    // Probability of inverse move (ikbt) choosing site where worm end is
    p_site = 1.0/total_nn;

    // Extract no. of particles in the flats adjacent to the new kink
    n_i = paths[prev_i].n;
    n_wi = n_i+1;
    n_wj = paths[prev_j].n;
    n_j = n_wj-1;                   // "w": segment with the extra particle

    // Calculate the diagonal energy difference on both sites
    dV_i = (U/2.0)*(n_wi*(n_wi-1)-n_i*(n_i-1)) - mu*(n_wi-n_i);
    dV_j = (U/2.0)*(n_wj*(n_wj-1)-n_j*(n_j-1)) - mu*(n_wj-n_j);

    // Calculate the weight ratio W'/W
    W = t * n_wj * exp((-dV_i+dV_j)*(tau_t-tau_kink));

    // Build the Metropolis ratio (R)
    p_dkbt = 0.5;
    p_ikbt = 0.5;
    R = W * (p_dkbt/p_ikbt) * (tau_t-tau_min)/p_site;
    R = 1.0/R;
    
    // Metropolis Sampling
    //boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    if (rng.rand() < R){ // Accept

        // Add to acceptance counter
        dkbt_accepts += 1;
        
        // Stage 1: Delete kink on i
        if (paths[num_kinks-1].next!=-1)
            paths[paths[num_kinks-1].next].prev = kink_idx_i;
        paths[paths[num_kinks-1].prev].next = kink_idx_i;
        
        swap(paths[kink_idx_i],paths[num_kinks-1]);
        
        if (prev_i==num_kinks-1){prev_i=kink_idx_i;}
        else if (next_i==num_kinks-1){next_i=kink_idx_i;}
        else if (prev_j==num_kinks-1){prev_j=kink_idx_i;}
        else if (next_j==num_kinks-1){next_j=kink_idx_i;}
        else if (kink_idx_j==num_kinks-1){kink_idx_j=kink_idx_i;}
        else if (tail_idx==num_kinks-1){tail_idx=kink_idx_i;}
        else {;}
        
        if (head_idx==num_kinks-1){head_idx=kink_idx_i;}
        
        if (paths[kink_idx_i].next==-1){
            last_kinks[paths[kink_idx_i].src]=kink_idx_i;
        }
        
        if (next_i!=-1)
            paths[next_i].prev = prev_i;
        paths[prev_i].next = next_i;
        
        if (next_i==-1){last_kinks[i]=prev_i;}

        // Stage 2: Delete worm tail on j
        if (paths[num_kinks-2].next!=-1)
            paths[paths[num_kinks-2].next].prev = tail_idx;
        paths[paths[num_kinks-2].prev].next = tail_idx;
        
        swap(paths[tail_idx],paths[num_kinks-2]);
        
        if (prev_i==num_kinks-2){prev_i=tail_idx;}
        else if (next_i==num_kinks-2){next_i=tail_idx;}
        else if (prev_j==num_kinks-2){prev_j=tail_idx;}
        else if (next_j==num_kinks-2){next_j=tail_idx;}
        else if (kink_idx_j==num_kinks-2){kink_idx_j=tail_idx;}
        else {;}
        
        if (head_idx==num_kinks-2){head_idx=tail_idx;}
        
        if (paths[tail_idx].next==-1){
            last_kinks[paths[tail_idx].src]=tail_idx;
        }
        
        if (next_j!=-1)
            paths[next_j].prev = kink_idx_j;
        paths[kink_idx_j].next = next_j;
        
        if (next_j==-1){last_kinks[j]=kink_idx_j;}
                
        // Stage 3: Delete kink on j
        if (paths[num_kinks-3].next!=-1)
            paths[paths[num_kinks-3].next].prev = kink_idx_j;
        paths[paths[num_kinks-3].prev].next = kink_idx_j;
        
        swap(paths[kink_idx_j],paths[num_kinks-3]);
        
        if (prev_i==num_kinks-3){prev_i=kink_idx_j;}
        else if (next_i==num_kinks-3){next_i=kink_idx_j;}
        else if (prev_j==num_kinks-3){prev_j=kink_idx_j;}
        else if (next_j==num_kinks-3){next_j=kink_idx_j;}
        else {;}
        
        if (head_idx==num_kinks-3){head_idx=kink_idx_j;}
        
        if (paths[kink_idx_j].next==-1){
            last_kinks[paths[kink_idx_j].src]=kink_idx_j;
        }
        
        if (next_j!=-1)
            paths[next_j].prev = prev_j;
        paths[prev_j].next = next_j;
        
        if (next_j==-1){last_kinks[j]=prev_j;}

        // Stage 4: Insert worm tail on i
        paths[num_kinks-3]=Kink(tau_t,n_wi,i,i,prev_i,next_i,
                                src_replica,dest_replica);

        tail_idx = num_kinks-3;

        paths[prev_i].next = tail_idx;
        if(next_i!=-1){paths[next_i].prev = tail_idx;}

        if (next_i==-1){last_kinks[i]=tail_idx;}

        // Update number of kinks tracker
        num_kinks -= 2;
        
        return;

    }
    else // Reject
        return;
    }

/*--------------------------------------------------------------------*/

void delete_kink_before_tail_2(vector<Kink> &paths, int &num_kinks,
                int &head_idx,int &tail_idx,
                int M, int N, double U, double mu, double t,
                vector<vector<int> > &adjacency_matrix, int total_nn,
                double beta, double eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                unsigned long long int &dkbt_attempts,
                unsigned long long int &dkbt_accepts,
                RNG &rng){
    
    // Variable declarations
    int prev,i,j,n_i,n_wi,n_j,n_wj,prev_i,prev_j,next_i,next_j,
    kink_idx_i,kink_idx_j,src_replica,dest_replica;
    double tau,tau_t,p_site,R,p_dkbt,p_ikbt,tau_prev_i,tau_prev_j,
    tau_kink,tau_min,dV_i,dV_j,tau_next_i,dV;
    
    // Update only possible if worm tail present
    if (tail_idx==-1){return;}
    
    // There has to be a regular kink after the worm tail
    if (paths[tail_idx].prev==head_idx ||
        paths[paths[tail_idx].prev].tau==0){return;}

    // Need at least two sites to perform a spaceshift
    if (M<2){return;}
    
    // Indices of: upper bound kink, kink before tail, lower bound kink ; site j
    next_j = paths[tail_idx].next;
    kink_idx_j = paths[tail_idx].prev;
    prev_j = paths[kink_idx_j].prev;

    // Times of: worm tail, kink before tail, lower bound kink; site j
    tau_t = paths[tail_idx].tau;
    tau_kink = paths[kink_idx_j].tau;
    tau_prev_j = paths[prev_j].tau;
    
    // Only kinks in which the particle hops from j TO i can be deleted
    if (paths[kink_idx_j].n-paths[prev_j].n>0){return;}
    
    // Retrieve worm tail site (j) and connecting site (i)
    j = paths[kink_idx_j].src;
    i = paths[kink_idx_j].dest;
    src_replica = paths[kink_idx_j].src_replica;
    dest_replica = paths[kink_idx_j].dest_replica;
    
    // Determine index of lower/upper bounds of flat where kink connects to (i)
    tau = 0.0;            // tau_prev_i candidate
    prev = i;           // prev_i candidate
    prev_i = i;         // this avoids "variable maybe not initialized" warning
    while (tau<tau_kink){
        // Set the lower bound index
        prev_i = prev;

        // Update lower bound index and tau candidates for next iteration
        prev = paths[prev].next;
        if (prev==-1){break;}
        tau = paths[prev].tau;
    }
    kink_idx_i = prev;
    next_i=paths[kink_idx_i].next;

    // Retrieve time of lower,upper bounds on connecting site (i)
    tau_prev_i = paths[prev_i].tau;
    if (next_i==-1){tau_next_i=beta;}
    else {tau_next_i=paths[next_i].tau;};
    
    // Deletion cannot interfere w/ kinks on other site
    if (tau_t>=tau_next_i){return;}

    // Add to proposal counter
    dkbt_attempts += 1;

    // Determine lowest time at which kink could've been inserted
    if (tau_prev_i>tau_prev_j){tau_min=tau_prev_i;}
    else {tau_min=tau_prev_j;}
    
    // Probability of inverse move (ikbt) choosing site where worm end is
    p_site = 1.0/total_nn;

    // Extract no. of particles in the flats adjacent to the new kink
    n_i = paths[prev_i].n;
    n_wi = n_i+1;
    n_wj = paths[prev_j].n;
    n_j = n_wj-1;                   // "w": segment with the extra particle

    // Calculate the diagonal energy difference on both sites
    dV_i = (U/2.0)*(n_wi*(n_wi-1)-n_i*(n_i-1)) - mu*(n_wi-n_i);
    dV_j = (U/2.0)*(n_wj*(n_wj-1)-n_j*(n_j-1)) - mu*(n_wj-n_j);
    dV = dV_i - dV_j;
    // if (dV == 0){cout << "NOOOO" << endl;}

    // Calculate the weight ratio W'/W
    // W = t * n_wj * exp((-dV_i+dV_j)*(tau_t-tau_kink));

    // inverse move (insert kink before tail) truncated exponential sampling
    double Z,a,b,c;
    a = tau_min;
    b = tau_t;
    c = dV;
    Z = 1.0 - exp(-c*(b-a));

    // Build the Metropolis ratio (R)
    p_dkbt = 0.5;
    p_ikbt = 0.5;
    if (dV != 0)
        R = t * n_wj * (p_dkbt/p_ikbt) * (Z/c) /p_site;
    else
        R = t * n_wj * (p_dkbt/p_ikbt) * (b-a) /p_site;
    // R = W * (p_dkbt/p_ikbt) * (tau_t-tau_min)/p_site;
    R = 1.0/R;
    
    // Metropolis Sampling
    //boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    if (rng.rand() < R){ // Accept

        // Add to acceptance counter
        dkbt_accepts += 1;
        
        // Stage 1: Delete kink on i
        if (paths[num_kinks-1].next!=-1)
            paths[paths[num_kinks-1].next].prev = kink_idx_i;
        paths[paths[num_kinks-1].prev].next = kink_idx_i;
        
        swap(paths[kink_idx_i],paths[num_kinks-1]);
        
        if (prev_i==num_kinks-1){prev_i=kink_idx_i;}
        else if (next_i==num_kinks-1){next_i=kink_idx_i;}
        else if (prev_j==num_kinks-1){prev_j=kink_idx_i;}
        else if (next_j==num_kinks-1){next_j=kink_idx_i;}
        else if (kink_idx_j==num_kinks-1){kink_idx_j=kink_idx_i;}
        else if (tail_idx==num_kinks-1){tail_idx=kink_idx_i;}
        else {;}
        
        if (head_idx==num_kinks-1){head_idx=kink_idx_i;}
        
        if (paths[kink_idx_i].next==-1){
            last_kinks[paths[kink_idx_i].src]=kink_idx_i;
        }
        
        if (next_i!=-1)
            paths[next_i].prev = prev_i;
        paths[prev_i].next = next_i;
        
        if (next_i==-1){last_kinks[i]=prev_i;}

        // Stage 2: Delete worm tail on j
        if (paths[num_kinks-2].next!=-1)
            paths[paths[num_kinks-2].next].prev = tail_idx;
        paths[paths[num_kinks-2].prev].next = tail_idx;
        
        swap(paths[tail_idx],paths[num_kinks-2]);
        
        if (prev_i==num_kinks-2){prev_i=tail_idx;}
        else if (next_i==num_kinks-2){next_i=tail_idx;}
        else if (prev_j==num_kinks-2){prev_j=tail_idx;}
        else if (next_j==num_kinks-2){next_j=tail_idx;}
        else if (kink_idx_j==num_kinks-2){kink_idx_j=tail_idx;}
        else {;}
        
        if (head_idx==num_kinks-2){head_idx=tail_idx;}
        
        if (paths[tail_idx].next==-1){
            last_kinks[paths[tail_idx].src]=tail_idx;
        }
        
        if (next_j!=-1)
            paths[next_j].prev = kink_idx_j;
        paths[kink_idx_j].next = next_j;
        
        if (next_j==-1){last_kinks[j]=kink_idx_j;}
                
        // Stage 3: Delete kink on j
        if (paths[num_kinks-3].next!=-1)
            paths[paths[num_kinks-3].next].prev = kink_idx_j;
        paths[paths[num_kinks-3].prev].next = kink_idx_j;
        
        swap(paths[kink_idx_j],paths[num_kinks-3]);
        
        if (prev_i==num_kinks-3){prev_i=kink_idx_j;}
        else if (next_i==num_kinks-3){next_i=kink_idx_j;}
        else if (prev_j==num_kinks-3){prev_j=kink_idx_j;}
        else if (next_j==num_kinks-3){next_j=kink_idx_j;}
        else {;}
        
        if (head_idx==num_kinks-3){head_idx=kink_idx_j;}
        
        if (paths[kink_idx_j].next==-1){
            last_kinks[paths[kink_idx_j].src]=kink_idx_j;
        }
        
        if (next_j!=-1)
            paths[next_j].prev = prev_j;
        paths[prev_j].next = next_j;
        
        if (next_j==-1){last_kinks[j]=prev_j;}

        // Stage 4: Insert worm tail on i
        paths[num_kinks-3]=Kink(tau_t,n_wi,i,i,prev_i,next_i,
                                src_replica,dest_replica);

        tail_idx = num_kinks-3;

        paths[prev_i].next = tail_idx;
        if(next_i!=-1){paths[next_i].prev = tail_idx;}

        if (next_i==-1){last_kinks[i]=tail_idx;}

        // Update number of kinks tracker
        num_kinks -= 2;
        
        return;

    }
    else // Reject
        return;
    }

/*--------------------------------------------------------------------*/

void insert_kink_after_tail(vector<Kink> &paths, int &num_kinks,
                int &head_idx,int &tail_idx,
                int M, int N, double U, double mu, double t,
                vector<vector<int> > &adjacency_matrix, int total_nn,
                double beta, double eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                unsigned long long int &ikat_attempts,
                unsigned long long int &ikat_accepts,
                RNG &rng){
    
    // Variable declarations
    int prev,i,j,n_i,n_wi,n_j,n_wj,prev_i,prev_j,next_i,next_j,
    src_replica,dest_replica;
    double tau,tau_t,p_site,W,R,p_dkat,p_ikat,
    tau_kink,tau_max,dV_i,dV_j,tau_next_i,tau_next_j;
    
    // Update only possible if worm tail present
    if (tail_idx==-1){return;}
    
    // Need at least two sites to perform a spaceshift
    if (M<2){return;}
    
    // Add to proposal counter
    ikat_attempts += 1;
    
    // Extract the worm tail site
    i = paths[tail_idx].src;
    src_replica = paths[tail_idx].src_replica;
    dest_replica = paths[tail_idx].dest_replica;
    
    // Randomly choose a nearest neighbor site
    //boost::random::uniform_int_distribution<> random_nn(0, total_nn-1);
    j = adjacency_matrix[i][rng.randInt(total_nn-1)];
    p_site = 1.0/total_nn;
    
    // Retrieve the time of the worm tail
    tau_t = paths[tail_idx].tau;
    
    // Determine index of lower/upper kinks of flat where tail is (site i)
    prev_i = paths[tail_idx].prev;
    next_i = paths[tail_idx].next;
    
    // Determine index of lower/upper kinks of flat where tail jumps to (site j)
    tau = 0.0;            // tau_prev_j candidate
    prev = j;           // prev_j candidate
    prev_j = j;         // this avoids "variable maybe not initialized" warning
    while (tau<tau_t){
        // Set the lower bound index
        prev_j = prev;
        
        // Update lower bound index and tau candidates for next iteration
        prev = paths[prev].next;
        if (prev==-1){break;}
        tau = paths[prev].tau;
    }
    next_j=prev;
    
    // Determine upper,lower bound times on both sites
    if (next_i!=-1)
        tau_next_i = paths[next_i].tau;
    else
        tau_next_i = beta;
    if (next_j!=-1)
        tau_next_j = paths[next_j].tau;
    else
        tau_next_j = beta;
    
    // Determine highest time at which kink could've been inserted
    if (tau_next_i<tau_next_j){tau_max=tau_next_i;}
    else {tau_max=tau_next_j;}
    
    // Randomly choose the time of the kink
    //boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    tau_kink = tau_t + rng.rand()*(tau_max-tau_t);
    if (tau_kink==tau_t){return;}
    
     // Extract no. of particles in the flats adjacent to the new kink
     n_i = paths[prev_i].n;
     n_wi = n_i+1;
     n_j = paths[prev_j].n;
     n_wj = n_j+1;                   // "w": segment with the extra particle
    
    // Calculate the diagonal energy difference on both sites
    dV_i = (U/2.0)*(n_wi*(n_wi-1)-n_i*(n_i-1)) - mu*(n_wi-n_i);
    dV_j = (U/2.0)*(n_wj*(n_wj-1)-n_j*(n_j-1)) - mu*(n_wj-n_j);
    
    // Calculate the weight ratio W'/W
    W = t * n_wj * exp((dV_i-dV_j)*(tau_kink-tau_t));

    // Build the Metropolis ratio (R)
    p_dkat = 0.5;
    p_ikat = 0.5;
    R = W * (p_dkat/p_ikat) * (tau_max-tau_t)/p_site;
    
    // Metropolis Sampling
    if (rng.rand() < R){ // Accept
        
        // Add to acceptance counter
        ikat_accepts += 1;
                
        // Change kink that stored tail information to a regular kink
        paths[tail_idx].tau = tau_kink;
        paths[tail_idx].n = n_wi;
        paths[tail_idx].src = i;
        paths[tail_idx].dest = j;
        paths[tail_idx].prev = prev_i;
        paths[tail_idx].next = next_i;
        
        // Create the kinks on the destination site
        paths[num_kinks]=Kink(tau_t,n_wj,j,j,prev_j,num_kinks+1,
                              src_replica,dest_replica);
        paths[num_kinks+1]=Kink(tau_kink,n_j,j,i,num_kinks,next_j,
                              src_replica,dest_replica);
        
        // Set new worm tail index
        tail_idx = num_kinks;
                
        // "Connect" next of lower bound kink to worm tail
        paths[prev_j].next = num_kinks;
        
        // "Connect" prev of next kink to new kink
        if(next_j!=-1){paths[next_j].prev = num_kinks+1;}
        
        // If new kink is last kink on site j, update last kinks tracker vector
        if (next_j==-1){last_kinks[j]=num_kinks+1;}
        
        // Update number of kinks tracker
        num_kinks += 2;
        
        return;
            
        }
        else // Reject
            return;
        }

/*--------------------------------------------------------------------*/

void insert_kink_after_tail_2(vector<Kink> &paths, int &num_kinks,
                int &head_idx,int &tail_idx,
                int M, int N, double U, double mu, double t,
                vector<vector<int> > &adjacency_matrix, int total_nn,
                double beta, double eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                unsigned long long int &ikat_attempts,
                unsigned long long int &ikat_accepts,
                RNG &rng){
    
    // Variable declarations
    int prev,i,j,n_i,n_wi,n_j,n_wj,prev_i,prev_j,next_i,next_j,
    src_replica,dest_replica;
    double tau,tau_t,p_site,R,p_dkat,p_ikat,
    tau_kink,tau_max,dV_i,dV_j,tau_next_i,tau_next_j,dV;
    
    // Update only possible if worm tail present
    if (tail_idx==-1){return;}
    
    // Need at least two sites to perform a spaceshift
    if (M<2){return;}
    
    // Add to proposal counter
    ikat_attempts += 1;
    
    // Extract the worm tail site
    i = paths[tail_idx].src;
    src_replica = paths[tail_idx].src_replica;
    dest_replica = paths[tail_idx].dest_replica;
    
    // Randomly choose a nearest neighbor site
    //boost::random::uniform_int_distribution<> random_nn(0, total_nn-1);
    j = adjacency_matrix[i][rng.randInt(total_nn-1)];
    p_site = 1.0/total_nn;
    
    // Retrieve the time of the worm tail
    tau_t = paths[tail_idx].tau;
    
    // Determine index of lower/upper kinks of flat where tail is (site i)
    prev_i = paths[tail_idx].prev;
    next_i = paths[tail_idx].next;
    
    // Determine index of lower/upper kinks of flat where tail jumps to (site j)
    tau = 0.0;            // tau_prev_j candidate
    prev = j;           // prev_j candidate
    prev_j = j;         // this avoids "variable maybe not initialized" warning
    while (tau<tau_t){
        // Set the lower bound index
        prev_j = prev;
        
        // Update lower bound index and tau candidates for next iteration
        prev = paths[prev].next;
        if (prev==-1){break;}
        tau = paths[prev].tau;
    }
    next_j=prev;
    
    // Determine upper,lower bound times on both sites
    if (next_i!=-1)
        tau_next_i = paths[next_i].tau;
    else
        tau_next_i = beta;
    if (next_j!=-1)
        tau_next_j = paths[next_j].tau;
    else
        tau_next_j = beta;
    
    // Determine highest time at which kink could've been inserted
    if (tau_next_i<tau_next_j){tau_max=tau_next_i;}
    else {tau_max=tau_next_j;}
    
    // Randomly choose the time of the kink
    //boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    // tau_kink = tau_t + rng.rand()*(tau_max-tau_t);
    // if (tau_kink==tau_t){return;}
    
     // Extract no. of particles in the flats adjacent to the new kink
     n_i = paths[prev_i].n;
     n_wi = n_i+1;
     n_j = paths[prev_j].n;
     n_wj = n_j+1;                   // "w": segment with the extra particle
    
    // Calculate the diagonal energy difference on both sites
    dV_i = (U/2.0)*(n_wi*(n_wi-1)-n_i*(n_i-1)) - mu*(n_wi-n_i);
    dV_j = (U/2.0)*(n_wj*(n_wj-1)-n_j*(n_j-1)) - mu*(n_wj-n_j);
    dV = dV_j - dV_i;
    if (dV == 0){dV = 0.0;}

    /* :::::::::::::::::::::::::: Truncated Sampling :::::::::::::::::::::::: */
    // Sample time on flat interval from truncated exponential for insertion
    double x,Z,a,b,c;

    a = tau_t;
    b = tau_max;
    c = dV; 

    x = rng.rand();
    Z = 1.0 - exp(-c*(b-a));
    if (dV !=0)
        tau_kink = a - log(1.0-Z*x)  / c;
    else // dV == 0
        tau_kink = a - x*(a-b); // L'Hopitale

    // if (!is_worm){cout << tau_new << endl;}
    /* :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */
    
    // Calculate the weight ratio W'/W
    // W = t * n_wj * exp((dV_i-dV_j)*(tau_kink-tau_t));

    // Build the Metropolis ratio (R)
    p_dkat = 0.5;
    p_ikat = 0.5;
    if (dV != 0)
        R = t * n_wj * (p_dkat/p_ikat) * (Z/dV) / p_site;
    else
        R = t * n_wj * (p_dkat/p_ikat) * (b-a) / p_site;

    // R = W * (p_dkat/p_ikat) * (tau_max-tau_t)/p_site;

    // cout << a << " " << b << " " << c << " " << dV << " " << Z << " " << R << " " << tau_kink << endl;

    // Metropolis Sampling
    if (rng.rand() < R){ // Accept
        
        // Add to acceptance counter
        ikat_accepts += 1;
                
        // Change kink that stored tail information to a regular kink
        paths[tail_idx].tau = tau_kink;
        paths[tail_idx].n = n_wi;
        paths[tail_idx].src = i;
        paths[tail_idx].dest = j;
        paths[tail_idx].prev = prev_i;
        paths[tail_idx].next = next_i;
        
        // Create the kinks on the destination site
        paths[num_kinks]=Kink(tau_t,n_wj,j,j,prev_j,num_kinks+1,
                              src_replica,dest_replica);
        paths[num_kinks+1]=Kink(tau_kink,n_j,j,i,num_kinks,next_j,
                              src_replica,dest_replica);
        
        // Set new worm tail index
        tail_idx = num_kinks;
                
        // "Connect" next of lower bound kink to worm tail
        paths[prev_j].next = num_kinks;
        
        // "Connect" prev of next kink to new kink
        if(next_j!=-1){paths[next_j].prev = num_kinks+1;}
        
        // If new kink is last kink on site j, update last kinks tracker vector
        if (next_j==-1){last_kinks[j]=num_kinks+1;}
        
        // Update number of kinks tracker
        num_kinks += 2;
        
        return;
            
        }
        else // Reject
            return;
        }

/*--------------------------------------------------------------------*/

void delete_kink_after_tail(vector<Kink> &paths, int &num_kinks,
                int &head_idx,int &tail_idx,
                int M, int N, double U, double mu, double t,
                vector<vector<int> > &adjacency_matrix, int total_nn,
                double beta, double eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                unsigned long long int &dkat_attempts,
                unsigned long long int &dkat_accepts,
                RNG &rng){
    
    // Variable declarations
    int prev,i,j,n_i,n_wi,n_j,n_wj,prev_i,prev_j,next_i,next_j,
    kink_idx_i,kink_idx_j,src_replica,dest_replica;
    double tau,tau_t,p_site,W,R,p_dkat,p_ikat,tau_prev_i,
    tau_kink,tau_max,dV_i,dV_j,tau_next_i,tau_next_j;
    
    // Update only possible if worm tail present
    if (tail_idx==-1){return;}
    
    // There has to be a regular kink after the worm tail
    if (paths[tail_idx].next==head_idx ||
        paths[tail_idx].next==-1){return;}

    // Need at least two sites to perform a spaceshift
    if (M<2){return;}

    // Indices of: upper bound kink, kink before tail, lower bound kink ; site j
    kink_idx_j = paths[tail_idx].next;
    next_j = paths[kink_idx_j].next;
    prev_j = paths[tail_idx].prev;

    // Times of: worm tail, kink before tail, lower bound kink; site j
    if (next_j!=-1)
        tau_next_j = paths[next_j].tau;
    else
        tau_next_j = beta;
    tau_kink = paths[kink_idx_j].tau;
    tau_t = paths[tail_idx].tau;
    
    // Only kinks in which the particle hops from j TO i can be deleted
    if ((paths[kink_idx_j].n-paths[tail_idx].n)>0){return;}
    
    // Retrieve worm tail site (j) and connecting site (i)
    j = paths[kink_idx_j].src;
    i = paths[kink_idx_j].dest;
    src_replica = paths[kink_idx_j].src_replica;
    dest_replica = paths[kink_idx_j].dest_replica;

    // Determine index of lower/upper bounds of flat where kink connects to (i)
    tau = 0.0;            // tau_prev_i candidate
    prev = i;           // prev_i candidate
    prev_i = i;         // this avoids "variable maybe not initialized" warning
    while (tau<tau_kink){
        // Set the lower bound index
        prev_i = prev;

        // Update lower bound index and tau candidates for next iteration
        prev = paths[prev].next;
        if (prev==-1){break;}
        tau = paths[prev].tau;
    }
    kink_idx_i = prev;
    next_i=paths[kink_idx_i].next;

    // Retrieve time of lower,upper bounds on connecting site (i)
    tau_prev_i = paths[prev_i].tau;
    if (next_i==-1){tau_next_i=beta;}
    else {tau_next_i = paths[next_i].tau;};
    
    // Deletion cannot interfere w/ kinks on other site
    if (tau_t <= tau_prev_i){return;}

    // Add to proposal counter
    dkat_attempts += 1;

    // Determine highest time at which kink could've been inserted
    if (tau_next_i<tau_next_j){tau_max=tau_next_i;}
    else {tau_max=tau_next_j;}

    // Probability of inverse move (ikah) choosing site where worm end is
    p_site = 1.0/total_nn;

    // Extract no. of particles in the flats adjacent to the new kink
    n_i = paths[prev_i].n;
    n_wi = n_i+1;
    n_j = paths[prev_j].n;
    n_wj = n_j+1;                   // "w": segment with the extra particle

    // Calculate the diagonal energy difference on both sites
    dV_i = (U/2.0)*(n_wi*(n_wi-1)-n_i*(n_i-1)) - mu*(n_wi-n_i);
    dV_j = (U/2.0)*(n_wj*(n_wj-1)-n_j*(n_j-1)) - mu*(n_wj-n_j);

    // Calculate the weight ratio W'/W
    W = t * n_wj * exp((dV_i-dV_j)*(tau_kink-tau_t));

    // Build the Metropolis ratio (R)
    p_dkat = 0.5;
    p_ikat = 0.5;
    R = W * (p_dkat/p_ikat) * (tau_max-tau_t)/p_site;
    R = 1.0/R;

    // Metropolis Sampling
    //boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    if (rng.rand() < R){ // Accept

        // Add to acceptance counter
        dkat_accepts += 1;
        
        // Stage 1: Delete kink on i
        if (paths[num_kinks-1].next!=-1)
            paths[paths[num_kinks-1].next].prev = kink_idx_i;
        paths[paths[num_kinks-1].prev].next = kink_idx_i;
        
        swap(paths[kink_idx_i],paths[num_kinks-1]);
        
        if (prev_i==num_kinks-1){prev_i=kink_idx_i;}
        else if (next_i==num_kinks-1){next_i=kink_idx_i;}
        else if (prev_j==num_kinks-1){prev_j=kink_idx_i;}
        else if (next_j==num_kinks-1){next_j=kink_idx_i;}
        else if (kink_idx_j==num_kinks-1){kink_idx_j=kink_idx_i;}
        else if (tail_idx==num_kinks-1){tail_idx=kink_idx_i;}
        else {;}
        
        if (head_idx==num_kinks-1){head_idx=kink_idx_i;}
        
        if (paths[kink_idx_i].next==-1){
            last_kinks[paths[kink_idx_i].src]=kink_idx_i;
        }
        
        if (next_i!=-1)
            paths[next_i].prev = prev_i;
        paths[prev_i].next = next_i;
        
        if (next_i==-1){last_kinks[i]=prev_i;}

        // Stage 2: Delete kink on j
        if (paths[num_kinks-2].next!=-1)
            paths[paths[num_kinks-2].next].prev = kink_idx_j;
        paths[paths[num_kinks-2].prev].next = kink_idx_j;
        
        swap(paths[kink_idx_j],paths[num_kinks-2]);
        
        if (prev_i==num_kinks-2){prev_i=kink_idx_j;}
        else if (next_i==num_kinks-2){next_i=kink_idx_j;}
        else if (prev_j==num_kinks-2){prev_j=kink_idx_j;}
        else if (next_j==num_kinks-2){next_j=kink_idx_j;}
        else if (tail_idx==num_kinks-2){tail_idx=kink_idx_j;}
        else {;}
        
        if (head_idx==num_kinks-2){head_idx=kink_idx_j;}
        
        if (paths[kink_idx_j].next==-1){
            last_kinks[paths[kink_idx_j].src]=kink_idx_j;
        }
        
        if (next_j!=-1)
            paths[next_j].prev = tail_idx;
        paths[tail_idx].next = next_j;
        
        if (next_j==-1){last_kinks[j]=tail_idx;}
        
        // Stage 3: Delete worm head on j
        if (paths[num_kinks-3].next!=-1)
            paths[paths[num_kinks-3].next].prev = tail_idx;
        paths[paths[num_kinks-3].prev].next = tail_idx;
        
        swap(paths[tail_idx],paths[num_kinks-3]);
        
        if (prev_i==num_kinks-3){prev_i=tail_idx;}
        else if (next_i==num_kinks-3){next_i=tail_idx;}
        else if (prev_j==num_kinks-3){prev_j=tail_idx;}
        else if (next_j==num_kinks-3){next_j=tail_idx;}
        else {;}
        
        if (head_idx==num_kinks-3){head_idx=tail_idx;}
        
        if (paths[tail_idx].next==-1){
            last_kinks[paths[tail_idx].src]=tail_idx;
        }
        
        if (next_j!=-1)
            paths[next_j].prev = prev_j;
        paths[prev_j].next = next_j;
        
        if (next_j==-1){last_kinks[j]=prev_j;}
        
        // Stage 4: Insert worm tail on i
        paths[num_kinks-3]=Kink(tau_t,n_wi,i,i,prev_i,next_i,
                                src_replica,dest_replica);
        
        tail_idx = num_kinks-3;
        
        paths[prev_i].next = tail_idx;
        if(next_i!=-1){paths[next_i].prev = tail_idx;}
        
        if (next_i==-1){last_kinks[i]=tail_idx;}
        
        // Update number of kinks tracker
        num_kinks -= 2;
        
        return;

    }
    else // Reject
        return;
    }

/*--------------------------------------------------------------------*/

void delete_kink_after_tail_2(vector<Kink> &paths, int &num_kinks,
                int &head_idx,int &tail_idx,
                int M, int N, double U, double mu, double t,
                vector<vector<int> > &adjacency_matrix, int total_nn,
                double beta, double eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                unsigned long long int &dkat_attempts,
                unsigned long long int &dkat_accepts,
                RNG &rng){
    
    // Variable declarations
    int prev,i,j,n_i,n_wi,n_j,n_wj,prev_i,prev_j,next_i,next_j,
    kink_idx_i,kink_idx_j,src_replica,dest_replica;
    double tau,tau_t,p_site,R,p_dkat,p_ikat,tau_prev_i,
    tau_kink,tau_max,dV_i,dV_j,tau_next_i,tau_next_j,dV;
    
    // Update only possible if worm tail present
    if (tail_idx==-1){return;}
    
    // There has to be a regular kink after the worm tail
    if (paths[tail_idx].next==head_idx ||
        paths[tail_idx].next==-1){return;}

    // Need at least two sites to perform a spaceshift
    if (M<2){return;}

    // Indices of: upper bound kink, kink before tail, lower bound kink ; site j
    kink_idx_j = paths[tail_idx].next;
    next_j = paths[kink_idx_j].next;
    prev_j = paths[tail_idx].prev;

    // Times of: worm tail, kink before tail, lower bound kink; site j
    if (next_j!=-1)
        tau_next_j = paths[next_j].tau;
    else
        tau_next_j = beta;
    tau_kink = paths[kink_idx_j].tau;
    tau_t = paths[tail_idx].tau;
    
    // Only kinks in which the particle hops from j TO i can be deleted
    if ((paths[kink_idx_j].n-paths[tail_idx].n)>0){return;}
    
    // Retrieve worm tail site (j) and connecting site (i)
    j = paths[kink_idx_j].src;
    i = paths[kink_idx_j].dest;
    src_replica = paths[kink_idx_j].src_replica;
    dest_replica = paths[kink_idx_j].dest_replica;

    // Determine index of lower/upper bounds of flat where kink connects to (i)
    tau = 0.0;            // tau_prev_i candidate
    prev = i;           // prev_i candidate
    prev_i = i;         // this avoids "variable maybe not initialized" warning
    while (tau<tau_kink){
        // Set the lower bound index
        prev_i = prev;

        // Update lower bound index and tau candidates for next iteration
        prev = paths[prev].next;
        if (prev==-1){break;}
        tau = paths[prev].tau;
    }
    kink_idx_i = prev;
    next_i=paths[kink_idx_i].next;

    // Retrieve time of lower,upper bounds on connecting site (i)
    tau_prev_i = paths[prev_i].tau;
    if (next_i==-1){tau_next_i=beta;}
    else {tau_next_i = paths[next_i].tau;};
    
    // Deletion cannot interfere w/ kinks on other site
    if (tau_t <= tau_prev_i){return;}

    // Add to proposal counter
    dkat_attempts += 1;

    // Determine highest time at which kink could've been inserted
    if (tau_next_i<tau_next_j){tau_max=tau_next_i;}
    else {tau_max=tau_next_j;}

    // Probability of inverse move (ikah) choosing site where worm end is
    p_site = 1.0/total_nn;

    // Extract no. of particles in the flats adjacent to the new kink
    n_i = paths[prev_i].n;
    n_wi = n_i+1;
    n_j = paths[prev_j].n;
    n_wj = n_j+1;                   // "w": segment with the extra particle

    // Calculate the diagonal energy difference on both sites
    dV_i = (U/2.0)*(n_wi*(n_wi-1)-n_i*(n_i-1)) - mu*(n_wi-n_i);
    dV_j = (U/2.0)*(n_wj*(n_wj-1)-n_j*(n_j-1)) - mu*(n_wj-n_j);
    dV = dV_j - dV_i;
    // if (dV == 0){cout << "NOOOO" << endl;}

    // Calculate the weight ratio W'/W
    // W = t * n_wj * exp((dV_i-dV_j)*(tau_kink-tau_t));

    // inverse move (insert kink after tail) truncated sampling
    double Z,a,b,c;
    a = tau_t;
    b = tau_max;
    c = dV; 
    Z = 1.0 - exp(-c*(b-a));

    // Build the Metropolis ratio (R)
    p_dkat = 0.5;
    p_ikat = 0.5;
    if (dV != 0)
        R = t * n_wj * (p_dkat/p_ikat) * (Z/dV) / p_site;
    else
        R = t * n_wj * (p_dkat/p_ikat) * (b-a) / p_site;
    // R = W * (p_dkat/p_ikat) * (tau_max-tau_t)/p_site;
    R = 1.0/R;

    // Metropolis Sampling
    //boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    if (rng.rand() < R){ // Accept

        // Add to acceptance counter
        dkat_accepts += 1;
        
        // Stage 1: Delete kink on i
        if (paths[num_kinks-1].next!=-1)
            paths[paths[num_kinks-1].next].prev = kink_idx_i;
        paths[paths[num_kinks-1].prev].next = kink_idx_i;
        
        swap(paths[kink_idx_i],paths[num_kinks-1]);
        
        if (prev_i==num_kinks-1){prev_i=kink_idx_i;}
        else if (next_i==num_kinks-1){next_i=kink_idx_i;}
        else if (prev_j==num_kinks-1){prev_j=kink_idx_i;}
        else if (next_j==num_kinks-1){next_j=kink_idx_i;}
        else if (kink_idx_j==num_kinks-1){kink_idx_j=kink_idx_i;}
        else if (tail_idx==num_kinks-1){tail_idx=kink_idx_i;}
        else {;}
        
        if (head_idx==num_kinks-1){head_idx=kink_idx_i;}
        
        if (paths[kink_idx_i].next==-1){
            last_kinks[paths[kink_idx_i].src]=kink_idx_i;
        }
        
        if (next_i!=-1)
            paths[next_i].prev = prev_i;
        paths[prev_i].next = next_i;
        
        if (next_i==-1){last_kinks[i]=prev_i;}

        // Stage 2: Delete kink on j
        if (paths[num_kinks-2].next!=-1)
            paths[paths[num_kinks-2].next].prev = kink_idx_j;
        paths[paths[num_kinks-2].prev].next = kink_idx_j;
        
        swap(paths[kink_idx_j],paths[num_kinks-2]);
        
        if (prev_i==num_kinks-2){prev_i=kink_idx_j;}
        else if (next_i==num_kinks-2){next_i=kink_idx_j;}
        else if (prev_j==num_kinks-2){prev_j=kink_idx_j;}
        else if (next_j==num_kinks-2){next_j=kink_idx_j;}
        else if (tail_idx==num_kinks-2){tail_idx=kink_idx_j;}
        else {;}
        
        if (head_idx==num_kinks-2){head_idx=kink_idx_j;}
        
        if (paths[kink_idx_j].next==-1){
            last_kinks[paths[kink_idx_j].src]=kink_idx_j;
        }
        
        if (next_j!=-1)
            paths[next_j].prev = tail_idx;
        paths[tail_idx].next = next_j;
        
        if (next_j==-1){last_kinks[j]=tail_idx;}
        
        // Stage 3: Delete worm head on j
        if (paths[num_kinks-3].next!=-1)
            paths[paths[num_kinks-3].next].prev = tail_idx;
        paths[paths[num_kinks-3].prev].next = tail_idx;
        
        swap(paths[tail_idx],paths[num_kinks-3]);
        
        if (prev_i==num_kinks-3){prev_i=tail_idx;}
        else if (next_i==num_kinks-3){next_i=tail_idx;}
        else if (prev_j==num_kinks-3){prev_j=tail_idx;}
        else if (next_j==num_kinks-3){next_j=tail_idx;}
        else {;}
        
        if (head_idx==num_kinks-3){head_idx=tail_idx;}
        
        if (paths[tail_idx].next==-1){
            last_kinks[paths[tail_idx].src]=tail_idx;
        }
        
        if (next_j!=-1)
            paths[next_j].prev = prev_j;
        paths[prev_j].next = next_j;
        
        if (next_j==-1){last_kinks[j]=prev_j;}
        
        // Stage 4: Insert worm tail on i
        paths[num_kinks-3]=Kink(tau_t,n_wi,i,i,prev_i,next_i,
                                src_replica,dest_replica);
        
        tail_idx = num_kinks-3;
        
        paths[prev_i].next = tail_idx;
        if(next_i!=-1){paths[next_i].prev = tail_idx;}
        
        if (next_i==-1){last_kinks[i]=tail_idx;}
        
        // Update number of kinks tracker
        num_kinks -= 2;
        
        return;

    }
    else // Reject
        return;
    }

/*--------------------------------------------------------------------*/

void insert_swap_kink(vector<vector<Kink> > &paths,
                      vector<int> &num_kinks,
                int num_replicas, int replica_idx,
                vector<int> &sub_sites, vector <int> &swapped_sites,
                vector<int> &swap_kinks, int &num_swaps,
                int l_A, int m_A,
                vector<int> &head_idx,vector<int> &tail_idx,
                int M, int N, double U, double mu, double t,
                vector<vector<int> > &adjacency_matrix, int total_nn,
                double beta,double eta,bool canonical,vector<double> &N_tracker,
                vector<int> &N_zero,vector<int> &N_beta,
                vector<vector<int> > &last_kinks,
                unsigned long long int &insert_swap_kink_attempts,
                unsigned long long int &insert_swap_kink_accepts,
                RNG &rng){
    
    // Note: ONE AND TWO DIMENSIONAL FOR NOW
    
    // Variable declarations
    int src_replica,dest_replica,n_src,n_dest,next,next_swap_site,prev_src,
    prev_dest,next_src,next_dest,num_kinks_src,
    num_kinks_dest;
    double tau;
    
    // Need at least two replicas to perform a spaceshift
    if (paths.size()<2){return;}

    // Can't perform update if SWAP region is full
    if (num_swaps==m_A){return;}
    
    insert_swap_kink_attempts+=1;
        
    // Retrieve source replica index and randomly choose destination replica
    src_replica = replica_idx;
    if (num_replicas==2){
        if (src_replica==0){dest_replica=1;}
        else {dest_replica=0;}
    }
    else{ // more than two replicas
        cout<<"ERROR: Only one or two replicas are valid at the moment."<<endl;
        exit(1);
    }

    // Propose the next site to swap
    next_swap_site = sub_sites[num_swaps];
     
    /*---------TEST THIS!!!*----------*/
    // Check if no. of particles at beta/2 is the same on both replicas
    // Source Replica
    tau = 0.0;
    next = next_swap_site; // next variable refers to "next kink" in worldline
    n_src = -1;
    prev_src = -1;
    while (tau<0.5*beta){
        n_src = paths[src_replica][next].n;
        prev_src = next;
        next = paths[src_replica][next].next;

        if (next==-1){break;}
        tau = paths[src_replica][next].tau;
    }
    next_src = next;
    // Destination Replica
    tau = 0.0;
    next = next_swap_site;
    n_dest=-1;
    prev_dest = -1;
    while (tau<0.5*beta){
        n_dest = paths[dest_replica][next].n;
        prev_dest = next;
        next = paths[dest_replica][next].next;

        if (next==-1){break;}
        tau = paths[dest_replica][next].tau;
    }
    next_dest = next;
    
    if (n_src!=n_dest){return;}
                    
    // Build and insert kinks to the paths of the src and the dest replica
    num_kinks_src = num_kinks[src_replica];
    num_kinks_dest = num_kinks[dest_replica];
    paths[src_replica][num_kinks_src] =
                              Kink(beta/2.0,n_src,next_swap_site,next_swap_site,
                              prev_src,next_src,
                              src_replica,dest_replica);
    paths[dest_replica][num_kinks_dest] =
                              Kink(beta/2.0,n_src,next_swap_site,next_swap_site,
                              prev_dest,next_dest,
                              dest_replica,src_replica);
    // Had to change the meaning of src_replica and dest_replica ATTRIBUTES
    // The now actually have directional meaning. From origin replica to other.

    // Connect next of prev_src to swap_kink
    paths[src_replica][prev_src].next = num_kinks_src;
    
    // Connect prev of next_src to swap_kink
    if (next_src!=-1)
        paths[src_replica][next_src].prev = num_kinks_src;
    
    // Connect next of prev_dest to swap_kink
    paths[dest_replica][prev_dest].next = num_kinks_dest;
    
    // Connect prev of next_dest to swap_kink
    if (next_dest!=-1)
        paths[dest_replica][next_dest].prev = num_kinks_dest;
    
    // Edit the last kinks vector of each replica if necessary
    if (next_src==-1){
        last_kinks[src_replica][next_swap_site] = num_kinks_src;
    }
    if (next_dest==-1){
        last_kinks[dest_replica][next_swap_site] = num_kinks_dest;
    }
    
    // Update number of swapped sites tracker
    num_swaps+=1;

    // Update number of kinks tracker of each replica
    num_kinks[src_replica] += 1;
    num_kinks[dest_replica] += 1;
     
    insert_swap_kink_accepts+=1;

    return;
}

/*--------------------------------------------------------------------*/

void delete_swap_kink(vector<vector<Kink> > &paths, vector<int> &num_kinks,
                int num_replicas, int replica_idx,
                vector<int> &sub_sites, vector <int> &swapped_sites,
                vector<int> &swap_kinks, int &num_swaps,
                int l_A, int m_A,
                vector<int> &head_idx,vector<int> &tail_idx,
                int M, int N, double U, double mu, double t,
                vector<vector<int> > &adjacency_matrix, int total_nn,
                double beta,double eta,bool canonical,vector<double> &N_tracker,
                vector<int> &N_zero,vector<int> &N_beta,
                vector<vector<int> > &last_kinks,
                unsigned long long int &delete_swap_kink_attempts,
                unsigned long long int &delete_swap_kink_accepts,
                RNG &rng){
    
    // Note: ONE AND TWO DIMENSIONAL FOR NOW
    
    // Variable declarations
    int src_replica,dest_replica,next,prev_src,prev_dest,next_src,next_dest,
    num_kinks_src,num_kinks_dest,site_to_unswap,kink_out_of_src,
    kink_out_of_dest,n_src_left,n_src_right,n_dest_left,n_dest_right;
    
    // Need at least two replicas to perform delete_swap_kink
    if (paths.size()<2){return;}
    
    // Need at least one swapped site to perform delete_swap_kink
    if (num_swaps==0){return;}
    
    delete_swap_kink_attempts+=1;

    // Retrieve source replica index and randomly choose destination replica
    src_replica = replica_idx;
    if (num_replicas==2){
        if (src_replica==0){dest_replica=1;}
        else {dest_replica=0;}
    }
    else{ // more than two replicas
        cout<<"ERROR: Only one or two replicas are valid at the moment."<<endl;
        exit(1);
    }
    
    // Get the number of kinks on each replica
    num_kinks_src = num_kinks[src_replica];
    num_kinks_dest = num_kinks[dest_replica];
    
    // Randomly choose a swapped site to unswap
    site_to_unswap = sub_sites[num_swaps-1];
     
    // Get swap kink indices
    // source replica
    next = site_to_unswap;
    while (paths[src_replica][next].dest_replica==
           paths[src_replica][next].src_replica){
        next = paths[src_replica][next].next;
    }
    kink_out_of_src = next;

    // destination replica
    next = site_to_unswap;
    while (paths[dest_replica][next].dest_replica==
           paths[dest_replica][next].src_replica){
        next = paths[dest_replica][next].next;
    }
    kink_out_of_dest = next;
    
    // Get lower and upper bounds of the flat interval on each replica
    prev_src = paths[src_replica][kink_out_of_src].prev;
    next_src = paths[src_replica][kink_out_of_src].next;
    prev_dest = paths[dest_replica][kink_out_of_dest].prev;
    next_dest = paths[dest_replica][kink_out_of_dest].next;
    
    // Get number of particles to the left and righ of the swap kinks
    n_src_left = paths[src_replica][prev_src].n;
    n_src_right = paths[src_replica][kink_out_of_src].n;
    n_dest_left = paths[dest_replica][prev_dest].n;
    n_dest_right = paths[dest_replica][kink_out_of_dest].n;
    
    // Can only delete swap kink if pre/post no. of particles is the same
    if (n_src_left!=n_src_right){return;}
    if (n_dest_left!=n_dest_right){return;}
     
    // Stage 1: delete kink coming out of source replica
    // Modify links to kink at end of paths vector that will be swapped
    if (paths[src_replica][num_kinks_src-1].next!=-1)
        paths[src_replica][paths[src_replica][num_kinks_src-1].next].prev=
                                                                kink_out_of_src;
    paths[src_replica][paths[src_replica][num_kinks_src-1].prev].next=
                                                                kink_out_of_src;
    
    swap(paths[src_replica][kink_out_of_src],
         paths[src_replica][num_kinks_src-1]);
    
    // Important kinks might've been swapped. Correct if so.
    if (prev_src==num_kinks_src-1)
        prev_src=kink_out_of_src;
    else if (next_src==num_kinks_src-1)
        next_src=kink_out_of_src;
    else{;}
    
    // Head & tail separately b.c they could've been one of the kinks above too
    if (head_idx[src_replica]==num_kinks_src-1)
        head_idx[src_replica]=kink_out_of_src;
    else if (tail_idx[src_replica]==num_kinks_src-1)
        tail_idx[src_replica]=kink_out_of_src;
    else {;}
    
    // The kink sent to where deleted kink was might be last on it's site
    if (paths[src_replica][kink_out_of_src].next==-1){
        last_kinks[src_replica][paths[src_replica][kink_out_of_src].src]=
                                                               kink_out_of_src;
    }
    
    // Reconnect upper and lower bounds of the flat
    if (next_src!=-1)
        paths[src_replica][next_src].prev = prev_src;
    paths[src_replica][prev_src].next = next_src;
    
    // Lower bound of flat could be the last kink in the site
    if (next_src==-1){last_kinks[src_replica][site_to_unswap]=prev_src;}
    
    // Stage 3: delete kink coming out of destination replica
    // Modify links to kink at end of paths vector that will be swapped
    if (paths[dest_replica][num_kinks_dest-1].next!=-1)
        paths[dest_replica][paths[dest_replica][num_kinks_dest-1].next].prev=
                                                               kink_out_of_dest;
    paths[dest_replica][paths[dest_replica][num_kinks_dest-1].prev].next=
                                                               kink_out_of_dest;
    
    swap(paths[dest_replica][kink_out_of_dest],
         paths[dest_replica][num_kinks_dest-1]);
    
    // Important kinks might've been swapped. Correct if so.
    if (prev_dest==num_kinks_dest-1)
        prev_dest=kink_out_of_dest;
    else if (next_dest==num_kinks_dest-1)
        next_dest=kink_out_of_dest;
    else{;}
    
    // Head & tail separately b.c they could've been one of the kinks above too
    if (head_idx[dest_replica]==num_kinks_dest-1)
        head_idx[dest_replica]=kink_out_of_dest;
    else if (tail_idx[dest_replica]==num_kinks_dest-1)
        tail_idx[dest_replica]=kink_out_of_dest;
    else{;}
    
    // The kink sent to where deleted kink was might be last on it's site
    if (paths[dest_replica][kink_out_of_dest].next==-1){
        last_kinks[dest_replica][paths[dest_replica][kink_out_of_dest].src]=
                                                               kink_out_of_dest;
    }
    
    // Reconnect upper and lower bounds of the flat
    if (next_dest!=-1)
        paths[dest_replica][next_dest].prev = prev_dest;
    paths[dest_replica][prev_dest].next = next_dest;
    
    // Lower bound of flat could be the last kink in the site
    if (next_dest==-1){last_kinks[dest_replica][site_to_unswap]=prev_dest;}
    
    // Modify number of swaps tracker
    num_swaps-=1;
    
    // Modify number of kinks trackers for each replica
    num_kinks[src_replica]-=1;
    num_kinks[dest_replica]-=1;
    
    delete_swap_kink_accepts+=1;
    
    return;
}

/*--------------------------------------------------------------------*/

void swap_timeshift_head(vector<vector<Kink> > &paths, vector<int> &num_kinks,
               int num_replicas, int replica_idx,
               vector<int> &sub_sites, vector <int> &swapped_sites,
               vector<int> &swap_kinks, int &num_swaps,
               int l_A, int m_A,
               vector<int> &head_idx,vector<int> &tail_idx,
               int M, int N, double U, double mu, double t,
               vector<vector<int> > &adjacency_matrix, int total_nn,
               double beta,double eta,bool canonical,vector<double> &N_tracker,
               vector<int> &N_zero,vector<int> &N_beta,
               vector<vector<int> > &last_kinks,
               unsigned long long int &swap_advance_head_attempts,
               unsigned long long int &swap_advance_head_accepts,
               unsigned long long int &swap_recede_head_attempts,
               unsigned long long int &swap_recede_head_accepts,
               RNG &rng){
    
    // Variable declarations
    int n,worm_end_idx,src_replica,dest_replica,head_idx_0,
    head_idx_1,prev_src,next_src,prev_dest,next_dest,worm_end_site,
    kink_out_of_dest,n_after_worm_end,n_after_swap_kink,
    num_kinks_src,num_kinks_dest,current_kink,n_before_swap_kink;
    double tau,tau_prev,tau_next,dV,R,tau_new,Z,l_path_src,
    l_path_dest,dN_src,dN_dest;
    bool swap_in_front,is_advance,is_over_swap;
    vector<Kink> paths_src,paths_dest;
    
    // Need at least two replicas to perform a spaceshift
    if (paths.size()<2){return;}
    
    // Need at least one swap kink to perform a timeshift through swap kink
    if (num_swaps==0){return;}
    
    // Retrieve worm head indices. -1 means no worm head
    head_idx_0 = head_idx[0];
    head_idx_1 = head_idx[1];
    
    // There had to be STRICTLY ONE worm head to timeshift over swap kink
    if (head_idx_0!=-1 && head_idx_1!=-1){return;} // two worms present
    if (head_idx_0==-1 && head_idx_1==-1){return;} // no worms present
    
    // Choose the "source" and "destination" replica
    if (head_idx_0!=-1){src_replica=0;}
    else{src_replica=1;}
    dest_replica = 1-src_replica;
    
    // Get index of the worm head to be moved
    worm_end_idx=head_idx[src_replica];
        
    // Extract worm head time and site
    tau = paths[src_replica][worm_end_idx].tau;
    worm_end_site = paths[src_replica][worm_end_idx].src;
    n = paths[src_replica][worm_end_idx].n;
    
    // Get lower and upper adjacent kinks of worm head to be moved
    // NOTE: One of the two bounds might be the swap kink.
    prev_src=paths[src_replica][worm_end_idx].prev;
    next_src=paths[src_replica][worm_end_idx].next;
    
    // Check if worm head is adjacent to a swap kink
    if (next_src!=-1){
        if (paths[src_replica][next_src].src_replica!=
            paths[src_replica][next_src].dest_replica){swap_in_front=true;}
        else if (paths[src_replica][prev_src].src_replica!=
                 paths[src_replica][prev_src].dest_replica){swap_in_front=false;}
        else {return;}
    }
    else{
        if (paths[src_replica][prev_src].src_replica!=
            paths[src_replica][prev_src].dest_replica){swap_in_front=false;}
        else {return;}
    }
                
    // Get index of central time slice at destination replica
    current_kink = worm_end_site; // next refers to index of next kink on site
    while (paths[dest_replica][current_kink].dest_replica==
           paths[dest_replica][current_kink].src_replica){
        current_kink = paths[dest_replica][current_kink].next;
    }
    kink_out_of_dest = current_kink;
    
    // Get indices of kinks before & after swap kink in destination replica
    prev_dest = paths[dest_replica][kink_out_of_dest].prev;
    next_dest = paths[dest_replica][kink_out_of_dest].next;
    
    // Determine the lower and upper bounds of the worm end to be timeshifted
    if (swap_in_front){
        if (next_dest!=-1)
            tau_next = paths[dest_replica][next_dest].tau;
        else
            tau_next = beta;
        tau_prev = paths[src_replica][prev_src].tau;
    }
    else{ // swap kink behind
        if (next_src!=-1)
            tau_next = paths[src_replica][next_src].tau;
        else
            tau_next=beta;
        tau_prev = paths[dest_replica][prev_dest].tau;
    }
    
    // Calculate change in diagonal energy
//    shift_head=true; // we are always moving head in this update. set to true.
//    dV=U*(n-!shift_head)-mu;
    dV=U*n-mu;
    
    // To make acceptance ratio unity,shift tail needs to sample w/ dV=eps-eps_w
//    if (!shift_head){dV *= -1;} // dV=eps-eps_w
    
    //boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    // Sample the new time of the worm end from truncated exponential dist.
    /*:::::::::::::::::::: Truncated Exponential RVS :::::::::::::::::::::::::*/
    Z = 1.0 - exp(-dV*(tau_next-tau_prev));
    tau_new = tau_prev - log(1.0-Z*rng.rand())  / dV;
    /*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
    
    // Check if advance or recede
    if (tau_new>tau){
        is_advance=true;
        swap_advance_head_attempts+=1;
    }
    else {
        is_advance=false;
        swap_recede_head_attempts+=1;
    }
    
    // Check if timeshift takes the worm end over the swap kink
    is_over_swap=false;
    if ( (is_advance && swap_in_front && tau_new>beta/2) ||
         (!is_advance && !swap_in_front && tau_new<beta/2) ){is_over_swap=true;}
        
    // Determine the length of path to be modified in each replica
    if (is_over_swap){
        l_path_src = beta/2 - tau;
        l_path_dest = tau_new - beta/2;
    }
    else{ // did not go over swap
        l_path_src = tau_new - tau;
        l_path_dest = 0;
    }
    
    // Determine the total particle change based on wormend to be shifted
    dN_src  = 1.0 * l_path_src/beta;
    dN_dest = 1.0 * l_path_dest/beta;
    
    // Canonical simulations: Restrict updates to interval N:(N-1,N+1)
    if (canonical){
        if ((N_tracker[src_replica]+dN_src)   < (N-1)  ||
            (N_tracker[src_replica]+dN_src)   > (N+1)  ||
            (N_tracker[dest_replica]+dN_dest) < (N-1)  ||
            (N_tracker[dest_replica]+dN_dest) > (N+1)){return;}
    }
    
    // Get number of particles after: worm end @ src & central kink @ dest
    n_after_worm_end = paths[src_replica][worm_end_idx].n;
    n_after_swap_kink = paths[dest_replica][kink_out_of_dest].n;
    n_before_swap_kink = paths[dest_replica][prev_dest].n;
        
    // Get number of kinks in source and destination replicas (before update)
    num_kinks_src = num_kinks[src_replica];
    num_kinks_dest = num_kinks[dest_replica];
        
    // Build the Metropolis condition (R)
    R = 1.0; // Sampling worm end time from truncated exponential makes R unity.

    // Metropolis sampling
    if (rng.rand() < R){
        if (!is_over_swap){ // worm end does not go over swap kink
            paths[src_replica][worm_end_idx].tau = tau_new;
            N_tracker[src_replica] += dN_src;
            if (is_advance){swap_advance_head_accepts+=1;}
            else {swap_recede_head_accepts+=1;}
        }
        else{ // We go Over Swap

            if (is_advance){ // advance OVER SWAP

                /*--------- Deletion of worm end from SOURCE replica ---------*/
                
                // num_kinks_src-1 will be swapped. Modify links to it
                if (paths[src_replica][num_kinks_src-1].next!=-1){
paths[src_replica][paths[src_replica][num_kinks_src-1].next].prev=worm_end_idx;
                }
paths[src_replica][paths[src_replica][num_kinks_src-1].prev].next=worm_end_idx;
                
                swap(paths[src_replica][worm_end_idx],
                     paths[src_replica][num_kinks_src-1]);
                
                // Upper or lower bound of flat could've been swapped. Correct.
                if (next_src==num_kinks_src-1){next_src=worm_end_idx;}
                else if (prev_src==num_kinks_src-1){prev_src=worm_end_idx;}
                else {;}
                
                // Tail could've been swapped. Correct if so.
                if (tail_idx[src_replica]==num_kinks_src-1)
                    tail_idx[src_replica]=worm_end_idx;
                
                // Whatever kink was swapped could've been the last on its site
                if (paths[src_replica][worm_end_idx].next==-1){
                  last_kinks[src_replica][paths[src_replica][worm_end_idx].src]
                    =worm_end_idx;
                }

                // Reconnect upper and lower bounds of the flat
                if (next_src!=-1)
                    paths[src_replica][next_src].prev = prev_src;
                paths[src_replica][prev_src].next = next_src;
                
                // Deactivate the worm end
                head_idx[src_replica]=-1;
                                
                // Update trackers for: no. of active kinks,total particles
                num_kinks[src_replica]-=1;
                N_tracker[src_replica] += dN_src;

                /*------- Insertion of worm end in DESTINATION replica -------*/
                
                // Activate first available kink
                paths[dest_replica][num_kinks_dest]=
                Kink(tau_new,n_after_swap_kink,worm_end_site,worm_end_site,
                     kink_out_of_dest,next_dest,dest_replica,dest_replica);
                
                // Save head index
                head_idx[dest_replica]=num_kinks_dest;
                
                // Update number of particles after swap kink in dest replica
                paths[dest_replica][kink_out_of_dest].n=n_after_swap_kink+1;
                
                // Add to acceptance counter
                swap_advance_head_accepts+=1;
                
                // Modify links of swap kink and next_dest kink
                if (next_dest!=-1)
                    paths[dest_replica][next_dest].prev=num_kinks_dest;
                paths[dest_replica][kink_out_of_dest].next=num_kinks_dest;
                
                // Update trackers for: no. of active kinks,total particles
                num_kinks[dest_replica] += 1;
                N_tracker[dest_replica] += dN_dest;

                // Created kink might be last on its site
                if (next_dest==-1)
                    last_kinks[dest_replica][worm_end_site]=
                    num_kinks_dest;
                                
            }
            else{ // Recede OVER SWAP

                // Cannot recede head over swap if no particles on destination
                if (n_before_swap_kink==0){return;}
                
                /*--------- Deletion of worm end from SOURCE replica ---------*/
                // num_kinks_src-1 will be swapped. Modify links to it
                if (paths[src_replica][num_kinks_src-1].next!=-1){
paths[src_replica][paths[src_replica][num_kinks_src-1].next].prev=worm_end_idx;
                }
paths[src_replica][paths[src_replica][num_kinks_src-1].prev].next=worm_end_idx;
                
                swap(paths[src_replica][worm_end_idx],
                     paths[src_replica][num_kinks_src-1]);
                
                // Upper or lower bound of flat could've been swapped. Correct.
                if (next_src==num_kinks_src-1){next_src=worm_end_idx;}
                else if (prev_src==num_kinks_src-1){prev_src=worm_end_idx;}
                else {;}
                
                // Tail could've been swapped. Correct if so.
                if (tail_idx[src_replica]==num_kinks_src-1)
                    tail_idx[src_replica]=worm_end_idx;
                
                // Whatever kink was swapped could've been the last on its site
                if (paths[src_replica][worm_end_idx].next==-1){
                  last_kinks[src_replica][paths[src_replica][worm_end_idx].src]
                    =worm_end_idx;
                }
                
                // Modify particles in prev_src
                paths[src_replica][prev_src].n=n_after_worm_end;

                // Reconnect upper and lower bounds of the flat
                if (next_src!=-1)
                    paths[src_replica][next_src].prev = prev_src;
                paths[src_replica][prev_src].next = next_src;
                
                // Deactivate the worm end
                head_idx[src_replica]=-1;
                                
                // Update trackers for: no. of active kinks,total particles
                num_kinks[src_replica]-=1;
                N_tracker[src_replica] += dN_src;
                
                // Swap kink on src might be last kink on it's site
                if (next_src==-1){
                    last_kinks[src_replica][worm_end_site]=prev_src;
                }
                

                /*------- Insertion of worm end in DESTINATION replica -------*/

                // Activate first available kink
                paths[dest_replica][num_kinks_dest]=
                Kink(tau_new,n_before_swap_kink-1,worm_end_site,worm_end_site,
                     prev_dest,kink_out_of_dest,dest_replica,dest_replica);
                
                // Save head index
                head_idx[dest_replica]=num_kinks_dest;
                
                // Add to acceptance counter
                swap_recede_head_accepts+=1;
                
                // Modify links to worm head in destination replica
                paths[dest_replica][kink_out_of_dest].prev=num_kinks_dest;
                paths[dest_replica][prev_dest].next=num_kinks_dest;
                
                // Update trackers for: no. of active kinks,total particles
                num_kinks[dest_replica] += 1;
                N_tracker[dest_replica] += dN_dest;
                
//                cout<<"Receded head over swap AFTER (left/right fock state)"<<endl;
//                for (int i=0; i<1; i++){
//                    cout<<paths[src_replica][paths[src_replica][prev_src].prev].n;
//                }
//                cout << " || ";
//                for (int i=0; i<1; i++){
//                    cout<<paths[src_replica][prev_src].n;
//                }
//
//                cout << "    ";
//
//                for (int i=0; i<1; i++){
//                    cout<<paths[dest_replica][num_kinks_dest].n;
//                }
//                cout << " || ";
//                for (int i=0; i<1; i++){
//                    cout<<paths[dest_replica][kink_out_of_dest].n;
//                }
//
//                cout << endl;
                
//                cout << "2 (recede)" << endl;
            }
        }
        return;
    }
    else // Reject
        return;
}

/*--------------------------------------------------------------------*/

void swap_timeshift_tail(vector<vector<Kink> > &paths, vector<int> &num_kinks,
               int num_replicas, int replica_idx,
               vector<int> &sub_sites, vector <int> &swapped_sites,
               vector<int> &swap_kinks, int &num_swaps,
               int l_A, int m_A,
               vector<int> &head_idx,vector<int> &tail_idx,
               int M, int N, double U, double mu, double t,
               vector<vector<int> > &adjacency_matrix, int total_nn,
               double beta,double eta,bool canonical,vector<double> &N_tracker,
               vector<int> &N_zero,vector<int> &N_beta,
               vector<vector<int> > &last_kinks,
               unsigned long long int &swap_advance_tail_attempts,
               unsigned long long int &swap_advance_tail_accepts,
               unsigned long long int &swap_recede_tail_attempts,
               unsigned long long int &swap_recede_tail_accepts,
               RNG &rng){
    
    // Variable declarations
    int n,worm_end_idx,src_replica,dest_replica,tail_idx_0,
    tail_idx_1,prev_src,next_src,prev_dest,next_dest,worm_end_site,
    kink_out_of_dest,n_after_worm_end,n_after_swap_kink,
    num_kinks_src,num_kinks_dest,current_kink,n_before_swap_kink;
    double tau,tau_prev,tau_next,dV,R,tau_new,Z,l_path_src,
    l_path_dest,dN_src,dN_dest;
    bool swap_in_front,is_advance,is_over_swap;
    vector<Kink> paths_src,paths_dest;
    
    // Need at least two replicas to perform a spaceshift
    if (paths.size()<2){return;}
    
    // Need at least one swap kink to perform a timeshift through swap kink
    if (num_swaps==0){return;}
    
    // Retrieve worm tail indices. -1 means no worm tail
    tail_idx_0 = tail_idx[0];
    tail_idx_1 = tail_idx[1];
    
    // There had to be STRICTLY ONE worm tail to timeshift over swap kink
    if (tail_idx_0!=-1 && tail_idx_1!=-1){return;} // two worms present
    if (tail_idx_0==-1 && tail_idx_1==-1){return;} // no worms present
    
    // Choose the "source" and "destination" replica
    if (tail_idx_0!=-1){src_replica=0;}
    else{src_replica=1;}
    dest_replica = 1-src_replica;
    
    // Get index of the worm tail to be moved
    worm_end_idx=tail_idx[src_replica];
        
    // Extract worm tail time and site
    tau = paths[src_replica][worm_end_idx].tau;
    worm_end_site = paths[src_replica][worm_end_idx].src;
    n = paths[src_replica][worm_end_idx].n;
    
    // Get lower and upper adjacent kinks of worm tail to be moved
    // NOTE: One of the two bounds might be the swap kink.
    prev_src=paths[src_replica][worm_end_idx].prev;
    next_src=paths[src_replica][worm_end_idx].next;
    
    // Check if worm tail is adjacent to a swap kink
    if (next_src!=-1){
        if (paths[src_replica][next_src].src_replica!=
            paths[src_replica][next_src].dest_replica){swap_in_front=true;}
        else if (paths[src_replica][prev_src].src_replica!=
                 paths[src_replica][prev_src].dest_replica){swap_in_front=false;}
        else {return;}
    }
    else{
        if (paths[src_replica][prev_src].src_replica!=
            paths[src_replica][prev_src].dest_replica){swap_in_front=false;}
        else {return;}
    }
                
    // Get index of central time slice at destination replica
    current_kink = worm_end_site; // next refers to index of next kink on site
    while (paths[dest_replica][current_kink].dest_replica==
           paths[dest_replica][current_kink].src_replica){
        current_kink = paths[dest_replica][current_kink].next;
    }
    kink_out_of_dest = current_kink;
    
    // Get indices of kinks before & after swap kink in destination replica
    prev_dest = paths[dest_replica][kink_out_of_dest].prev;
    next_dest = paths[dest_replica][kink_out_of_dest].next;
    
    // Determine the lower and upper bounds of the worm end to be timeshifted
    if (swap_in_front){
        if (next_dest!=-1)
            tau_next = paths[dest_replica][next_dest].tau;
        else
            tau_next = beta;
        tau_prev = paths[src_replica][prev_src].tau;
    }
    else{ // swap kink behind
        if (next_src!=-1)
            tau_next = paths[src_replica][next_src].tau;
        else
            tau_next=beta;
        tau_prev = paths[dest_replica][prev_dest].tau;
    }
    
    // Calculate change in diagonal energy
//    shift_tail=true; // we are always moving tail in this update. set to true.
//    dV=U*(n-!shift_head)-mu;
    dV=U*(n-1)-mu;
    
    // To make acceptance ratio unity,shift tail needs to sample w/ dV=eps-eps_w
    dV *= -1; // dV=eps-eps_w
    
    //boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    // Sample the new time of the worm end from truncated exponential dist.
    /*:::::::::::::::::::: Truncated Exponential RVS :::::::::::::::::::::::::*/
    Z = 1.0 - exp(-dV*(tau_next-tau_prev));
    tau_new = tau_prev - log(1.0-Z*rng.rand())  / dV;
    /*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
    
    // Check if advance or recede
    if (tau_new>tau){
        is_advance=true;
        swap_advance_tail_attempts+=1;
    }
    else {
        is_advance=false;
        swap_recede_tail_attempts+=1;
    }
    
    // Check if timeshift takes the worm end over the swap kink
    is_over_swap=false;
    if ( (is_advance && swap_in_front && tau_new>beta/2) ||
         (!is_advance && !swap_in_front && tau_new<beta/2) ){is_over_swap=true;}
        
    // Determine the length of path to be modified in each replica
    if (is_over_swap){
        l_path_src = beta/2 - tau;
        l_path_dest = tau_new - beta/2;
    }
    else{ // did not go over swap
        l_path_src = tau_new - tau;
        l_path_dest = 0;
    }
    
    // Determine the total particle change based on wormend to be shifted
    dN_src  = -1.0 * l_path_src/beta;
    dN_dest = -1.0 * l_path_dest/beta;
    
    // Canonical simulations: Restrict updates to interval N:(N-1,N+1)
    if (canonical){
        if ((N_tracker[src_replica]+dN_src)   < (N-1)  ||
            (N_tracker[src_replica]+dN_src)   > (N+1)  ||
            (N_tracker[dest_replica]+dN_dest) < (N-1)  ||
            (N_tracker[dest_replica]+dN_dest) > (N+1)){return;}
    }
    
    // Get number of particles after: worm end @ src & central kink @ dest
    n_after_worm_end = paths[src_replica][worm_end_idx].n;
    n_after_swap_kink = paths[dest_replica][kink_out_of_dest].n;
    n_before_swap_kink = paths[dest_replica][prev_dest].n;
        
    // Get number of kinks in source and destination replicas (before update)
    num_kinks_src = num_kinks[src_replica];
    num_kinks_dest = num_kinks[dest_replica];
        
    // Build the Metropolis condition (R)
    R = 1.0; // Sampling worm end time from truncated exponential makes R unity.

    // Metropolis sampling
    if (rng.rand() < R){
        if (!is_over_swap){ // worm end does not go over swap kink
            paths[src_replica][worm_end_idx].tau = tau_new;
            N_tracker[src_replica] += dN_src;
            if (is_advance){swap_advance_tail_accepts+=1;}
            else {swap_recede_tail_accepts+=1;}
        }
        else{ // We go Over Swap

            if (is_advance){ // advance OVER SWAP
                
                // Cannot advance tail over swap if no particles on destination
                if (n_after_swap_kink==0){return;}

                /*--------- Deletion of worm end from SOURCE replica ---------*/
                
                // num_kinks_src-1 will be swapped. Modify links to it
                if (paths[src_replica][num_kinks_src-1].next!=-1){
paths[src_replica][paths[src_replica][num_kinks_src-1].next].prev=worm_end_idx;
                }
paths[src_replica][paths[src_replica][num_kinks_src-1].prev].next=worm_end_idx;
                
                swap(paths[src_replica][worm_end_idx],
                     paths[src_replica][num_kinks_src-1]);
                
                // Upper or lower bound of flat could've been swapped. Correct.
                if (next_src==num_kinks_src-1){next_src=worm_end_idx;}
                else if (prev_src==num_kinks_src-1){prev_src=worm_end_idx;}
                else {;}
                
                // Head could've been swapped. Correct if so.
                if (head_idx[src_replica]==num_kinks_src-1)
                    head_idx[src_replica]=worm_end_idx;
                
                // Whatever kink was swapped could've been the last on its site
                if (paths[src_replica][worm_end_idx].next==-1){
                  last_kinks[src_replica][paths[src_replica][worm_end_idx].src]
                    =worm_end_idx;
                }

                // Reconnect upper and lower bounds of the flat
                if (next_src!=-1)
                    paths[src_replica][next_src].prev = prev_src;
                paths[src_replica][prev_src].next = next_src;
                
                // Deactivate the worm end
                tail_idx[src_replica]=-1;
                                
                // Update trackers for: no. of active kinks,total particles
                num_kinks[src_replica]-=1;
                N_tracker[src_replica] += dN_src;

                /*------- Insertion of worm end in DESTINATION replica -------*/
                
                // Activate first available kink
                paths[dest_replica][num_kinks_dest]=
                Kink(tau_new,n_after_swap_kink,worm_end_site,worm_end_site,
                     kink_out_of_dest,next_dest,dest_replica,dest_replica);
                
                // Save tail index
                tail_idx[dest_replica]=num_kinks_dest;
                
                // Update number of particles after swap kink in dest replica
                paths[dest_replica][kink_out_of_dest].n=n_after_swap_kink-1;
                
                // Add to acceptance counter
                swap_advance_tail_accepts+=1;
                
                // Modify links of swap kink and next_dest kink
                if (next_dest!=-1)
                    paths[dest_replica][next_dest].prev=num_kinks_dest;
                paths[dest_replica][kink_out_of_dest].next=num_kinks_dest;
                
                // Update trackers for: no. of active kinks,total particles
                num_kinks[dest_replica] += 1;
                N_tracker[dest_replica] += dN_dest;

                // Created kink might be last on its site
                if (next_dest==-1)
                    last_kinks[dest_replica][worm_end_site]=
                    num_kinks_dest;
                                
            }
            else{ // Recede OVER SWAP
                
                /*--------- Deletion of worm end from SOURCE replica ---------*/
                // num_kinks_src-1 will be swapped. Modify links to it
                if (paths[src_replica][num_kinks_src-1].next!=-1){
paths[src_replica][paths[src_replica][num_kinks_src-1].next].prev=worm_end_idx;
                }
paths[src_replica][paths[src_replica][num_kinks_src-1].prev].next=worm_end_idx;
                
                swap(paths[src_replica][worm_end_idx],
                     paths[src_replica][num_kinks_src-1]);
                
                // Upper or lower bound of flat could've been swapped. Correct.
                if (next_src==num_kinks_src-1){next_src=worm_end_idx;}
                else if (prev_src==num_kinks_src-1){prev_src=worm_end_idx;}
                else {;}
                
                // Head could've been swapped. Correct if so.
                if (head_idx[src_replica]==num_kinks_src-1)
                    head_idx[src_replica]=worm_end_idx;
                
                // Whatever kink was swapped could've been the last on its site
                if (paths[src_replica][worm_end_idx].next==-1){
                  last_kinks[src_replica][paths[src_replica][worm_end_idx].src]
                    =worm_end_idx;
                }
                
                // Modify particles in prev_src
                paths[src_replica][prev_src].n=n_after_worm_end;

                // Reconnect upper and lower bounds of the flat
                if (next_src!=-1)
                    paths[src_replica][next_src].prev = prev_src;
                paths[src_replica][prev_src].next = next_src;
                
                // Deactivate the worm end
                tail_idx[src_replica]=-1;
                                
                // Update trackers for: no. of active kinks,total particles
                num_kinks[src_replica]-=1;
                N_tracker[src_replica] += dN_src;
                
                // Swap kink on src might be last kink on it's site
                if (next_src==-1){
                    last_kinks[src_replica][worm_end_site]=prev_src;
                }
                

                /*------- Insertion of worm end in DESTINATION replica -------*/

                // Activate first available kink
                paths[dest_replica][num_kinks_dest]=
                Kink(tau_new,n_before_swap_kink+1,worm_end_site,worm_end_site,
                     prev_dest,kink_out_of_dest,dest_replica,dest_replica);
                
                // Save tail index
                tail_idx[dest_replica]=num_kinks_dest;
                
                // Add to acceptance counter
                swap_recede_tail_accepts+=1;
                
                // Modify links to worm head in destination replica
                paths[dest_replica][kink_out_of_dest].prev=num_kinks_dest;
                paths[dest_replica][prev_dest].next=num_kinks_dest;
                
                // Update trackers for: no. of active kinks,total particles
                num_kinks[dest_replica] += 1;
                N_tracker[dest_replica] += dN_dest;
            }
        }
        return;
    }
    else // Reject
        return;
}
/*--------------------------------------------------------------------*/

/*------------------------------- Estimators ---------------------------------*/

// For diagonal estimators around a time slice, Fock State will be needed
void get_fock_state(double measurement_center, int M,
                    vector<int> &fock_state_at_slice,
                    vector<Kink> &paths){
    
    double tau;
    int current,n_i;
    
    for (int i=0; i<M; i++){
        current=i;
        tau = paths[current].tau;
        while (tau<measurement_center+1.0E-12 && current!=-1){
            n_i=paths[current].n;
            fock_state_at_slice[i]=n_i;
            
            current=paths[current].next;
            if (current!=-1)
                tau=paths[current].tau;
        }
    }
    return;
}

/*--------------------------------------------------------------------*/

vector<double> get_measurement_centers(double beta){
    
    double tau_center;
    vector<double> measurement_centers;
    int num_centers;
    
    num_centers=25; // use odd number please.
    tau_center=beta/(2*num_centers);
    for (int i=0; i<num_centers; i++){
        measurement_centers.push_back(tau_center);
        tau_center+=(beta/(num_centers));
    }
    return measurement_centers;
}

/*-------------------------------- Diagonal ----------------------------------*/

double pimc_diagonal_energy(vector<int> &fock_state_at_slice, int M,
                            bool canonical, double U, double mu){
    
    double diagonal_energy;
    int n_i;
        
    diagonal_energy=0.0;
    for (int i=0; i<M; i++){
        n_i = fock_state_at_slice[i];
        if (canonical)
            diagonal_energy += (U/2.0*n_i*(n_i-1));
        else
            diagonal_energy += (U/2.0*n_i*(n_i-1)-mu*n_i);
    }
    return diagonal_energy;
}

/*--------------------------------------------------------------------*/


/*---------TEST THIS!!!*----------*/
void tau_resolved_diagonal_energy(vector<Kink> &paths,
                                           int num_kinks, int M, bool canonical,
                                           double U, double mu, double beta,
                                           vector<double> &measurement_centers,
                                           vector<double> &tr_diagonal_energy){
    double tau,measurement_center;
    int current,n_i;
    
    for (int i=0; i<M; i++){
        current=i;
        tau=paths[current].tau;
        n_i=paths[current].n;
        for (size_t j=0; j<measurement_centers.size(); j++){
            measurement_center=measurement_centers[j];
            while (tau<=measurement_center && current!=-1){
                n_i=paths[current].n;
                
                current=paths[current].next;
                if (current!=-1)
                    tau=paths[current].tau;
            }
            if (canonical)
                tr_diagonal_energy[j]+=(U/2.0*n_i*(n_i-1.0));
            else
                tr_diagonal_energy[j]+=(U/2.0*n_i*(n_i-1.0)-mu*n_i);
        }
    }
    
    return;
}

/*------------------------------ Off-Diagonal --------------------------------*/

double pimc_kinetic_energy(vector<Kink> &paths, int num_kinks,
                      double measurement_center,double measurement_plus_minus,
                      int M, double t, double beta){
    
    int kinks_in_window=0;
    
    for (int k=0; k<num_kinks; k++){
        if (paths[k].tau>=measurement_center-measurement_plus_minus
            && paths[k].tau<=measurement_center+measurement_plus_minus){
            kinks_in_window+=1;
        }
    }
    
    return (-kinks_in_window/2.0)/(2.0*measurement_plus_minus);
}

/*--------------------------------------------------------------------*/

void tau_resolved_kinetic_energy(vector<Kink> &paths,
                                           int num_kinks, int M,
                                           double t, double beta,
                                           vector<double> &measurement_centers,
                                           vector<double> &tr_kinetic_energy){
    double tau,measurement_center,window_width;

    window_width = measurement_centers[2]-measurement_centers[1];
    
    for (int i=M; i<num_kinks; i++){ // Note: the tau=0 kinks not counted
        tau = paths[i].tau;

        for (size_t j=0; j<measurement_centers.size(); j++){
            measurement_center=measurement_centers[j];

            if (tau>=measurement_center-window_width/2.0 &&
                tau<measurement_center+window_width/2.0){
                // add kink to bin
                tr_kinetic_energy[j]+=(-1.0/(2.0*window_width));
                break;
            }
        }
    }
    return;
}
/*--------------------------------------------------------------------*/

#endif /* pimc_hpp */
