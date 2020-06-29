//
//  main.cpp
//  pimc
//
//  Created by ecasiano on 6/19/20.
//  Copyright © 2020 Emanuel Casiano-Díaz. All rights reserved.
//

// A complete working C++ program to
// demonstrate all insertion methods
#include<iostream>
#include<vector>
#include<boost/random.hpp>
#include<cmath>

using namespace std;


// Set the random number generator
boost::random::mt19937 rng;

class Kink
{
    public:
    // Attribute declarations
    float tau;
    int n,site,dir,prev,next;
    
    // Member function declarations (prototypes)
    Kink (float,int,int,int,int,int); // Kink constructor
    
    // Make "<<" a friend of the Kink class
    friend ostream& operator<<(ostream& os, const Kink& dt);
};

// Member function definitions
Kink::Kink (float a,int b,int c,int d,int e,int f){
    tau = a;
    n = b;
    site = c;
    dir = d;
    prev = e;
    next = f;
}

// Overload "<<" operator
ostream& operator<<(ostream& os, const Kink& dt)
{
    os << '<' << dt.tau << ',' << dt.n << ',' << dt.site << ','
    << dt.dir << ',' << dt.prev << ',' << dt.next << '>' << endl;
    
    return os;
}

// Function definitions
vector<unsigned int> random_boson_config(int M,int N){
    // Generates random Fock state of N bosons in M=L^D sites
    
    vector<unsigned int> alpha (M,0);
    int site;
    
    // Initialize the distribution object. Note: Support is fully closed [0,M-1]
    boost::random::uniform_int_distribution<> sites(0, M-1);
    
    // Randomly sprinkle the N particles among sites
    for (int n=1; n<N; n++){
        site = sites(rng);
        alpha[site] += 1;
    }
    return alpha;
}

// Main
int main(){
    
    // Create an uniform distribution with support: [0.0,1.0)
    boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    
    // Bose-Hubbard parameters
    int L = 4, D = 1, N = 4;
    float t = 1.0, U = 1.0, mu = 0.5;
    vector<unsigned int> alpha;
    int M = pow(L,D); // total sites
    
    // Simulation parameters
    float eta = 1.0, beta = 1.0;
    
    // Counters and trackers
    int num_kinks = M;
    int worm_head_idx = -1, worm_tail_idx = -1;

    // Initialize Kink object
    Kink kink (-1,-1,-1,-1,-1,-1);
    
    // Pre-allocate vector that will store all kinks in the worldline
    vector<Kink> kinks_vector (100,kink);
                
    // Initialize the first L^D kinks
    for (int i=0; i<L; i++){
        // Modify kink attributes
        kinks_vector[i].tau = 0;
        kinks_vector[i].n = 1;
        kinks_vector[i].site = i;
        kinks_vector[i].dir = 0;
        kinks_vector[i].prev = -1;
        kinks_vector[i].next = -1;
    }
    
    // Total number of lattice points
    cout << "total_sites: " << M << endl;
    
    // Print out the number of kinks
    cout << "num_kinks: " << num_kinks << endl;

    // Generate some random numbers to test the generator
    for (int i = 0; i < 10; ++i) {
        std::cout << rnum(rng) << "\n";
    }
    
    // Generate a random fock state
    alpha = random_boson_config(M,N);
    for (int i=0;i<M;i++){
        cout << alpha[i];
    }
    cout << endl;
    
    return 0;
}
