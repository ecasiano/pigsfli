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
vector<int> random_boson_config(int M,int N){
    // Generates random Fock state of N bosons in M=L^D sites
    
    vector<int> alpha (M,0);
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

vector<Kink> create_kinks_vector(vector<int> &alpha, int M){

    // Pre-allocate kinks. Recall: (tau,n,site,dir,prev,next)
    vector<Kink> kinks_vector(100,Kink(-1,-1,-1,-1,-1,-1));

    // Initialize the first M=L^D kinks
    for (int i=0; i<M; i++){
        kinks_vector[i] = Kink(0,alpha[i],i,0,-1,-1);
    }
    return kinks_vector;
}

void insert_worm(vector<Kink> &kinks_vector, int &num_kinks, int &worm_head_idx,
                 int &worm_tail_idx, int M, int N, float U, float mu, float t,
                 float beta, float eta, bool canonical, double &N_tracker,
                 int &insert_worm_attempts, int &insert_worm_accepts,
                 int &insert_anti_attempts, int &insert_anti_accepts){
    
    // Variable declarations
    int k,n,site,dir,prev,next,n_head,n_tail;
    double tau,tau_h,tau_t,tau_prev,tau_next,tau_flat,l_path,dN,dV,p_iw,p_dw,R;
    bool is_worm;
    
    // Can only perform update if there are no worm ends
    if (worm_head_idx != -1 || worm_tail_idx != -1){return;}
    
    // Randomly sample a flat interval (or kink if you like)
    boost::random::uniform_int_distribution<> flats(0, num_kinks-1);
    k = flats(rng);
    
    // Extract the attributes of the kink at the bottom of the flat interval
    tau = kinks_vector[k].tau;
    n = kinks_vector[k].n;
    site = kinks_vector[k].site;
    dir = kinks_vector[k].dir;
    prev = kinks_vector[k].prev;
    next = kinks_vector[k].next;
    
    // Calculate the length of the flat interval
    tau_prev = tau;
    if (next != -1) // tau_next extractable iff sampled kink is not the last
        tau_next = kinks_vector[next].tau;
    else
        tau_next = beta;
    tau_flat = tau_next - tau_prev;
    
    // Randomly choose where to insert worm ends in the flat interval
    boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    tau_h = tau_prev + tau_flat*rnum(rng);
    tau_t = tau_prev + tau_flat*rnum(rng);
    
    // Based on worm end time, determine worm type: antiworm or worm
    if (tau_h > tau_t){
        is_worm = true;
        insert_worm_accepts += 1; // Attempts counter
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
    if (n == 0 && !(is_worm)){return;}
    if (tau_h == tau_prev || tau_t == tau_prev){return;}
    if (tau_h == tau_t){return;}
    
    // Determine length of modified path and particle change
    l_path = tau_h - tau_t;
    dN = l_path/beta;
    
    // Canonical simulations: Restrict updates to interval N:(N-1,N+1)
    if (canonical)
        if ((N_tracker+dN) <= (N-1) || (N_tracker+dN) >= (N+1)){return;}
    
    // Calculate the difference in diagonal energy dV = \epsilon_w - \epsilon
    dV = (U/2)*(n_tail*(n_tail-1)-n_head*(n_head-1)) - mu*(n_tail-n_head);
    
    // Build the Metropolis ratio (R)
    p_dw = 1;
    p_iw = 1;
    R = eta * eta * n_tail * exp(-dV*(tau_h-tau_t))* (p_dw/p_iw) *
    num_kinks * tau_flat * tau_flat;

    // Metropolis sampling
    if (rnum(rng) < R){ // Accept

        // Activate the first two available kinks
        if (is_worm){
            kinks_vector[num_kinks]= Kink(tau_t,n_tail,site,0,k,num_kinks+1);
            kinks_vector[num_kinks+1]= Kink(tau_h,n_head,site,0,num_kinks,next);
            
            // Save indices of head & tail kinks
            worm_head_idx = num_kinks + 1;
            worm_tail_idx = num_kinks;
            
            // Add to Acceptance counter
            insert_worm_accepts += 1;
        }
        else{ // Antiworm
            kinks_vector[num_kinks]= Kink(tau_h,n_head,site,0,k,num_kinks+1);
            kinks_vector[num_kinks+1]= Kink(tau_t,n_tail,site,0,num_kinks,next);
            
            // Save indices of head & tail kinks
            worm_head_idx = num_kinks;
            worm_tail_idx = num_kinks+1;
            
            // Add to Acceptance counter
            insert_anti_accepts += 1;
        }
        
        // "Connect" next of lower bound kink to nearest worm end
        kinks_vector[k].next = num_kinks;
        
        // "Connect" prev of next kink to nearest worm end
        kinks_vector[next].prev = num_kinks+1;
        
        // Update trackers for: num of active kinks, total particles
        num_kinks += 2;
        N_tracker += dN;
        
        return;
    }
    else // Reject
        return;
}

void delete_worm(vector<Kink> &kinks_vector, int &num_kinks, int &worm_head_idx,
                 int &worm_tail_idx, int M, int N, float U, float mu, float t,
                 float beta, float eta, bool canonical, double &N_tracker,
                 int &delete_worm_attempts, int &delete_worm_accepts,
                 int &delete_anti_attempts, int &delete_anti_accepts){
}

// Main
int main(){
    
    // Create a uniform distribution with support: [0.0,1.0)
    boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    
    // Bose-Hubbard parameters
    int L = 4, D = 1, N = 4;
    float t = 1.0, U = 1.0, mu = 0.5;
    vector<int> alpha;
    int M = pow(L,D); // total sites
    
    // Simulation parameters
    float eta = 1.0, beta = 1.0;
    bool canonical = true;
    
    // Trackers
    int num_kinks = M;
    double N_tracker = N;
    int worm_head_idx = -1, worm_tail_idx = -1;
    
    // Attempt/Acceptance counters
    int insert_worm_attempts=0, insert_worm_accepts=0;
    int insert_anti_attempts=0, insert_anti_accepts=0;
    int delete_worm_attempts=0, delete_worm_accepts=0;
    int delete_anti_attempts=0, delete_anti_accepts=0;
    
    // Generate a random fock state
    alpha = random_boson_config(M,N);
    
    // Print out the generated random Fock state
    for (int i=0;i<M;i++){
        cout << alpha[i];
    }
    cout << endl;

    // Generate the data structure (vector of kinks)
    vector<Kink> kinks_vector = create_kinks_vector(alpha, M);

    // Print out the data structure
    for (int i=0;i<8;i++){
        cout << kinks_vector[i] << endl;
    }

    // Perform an insert_worm
    insert_worm(kinks_vector,num_kinks,worm_head_idx,worm_tail_idx,
                M,N,U,mu,t,beta,eta,canonical,N_tracker,
                insert_worm_attempts,insert_worm_accepts,
                insert_anti_attempts,insert_anti_accepts);

    // Print out the data structure
    for (int i=0;i<8;i++){
        cout << kinks_vector[i] << endl;
    }
    
    // Print out the head and tail indices
    cout << "worm_head_idx: " << worm_head_idx << endl;
    cout << "worm_tail_idx: " << worm_tail_idx << endl;
    
    // Print out the N_tracker
    cout << "N_tracker: " << N_tracker << endl;
    
    // Print out number of active kinks
    cout << "num_kinks: " << num_kinks << endl;
    
    // Print out accept/reject statistics
    cout<<"Insert Worm: "<<insert_worm_accepts<<"/"<<insert_worm_attempts<<endl;
    cout<<"Insert Anti: "<<insert_anti_accepts<<"/"<<insert_anti_attempts<<endl;
    
    
    return 0;
}
