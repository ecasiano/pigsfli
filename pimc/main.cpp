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

void insert_worm(vector<Kink> &kinks_vector, int &num_kinks, int &head_idx,
                 int &tail_idx, int M, int N, float U, float mu, float t,
                 float beta, float eta, bool canonical, double &N_tracker,
                 int &insert_worm_attempts, int &insert_worm_accepts,
                 int &insert_anti_attempts, int &insert_anti_accepts){
    
    // Variable declarations
    int k,n,site,dir,prev,next,n_head,n_tail;
    double tau,tau_h,tau_t,tau_prev,tau_next,tau_flat,l_path,dN,dV,p_iw,p_dw,R;
    bool is_worm;
    
    // Can only perform update if there are no worm ends
    if (head_idx != -1 || tail_idx != -1){return;}
    
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
            head_idx = num_kinks + 1;
            tail_idx = num_kinks;
            
            // Add to Acceptance counter
            insert_worm_accepts += 1;
        }
        else{ // Antiworm
            kinks_vector[num_kinks]= Kink(tau_h,n_head,site,0,k,num_kinks+1);
            kinks_vector[num_kinks+1]= Kink(tau_t,n_tail,site,0,num_kinks,next);
            
            // Save indices of head & tail kinks
            head_idx = num_kinks;
            tail_idx = num_kinks+1;
            
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

void delete_worm(vector<Kink> &kinks_vector, int &num_kinks, int &head_idx,
                 int &tail_idx, int M, int N, float U, float mu, float t,
                 float beta, float eta, bool canonical, double &N_tracker,
                 int &delete_worm_attempts, int &delete_worm_accepts,
                 int &delete_anti_attempts, int &delete_anti_accepts){
    
    // Variable declarations
    int k,n,site,dir,prev,next,n_head,n_tail;
    int head_prev,head_next,tail_prev,tail_next;
    double tau,tau_h,tau_t,tau_prev,tau_next,tau_flat,l_path,dN,dV,p_iw,p_dw,R;
    bool is_worm;
    
    // Can only propose worm deletion if both worm ends are present
    if (head_idx == -1 || tail_idx == -1){return;}
    
    // Can only delete worm if wormends are on same flat interval
    if (kinks_vector[head_idx].next != kinks_vector[tail_idx].prev &&
        kinks_vector[tail_idx].next != kinks_vector[head_idx].prev)
        return;
        
    // Extract worm end attributes
    tau_h = kinks_vector[head_idx].tau; // Head attributes
    n_head = kinks_vector[head_idx].n;
    site = kinks_vector[head_idx].site;
    dir = kinks_vector[head_idx].dir;
    head_prev = kinks_vector[head_idx].prev;
    head_next = kinks_vector[head_idx].next;
    
    tau_t = kinks_vector[head_idx].tau; // Tail attributes
    n_tail = kinks_vector[head_idx].n;
    site = kinks_vector[head_idx].site;
    dir = kinks_vector[head_idx].dir;
    tail_prev = kinks_vector[head_idx].prev;
    tail_next = kinks_vector[head_idx].next;
    
    // Identify the type of worm
    if (tau_h > tau_t)
        is_worm = true;
    else
        is_worm = false; // antiworm
    
    // Identify lower and upper bound of flat interval where worm lives
    if(is_worm){
        tau_prev = kinks_vector[tail_prev].tau;
        if(kinks_vector[head_idx].next == -1)
            tau_next = beta;
        else
            tau_next = kinks_vector[int(head_next)].tau;
        n = kinks_vector[tail_prev].n; // particles originally in the flat
            }
    else{ // antiworm
        tau_prev = kinks_vector[head_prev].tau;
        if (kinks_vector[tail_idx].next == -1)
            tau_next = beta;
        else
            tau_next = kinks_vector[tail_next].tau;
        n = kinks_vector[head_prev].n;
            }
    
    // Calculate the length of the flat interval
    tau_flat = tau_next - tau_prev;
    
    // Determine length of modified path and particle change
    l_path = tau_h - tau_t;
    dN = l_path/beta;
    
    // Canonical simulaton: Restrict updates to interval  N: (N-1,N+1)
    if (canonical)
        if ((N_tracker+dN) <= (N-1) || (N_tracker+dN) >= (N+1)){return;}
    
    // Calculate the difference in diagonal energy dV = \epsilon_w - \epsilon
    dV = (U/2)*(n_tail*(n_tail-1)-n_head*(n_head-1)) - mu*(n_tail-n_head);
    
    // Build the Metropolis ratio (R)
    p_dw = 1;
    p_iw = 1;
    R = eta * eta * n_tail * exp(-dV*(tau_h-tau_t))* (p_dw/p_iw) *
    (num_kinks-2) * tau_flat * tau_flat;
    R = 1/R;
    
    // Metropolis sampling
    if (rnum(rng) < R){ // Accept
        // Reconnect the lower,upper bounds of the flat interval
        if (is_worm){
            kinks_vector[tail_prev].next = head_next;
            kinks_vector[head_next].prev = tail_prev;
        }
        else{ // antiworm
            kinks_vector[head_prev].next = tail_next;
            kinks_vector[tail_next].prev = head_prev;
        }
        
        // Swap head kink with last active kink
        kinks_vector[head_idx] = kinks_vector[num_kinks-1]
        kinks_vector[num_kinks-1] = kinks_vector[head_idx]
        prv = int(kinks[head_idx[0]][4])
        nxt = int(kinks[head_idx[0]][5])
        kinks[prv][5] = head_idx[0]
        kinks[nxt][4] = head_idx[0]
        
        // Swap tail kink with second-to-last active kink
        kinks[tail_idx],kinks[num_kinks-2] = kinks[num_kinks-2],kinks[tail_idx]
        prv = int(kinks[tail_idx[0]][4])
        nxt = int(kinks[tail_idx[0]][5])
        kinks[prv][5] = tail_idx[0]
        kinks[nxt][4] = tail_idx[0]
        
        // Reduce number of active kinks tracker
        num_kinks[0] -= 2;
        
        // Deactivate the head,tail indices
        head_idx[0],tail_idx[0] = -1,-1;
        
        return;
    }
        else // Reject
            return;
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
    int head_idx = -1, tail_idx = -1;
    
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
    insert_worm(kinks_vector,num_kinks,head_idx,tail_idx,
                M,N,U,mu,t,beta,eta,canonical,N_tracker,
                insert_worm_attempts,insert_worm_accepts,
                insert_anti_attempts,insert_anti_accepts);

    // Print out the data structure
    for (int i=0;i<8;i++){
        cout << kinks_vector[i] << endl;
    }
    
    // Print out the head and tail indices
    cout << "head_idx: " << head_idx << endl;
    cout << "tail_idx: " << tail_idx << endl;
    
    // Print out the N_tracker
    cout << "N_tracker: " << N_tracker << endl;
    
    // Print out number of active kinks
    cout << "num_kinks: " << num_kinks << endl;
    
    // Print out accept/reject statistics
    cout<<"Insert Worm: "<<insert_worm_accepts<<"/"<<insert_worm_attempts<<endl;
    cout<<"Insert Anti: "<<insert_anti_accepts<<"/"<<insert_anti_attempts<<endl;
    
    return 0;
}
