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
#include<chrono>

using namespace std;
using namespace std::chrono;

// Set the random number generator
boost::random::mt19937 rng;

// Time code execution
auto start = high_resolution_clock::now();

class Kink
{
    public:
    // Attribute declarations
    double tau;
    int n,site,dir,prev,next;
    
    // Member function declarations (prototypes)
    Kink (double,int,int,int,int,int); // Kink constructor
    
    // Make "<<" a friend of the Kink class
    friend ostream& operator<<(ostream& os, const Kink& dt);
};

// Member function definitions
Kink::Kink (double a,int b,int c,int d,int e,int f){
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
    << dt.dir << ',' << dt.prev << ',' << dt.next << '>';
    
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
    vector<Kink> kinks_vector(10000000,Kink(-1,-1,-1,-1,-1,-1));

    // Initialize the first M=L^D kinks
    for (int i=0; i<M; i++){
        kinks_vector[i] = Kink(0,alpha[i],i,0,-1,-1);
    }
    return kinks_vector;
}

/*************************************************************************************************/

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
            kinks_vector[num_kinks]=Kink(tau_t,n_tail,site,0,k,num_kinks+1);
            kinks_vector[num_kinks+1]=Kink(tau_h,n_head,site,0,num_kinks,next);
            
            // Save indices of head & tail kinks
            head_idx = num_kinks + 1;
            tail_idx = num_kinks;
            
            // Add to Acceptance counter
            insert_worm_accepts += 1;
        }
        else{ // Antiworm
            kinks_vector[num_kinks]=Kink(tau_h,n_head,site,0,k,num_kinks+1);
            kinks_vector[num_kinks+1]=Kink(tau_t,n_tail,site,0,num_kinks,next);
            
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
        
        // Update trackers for: no of active kinks, total particles
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
    int prev_h,next_h,prev_t,next_t;
    double tau,tau_h,tau_t,tau_prev,tau_next,tau_flat,l_path,dN,dV,p_iw,p_dw,R;
    bool is_worm;
    
    // Can only propose worm deletion if both worm ends are present
    if (head_idx == -1 || tail_idx == -1){return;}
    
    // Can only delete worm if wormends are on same flat interval
    if (head_idx != kinks_vector[tail_idx].prev &&
        tail_idx != kinks_vector[head_idx].prev)
        return;
    
    // Extract worm end attributes
    tau_h = kinks_vector[head_idx].tau; // Head attributes
    n_head = kinks_vector[head_idx].n;
    site = kinks_vector[head_idx].site;
    dir = kinks_vector[head_idx].dir;
    prev_h = kinks_vector[head_idx].prev;
    next_h = kinks_vector[head_idx].next;
    
    tau_t = kinks_vector[tail_idx].tau; // Tail attributes
    n_tail = kinks_vector[tail_idx].n;
    site = kinks_vector[tail_idx].site;
    dir = kinks_vector[tail_idx].dir;
    prev_t = kinks_vector[tail_idx].prev;
    next_t = kinks_vector[tail_idx].next;

    // Identify the type of worm
    if (tau_h > tau_t)
        is_worm = true;
    else
        is_worm = false; // antiworm
    
    // Identify lower and upper bound of flat interval where worm lives
    if(is_worm){
        tau_prev = kinks_vector[prev_t].tau;
        delete_worm_attempts += 1; // Attempts counter
        if(kinks_vector[head_idx].next == -1)
            tau_next = beta;
        else
            tau_next = kinks_vector[int(next_h)].tau;
        n = kinks_vector[prev_t].n; // particles originally in the flat
            }
    else{ // antiworm
        tau_prev = kinks_vector[prev_h].tau;
        delete_anti_attempts += 1; // Attempts counter
        if (kinks_vector[tail_idx].next == -1)
            tau_next = beta;
        else
            tau_next = kinks_vector[next_t].tau;
        n = kinks_vector[prev_h].n;
            }
    
    // Calculate the length of the flat interval
    tau_flat = tau_next - tau_prev;
    
    // Determine length of modified path and particle change
    l_path = tau_h - tau_t;
    dN = -l_path/beta;
    
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
    boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    if (rnum(rng) < R){ // Accept
        
        // num_kinks-1,num_kinks-2 will be swapped. Modify links to these.
        kinks_vector[kinks_vector[num_kinks-1].next].prev = head_idx;
        kinks_vector[kinks_vector[num_kinks-1].prev].next = head_idx;
                
        swap(kinks_vector[head_idx],kinks_vector[num_kinks-1]);

        kinks_vector[kinks_vector[num_kinks-2].next].prev = tail_idx;
        kinks_vector[kinks_vector[num_kinks-2].prev].next = tail_idx;
        
        swap(kinks_vector[tail_idx],kinks_vector[num_kinks-2]);

        // Upper,lower bounds of flat could've been swapped. Correct if so.
        if (prev_h==num_kinks-1){prev_h=head_idx;}
        if (next_h==num_kinks-1){next_h=head_idx;}

        if (prev_t==num_kinks-2){prev_t=tail_idx;}
        if (next_t==num_kinks-2){next_t=tail_idx;}

        // Reconnect the lower,upper bounds of the flat interval
        if (is_worm){
            kinks_vector[prev_t].next = next_h;
            kinks_vector[next_h].prev = prev_t;
            
            // Add to Acceptance counter
            delete_worm_accepts += 1;
        }
        else{ // antiworm
            kinks_vector[prev_h].next = next_t;
            kinks_vector[next_t].prev = prev_h;
            
            // Add to Acceptance counter
            delete_anti_accepts += 1;
        }
        
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

/*************************************************************************************************/

void insertZero(vector<Kink> &kinks_vector, int &num_kinks, int &head_idx,
                int &tail_idx, int M, int N, float U, float mu, float t,
                float beta, float eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta,
                int &insertZero_worm_attempts, int &insertZero_worm_accepts,
                int &insertZero_anti_attempts, int &insertZero_anti_accepts){
    
    // Variable declarations
    int k,n,site,dir,prev,next,n_head,n_tail,i,N_b;
    double tau,tau_h,tau_t,tau_prev,tau_next,tau_flat,
    l_path,dN,dV,p_iw,p_dw,R,p_type,tau_new,p_wormend,C,W,p_dz,p_iz;
    bool is_worm;

    // Cannot insert if there's two worm ends present
    if (head_idx != -1 and tail_idx != -1){return;}
        
    // Randomly select site on which to insert worm/antiworm from tau=0
    boost::random::uniform_int_distribution<> sites(0, M-1);
    i = sites(rng);
    
    // Extract the flat interval where insertion is proposed & its attributes
    Kink insertion_flat = kinks_vector[i];
    // tau = insertion_flat.tau; tau is just zero
    n = insertion_flat.n;
    site = insertion_flat.site;
    dir = insertion_flat.dir;
    prev = insertion_flat.prev;
    next = insertion_flat.next;
    
    // Determine the length of insertion flat interval
    if (next != -1)
        tau_flat = kinks_vector[next].tau;
    else
        tau_flat = beta;
    
    // Choose worm/antiworm insertion based on worm ends present
    boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    if (head_idx==-1 and tail_idx==-1){ // no worm ends present
        if (n==0){ // can only insert worm, not antiworm
            is_worm = true;
            p_type = 1;
        }
        else{ // choose worm or antiworm insertion with equal probability
            if (rnum(rng) < 0.5)
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
            p_type = 1;
        }
    }
    else{ // only tail present, can insert worm only
        is_worm = true;
        p_type = 1;
    }
    
    // Add to worm/antiworm insertion attempt counters
    if (is_worm){insertZero_worm_attempts += 1;}
    else {insertZero_anti_attempts += 1;}
    
    // Randomly choose where to insert worm end on the flat interval
    tau_new = tau_flat*rnum(rng);
    
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
    dV = (U/2)*(n_tail*(n_tail-1)-n_head*(n_head-1)) - mu*(n_tail-n_head);
    
    // deleteZero (reverse update) might've had to choose either head or tail
    if (head_idx==-1 and tail_idx==-1) // only one end after insertZero
        p_wormend = 1;
    else{                              // two worm ends after insertZero
        if (is_worm){
            if (kinks_vector[kinks_vector[tail_idx].prev].tau != 0)
                p_wormend = 1; // cannot choose tail. was not coming from tau=0.
            else
                p_wormend = 0.5; // deleteZero could choose either head or tail
        }
        else{ // if insert anti (i.e, a tail) the end present was a head
            if (kinks_vector[kinks_vector[head_idx].prev].tau != 0)
                p_wormend = 1;
            else
                p_wormend = 0.5;
        }
    }
    
    // Determine the length of the path to be modified
    l_path = tau_new;
    
    // Determine the total particle change based on worm type
    if (is_worm){
        dN = +1 * l_path/beta;
    }
    else{
        dN = -1 * l_path/beta;
    }
    
    // Canonical simulations: Restrict updates to interval N:(N-1,N+1)
    if (canonical)
        if ((N_tracker+dN) <= (N-1) || (N_tracker+dN) >= (N+1)){return;}
    
    // Count the TOTAL number of particles at tau=0
    N_b = N_zero;
    
    // Build the weight ratio W'/W
    // C = 1;
    if (is_worm){
        C = sqrt(N_b+1)/sqrt(n+1);
        W = eta * sqrt(n_tail) * C * exp(-dV*tau_new);
    }
    else{
        C = sqrt(n)/sqrt(N_b);
        W = eta * sqrt(n_tail) * C * exp(dV*tau_new);
    }
    
    // Build the Metropolis Ratio (R)
    p_dz = 0.5;
    p_iz = 0.5;
    R = W * (p_dz/p_iz) * M * p_wormend * tau_flat / p_type;

    // Metropolis sampling
    if (rnum(rng) < R){ // Accept
        
        // Activate the first available kink
        if (is_worm){
            kinks_vector[num_kinks] = Kink (tau_new,n_head,i,0,i,next);
            
            // Save head index
            head_idx = num_kinks;
            
            // Update the number of particles in the initial kink at site i
            kinks_vector[i].n = n_tail;
            
            // Add to Acceptance counter
            insertZero_worm_accepts += 1;
        }
        else{ // antiworm
            kinks_vector[num_kinks] = Kink (tau_new,n_tail,i,0,i,next);
            
            // Save head index
            tail_idx = num_kinks;
            
            // Update the number of particles in the initial kink at site i
            kinks_vector[i].n = n_head;
            
            // Add to Acceptance counter
            insertZero_anti_accepts += 1;
        }
        
        // "Connect" next of lower bound kink to the new worm end
        kinks_vector[i].next = num_kinks;
        
        // "Connect" prev of next kink to the new worm end
        kinks_vector[next].prev = num_kinks;
        
        // Update trackers for: no. of active kinks, total particles
        num_kinks += 1;
        N_tracker += dN;
        
        return;
    }
    else // Reject
        return;
}

/*************************************************************************************************/

void deleteZero(vector<Kink> &kinks_vector, int &num_kinks, int &head_idx,
                int &tail_idx, int M, int N, float U, float mu, float t,
                float beta, float eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta,
                int &deleteZero_worm_attempts, int &deleteZero_worm_accepts,
                int &deleteZero_anti_attempts, int &deleteZero_anti_accepts){
    
    // Variable declarations
    int k,n,site,dir,prev,next,n_head,n_tail,i,N_b,worm_end_idx;
    double tau,tau_h,tau_t,tau_prev,tau_next,tau_flat,
    l_path,dN,dV,p_iw,p_dw,R,p_type,tau_new,p_wormend,C,W,p_dz,p_iz;
    bool is_worm,delete_head;

    // Cannot delete if there are no worm ends present
    if (head_idx==-1 and tail_idx==-1){return;}

    // Cannot delete if there are no worm ends coming from tau=0
    if (head_idx!=-1 and tail_idx!=-1){
        if (kinks_vector[kinks_vector[head_idx].prev].tau != 0 and
            kinks_vector[kinks_vector[tail_idx].prev].tau != 0){return;}
    }
    else if (head_idx!=-1){ // only head present
        if (kinks_vector[kinks_vector[head_idx].prev].tau != 0){return;}
    }
    else{ // only tail present
        if (kinks_vector[kinks_vector[head_idx].prev].tau != 0){return;}
    }

    // Decide which worm end to delete
    boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    if (head_idx==-1 and tail_idx==-1){ // both wormends present
        if (kinks_vector[kinks_vector[head_idx].prev].tau == 0 and
            kinks_vector[kinks_vector[tail_idx].prev].tau == 0){ //both near 0
            if (rnum(rng) < 0.5)
                delete_head = true;
            else
                delete_head = false;
            p_wormend = 0.5;
        }
        else if (kinks_vector[kinks_vector[head_idx].prev].tau == 0){
            delete_head = true;
            p_wormend = 1;
        }
        else{ // only the tail is near zero
            delete_head = false;
            p_wormend = 1;
        }
    }
    else if (head_idx!=-1){ // only head present
        delete_head = true;
        p_wormend = 1;
    }
    else{ // only tail present
        delete_head = false;
        p_wormend = 1;
    }

    // Get index of worm end to be deleted
    if (delete_head==true)
        worm_end_idx = head_idx;
    else
        worm_end_idx = tail_idx;

    // Extract worm end attributes
    tau = kinks_vector[worm_end_idx].tau;
    n = kinks_vector[worm_end_idx].n;
    site = kinks_vector[worm_end_idx].site;
    dir = kinks_vector[worm_end_idx].dir;
    prev = kinks_vector[worm_end_idx].prev;
    next = kinks_vector[worm_end_idx].next;

    // Calculate the length of the flat interval
    if (next == -1)
        tau_next = beta;
    else
        tau_next = kinks_vector[next].tau;
    tau_prev = kinks_vector[prev].prev;
    tau_flat = tau_next - tau_prev;

    // No. of particles before,after the worm end to be deleted
    if (delete_head==true) // delete worm
        n_tail = n+1;
    else //  delete antiworm
        n_tail = n;
    n_head = n_tail-1;
    
    // Worm insert (reverse update) probability of choosing worm or antiworm
    if (head_idx!=-1 and tail_idx!=-1) // worm end present before insertion
        p_type = 1;
    else{ // no worm ends present before insertion
        if (n==0)
            p_type = 1; // only worm can be inserted if no particles on flat
        else
            p_type = 0.5;
    }

    // Add to deleteZero PROPOSAL counters
    if (delete_head==true) // delete head (delete worm)
        deleteZero_worm_attempts += 1;
    else                   // delete tail (delete antiworm)
        deleteZero_anti_attempts += 1;
    
    // Determine the length of path to be modified
    l_path = tau;
    
    // Determine the total particle change based on worm type
    if (delete_head==true) // delete worm
        dN = -1 * l_path / beta;
    else // delete anti
        dN = +1 * l_path / beta;
    
    // Canonical simulaton: Restrict updates to interval  N: (N-1,N+1)
    if (canonical)
        if ((N_tracker+dN) <= (N-1) || (N_tracker+dN) >= (N+1)){return;}
    
    // Calculate diagonal energy difference
    dV = (U/2)*(n_tail*(n_tail-1)-n_head*(n_head-1)) - mu*(n_tail-n_head);
    
    // Determine the number of total bosons before the worm/anti was inserted
    if (delete_head)   // delete worm
        N_b = N_zero-1;
    else               // delete antiworm
        N_b = N_zero+1;
    
    // Build the weigh ratio W'/W
    // C = 1;
    if (delete_head){ // delete worm
        C = sqrt(N_b+1)/sqrt(n+1);
        W = eta * sqrt(n_tail) * C * exp(-dV*tau);
    }
    else{ // delete antiworm
        C = sqrt(n)/sqrt(N_b);
        W = eta * sqrt(n_tail) * C * exp(dV*tau);
    }
    
    // Build the Metropolis Ratio  (R)
    p_dz = 0.5;
    p_iz = 0.5;
    R = W * (p_dz/p_iz) * M * p_wormend * tau_flat / p_type;
    R = 1/R;
    
    // Metropolis sampling
    if (rnum(rng) < R){ // accept
        return;
    }
    else // reject
        return;
}
/*************************************************************************************************/

void timeshift(vector<Kink> &kinks_vector, int &num_kinks, int &head_idx,
int &tail_idx, int M, int N, float U, float mu, float t,
float beta, float eta, bool canonical, double &N_tracker,
int &delete_worm_attempts, int &delete_worm_accepts,
               int &delete_anti_attempts, int &delete_anti_accepts){
    
    // Variable declarations
    int k,n,site,dir,prev,next,n_head,n_tail;
    int prev_h,next_h,prev_t,next_t;
    double tau,tau_h,tau_t,tau_prev,tau_next,tau_flat,l_path,dN,dV,p_iw,p_dw,R;
    bool is_worm;
    
    // Metropolis sampling
    boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    if (rnum(rng) < R){
        return;
    }
    else // Reject
        return;
}

/*************************************************************************************************/

//class new_finite_exponential_random {
//public:
//  new_finite_exponential_random(double x, double l, double t)
//   : random_x(x), lambda(l), time(t) {}
//
//  double operator()() {
//    if(fabs(lambda*time)<1e-8)
//      return time*random_x;
//    else if(lambda<0)
//    {
//      if(lambda*time > std::log(std::numeric_limits<double>::min()))
//      {
//        double factor=std::exp(lambda*time);
//        return 1./lambda*std::log(factor+(1.-factor)*random_x);
//      }
//      else
//      {
//        double t= 1./lambda*std::log(random_x);
//        return (t < time ? t : time-std::numeric_limits<double>::epsilon());
//       }
//    }
//    else
//    {
//      if(-lambda*time > std::log(std::numeric_limits<double>::min()))
//      {
//        double factor=std::exp(-lambda*time);
//        return time - 1./(-lambda)*std::log(factor+(1.-factor)*random_x);
//        }
//      else
//      {
//        double t = time -1./(-lambda)*std::log(random_x);
//        return (t>0 ? t : std::numeric_limits<double>::epsilon());
//      }
//    }
//  }
//
//private:
//  double random_x;
//  double lambda;
//  double time;
//};   // new_finite_exponential_random
//
//template <class RNG>
//class finite_exponential_random {
//public:
//  finite_exponential_random(RNG& r, double l, double t)
//   : rng(r,boost::uniform_real<>()), lambda(l), time(t) {}
//
//  double operator()() {
//#ifdef SIMPLE
//    return time*rng();
//#else
//  return new_finite_exponential_random(rng(), lambda, time)();
//#endif
//  }
//
//private:
// boost::variate_generator<RNG&, boost::uniform_real<> > rng;
// double lambda;
// double time;
//};   // finite_exponential_random

/*************************************************************************************************/

// Main
int main(){
    
//    auto start = high_resolution_clock::now();

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
    int N_zero=N, N_beta=N;
    
    // Attempt/Acceptance counters
    int insert_worm_attempts=0, insert_worm_accepts=0;
    int insert_anti_attempts=0, insert_anti_accepts=0;
    int delete_worm_attempts=0, delete_worm_accepts=0;
    int delete_anti_attempts=0, delete_anti_accepts=0;
    int insertZero_worm_attempts=0, insertZero_worm_accepts=0;
    int insertZero_anti_attempts=0, insertZero_anti_accepts=0;

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
//    insert_worm(kinks_vector,num_kinks,head_idx,tail_idx,
//                M,N,U,mu,t,beta,eta,canonical,N_tracker,
//                insert_worm_attempts,insert_worm_accepts,
//                insert_anti_attempts,insert_anti_accepts);

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

    // Print out the data structure
    for (int i=0;i<8;i++){
        cout << kinks_vector[i] << endl;
    }

    // Perform a delete_worm
//    delete_worm(kinks_vector,num_kinks,head_idx,tail_idx,
//                M,N,U,mu,t,beta,eta,canonical,N_tracker,
//                delete_worm_attempts,delete_worm_accepts,
//                delete_anti_attempts,delete_anti_accepts);

    // Print out the data structure
    cout << endl;
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
    cout<<"Delete Worm: "<<delete_worm_accepts<<"/"<<delete_worm_attempts<<endl;
    cout<<"Delete Anti: "<<delete_anti_accepts<<"/"<<delete_anti_attempts<<endl;
    
    // Perform a delete_worm
    insertZero(kinks_vector,num_kinks,head_idx,tail_idx,
               M,N,U,mu,t,beta,eta,canonical,N_tracker,
               N_zero,N_beta,
               insertZero_worm_attempts,insertZero_worm_accepts,
               insertZero_anti_attempts,insertZero_anti_accepts);

    // Print out the data structure
    cout << endl;
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
    cout<<"InsertZero Worm: "<<insertZero_worm_accepts<<"/"<<
                               insertZero_worm_attempts<<endl;
    cout<<"InsertZero Anti: "<<insertZero_anti_accepts<<"/"<<
                               insertZero_anti_attempts<<endl;
        
    /************************************* insert/delete benchmarking ******************************************/

//    // Perform many insert/deletes back-to-back
//    for (int i=0; i<10000000; i++){
//
//        // Perform an insert_worm
//        insert_worm(kinks_vector,num_kinks,head_idx,tail_idx,
//                    M,N,U,mu,t,beta,eta,canonical,N_tracker,
//                    insert_worm_attempts,insert_worm_accepts,
//                    insert_anti_attempts,insert_anti_accepts);
//
//        // Perform a delete_worm
//        delete_worm(kinks_vector,num_kinks,head_idx,tail_idx,
//                    M,N,U,mu,t,beta,eta,canonical,N_tracker,
//                    delete_worm_attempts,delete_worm_accepts,
//                    delete_anti_attempts,delete_anti_accepts);
//    }
//
//
//    // Print out the data structure
//    cout << endl;
//    for (int i=0;i<8;i++){
//        cout << kinks_vector[i] << endl;
//    }
//
//    // Print out the head and tail indices
//    cout << "head_idx: " << head_idx << endl;
//    cout << "tail_idx: " << tail_idx << endl;
//
//    // Print out the N_tracker
//    cout << "N_tracker: " << N_tracker << endl;
//
//    // Print out number of active kinks
//    cout << "num_kinks: " << num_kinks << endl;
//
//    // Print out accept/reject statistics
//    cout<<"Insert Worm: "<<insert_worm_accepts<<"/"<<insert_worm_attempts<<endl;
//    cout<<"Delete Worm: "<<delete_worm_accepts<<"/"<<delete_worm_attempts<<endl;
//    cout<<"Insert Anti: "<<insert_anti_accepts<<"/"<<insert_anti_attempts<<endl;
//    cout<<"Delete Anti: "<<delete_anti_accepts<<"/"<<delete_anti_attempts<<endl;
//
//    auto end = high_resolution_clock::now();
//
//    auto elapsed_time = duration_cast<nanoseconds>(end - start);
//    float duration = elapsed_time.count() * 1e-9;
//
//    //elapsed_time *= 1
//
//    cout << "Elapsed time: " << duration << " seconds" << endl;
    
    /*************************************************************************************************/


    
    return 0;

}
