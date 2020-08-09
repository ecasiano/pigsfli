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
//auto start = high_resolution_clock::now();

class Kink
{
    public:
    // Attribute declarations
    double tau;
    int n,src,dest,prev,next;
    
    // Member function declarations (prototypes)
    Kink (double,int,int,int,int,int); // Kink constructor
    
    // Make "<<" a friend of the Kink class
    friend ostream& operator<<(ostream& os, const Kink& dt);
};

// Member function definitions
Kink::Kink (double time,int particles,int source_site,int destination_site,
            int prev_kink_idx,int next_kink_idx){
    tau = time;
    n = particles;
    src = source_site;
    dest = destination_site;
    prev = prev_kink_idx;
    next = next_kink_idx;
//    pair = kink_pair_idx;
}

// Overload "<<" operator
ostream& operator<<(ostream& os, const Kink& dt)
{
    os << '<' << dt.tau << ',' << dt.n << ',' << dt.src << ','
    << dt.dest << ',' << dt.prev << ',' << dt.next << '>';
    
    return os;
}

/*-------------------------- Function Definitions ----------------------------*/

vector<int> random_boson_config(int M,int N){
    // Generates random Fock state of N bosons in M=L^D sites
    
    vector<int> alpha (M,0);
    int src;
    
    // Initialize the distribution object. Note: Support is fully closed [0,M-1]
    boost::random::uniform_int_distribution<> sites(0, M-1);
    
    // Randomly sprinkle the N particles among sites
    for (int n=1; n<=N; n++){
        src = sites(rng);
        alpha[src] += 1;
    }
    return alpha;
}

/*----------------------------------------------------------------------------*/

vector<Kink> create_kinks_vector(vector<int> &fock_state, int M){

    // Pre-allocate kinks. Recall: (tau,n,src,dest,prev,next)
    vector<Kink> kinks_vector(10000000,Kink(-1,-1,-1,-1,-1,-1));

    // Initialize the first M=L^D kinks
    for (int site=0; site<M; site++){
        kinks_vector[site] = Kink(0,fock_state[site],site,site,-1,-1);
    }
    return kinks_vector;
}

/*----------------------------------------------------------------------------*/

// Create function that calculates vector norm given array
double norm(vector<double> point){
    
    double squared_sum=0;
        
    for (int i=0; i<point.size(); i++){
        squared_sum += point[i]*point[i];
    }
    
    return sqrt(squared_sum);
}


/*----------------------------------------------------------------------------*/

int count_hypercube_nearest_neighbors(int L,int D,string boundary_condition){
    
    // NOTE: THESE OUTPUTS ARE VALID FOR PERIODIC BOUNDARY CONDITIONS.
    if (L>2){return D*2;}
    else {return D;}
}

/*----------------------------------------------------------------------------*/

void build_hypercube_adjacency_matrix(int L,int D, string boundary_condition,
                                      vector<vector<int>> &adjacency_matrix){
    
    int M = pow(L,D);
    int site_left,site_right,site_up,site_down,site_high,site_low;
    int top_row_ctr=L,bottom_row_ctr=0;
    
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
    
/*----------------------------------------------------------------------------*/


void build_adjacency_matrix(int L,int D,string boundary_condition,
                            vector<vector<bool>>&adjacency_matrix){

    // Variable declarations
    int M = pow(L,D); // Number of lattice points
    int ctr,a1,a2,a3;
    double r_NN;
    vector<double> points_difference (D,0);

    // Initialize normalized basis vectors
    vector<double> a1_vec {1,0,0};
    vector<double> a2_vec {0,1,0};
    vector<double> a3_vec {0,0,1};
    
    // Initialize array that will store all the points
    vector<double> empty_point (D,0);
    vector<vector<double>> points (M,empty_point);
    
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

/*----------------------------------------------------------------------------*/

void insert_worm(vector<Kink> &kinks_vector, int &num_kinks, int &head_idx,
                 int &tail_idx, int M, int N, float U, float mu, float t,
                 float beta, float eta, bool canonical, double &N_tracker,
                 int &N_zero, int &N_beta, vector<int> &last_kinks,
                 int &insert_worm_attempts, int &insert_worm_accepts,
                 int &insert_anti_attempts, int &insert_anti_accepts){
    
    // Variable declarations
    int k,n,src,dest,prev,next,n_head,n_tail;
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
    src = kinks_vector[k].src;
    dest = kinks_vector[k].dest;
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
            kinks_vector[num_kinks]=Kink(tau_t,n_tail,src,dest,
                                         k,num_kinks+1);
            kinks_vector[num_kinks+1]=Kink(tau_h,n_head,src,dest,
                                           num_kinks,next);
            
            // Save indices of head & tail kinks
            head_idx = num_kinks + 1;
            tail_idx = num_kinks;
            
            // Add to Acceptance counter
            insert_worm_accepts += 1;
        }
        else{ // Antiworm
            kinks_vector[num_kinks]=Kink(tau_h,n_head,src,dest,
                                         k,num_kinks+1);
            kinks_vector[num_kinks+1]=Kink(tau_t,n_tail,src,dest,
                                           num_kinks,next);
            
            // Save indices of head & tail kinks
            head_idx = num_kinks;
            tail_idx = num_kinks+1;
            
            // Add to Acceptance counter
            insert_anti_accepts += 1;
        }
        
        // "Connect" next of lower bound kink to nearest worm end
        kinks_vector[k].next = num_kinks;
        
        // "Connect" prev of next kink to nearest worm end
        if(next!=-1){kinks_vector[next].prev = num_kinks+1;}
        
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

/*----------------------------------------------------------------------------*/

void delete_worm(vector<Kink> &kinks_vector, int &num_kinks, int &head_idx,
                 int &tail_idx, int M, int N, float U, float mu, float t,
                 float beta, float eta, bool canonical, double &N_tracker,
                 int &N_zero, int &N_beta, vector<int> &last_kinks,
                 int &delete_worm_attempts, int &delete_worm_accepts,
                 int &delete_anti_attempts, int &delete_anti_accepts){
    
    // Variable declarations
    int k,n,src,dest,prev,next,n_head,n_tail;
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
    src = kinks_vector[head_idx].src;
    dest = kinks_vector[head_idx].dest;
    prev_h = kinks_vector[head_idx].prev;
    next_h = kinks_vector[head_idx].next;
    
    tau_t = kinks_vector[tail_idx].tau; // Tail attributes
    n_tail = kinks_vector[tail_idx].n;
    src = kinks_vector[tail_idx].src;
    dest = kinks_vector[tail_idx].dest;
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
        else if (next_h==num_kinks-1){next_h=head_idx;}
        else {;}

        if (prev_t==num_kinks-2){prev_t=tail_idx;}
        else if (next_t==num_kinks-2){next_t=tail_idx;}
        else {;}

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
        
        // Deactivate the deleted kinks
        kinks_vector[num_kinks-1].tau = -1;
        kinks_vector[num_kinks-1].n = -1;
        kinks_vector[num_kinks-1].src = -1;
        kinks_vector[num_kinks-1].dest = -1;
        kinks_vector[num_kinks-1].prev = -1;
        kinks_vector[num_kinks-1].next = -1;

        kinks_vector[num_kinks-2].tau = -1;
        kinks_vector[num_kinks-2].n = -1;
        kinks_vector[num_kinks-2].src = -1;
        kinks_vector[num_kinks-2].dest = -1;
        kinks_vector[num_kinks-2].prev = -1;
        kinks_vector[num_kinks-2].next = -1;
        
        // Update trackers for: num of active kinks, total particles
        num_kinks -= 2;
        N_tracker += dN;
        
        // If later worm end is last kink on site, update last kinks tracker vec
        if (is_worm){
            if (next_h==-1){last_kinks[src]=prev_t;}
        }
        else
            if (next_t==-1){last_kinks[src]=prev_h;}
        
        return;
    }
        
    else // Reject
            return;
}

/*----------------------------------------------------------------------------*/

void insertZero(vector<Kink> &kinks_vector, int &num_kinks, int &head_idx,
                int &tail_idx, int M, int N, float U, float mu, float t,
                float beta, float eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                int &insertZero_worm_attempts, int &insertZero_worm_accepts,
                int &insertZero_anti_attempts, int &insertZero_anti_accepts){
    
    // Variable declarations
    int k,n,src,dest,prev,next,n_head,n_tail,i,N_b;
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
    tau_prev = insertion_flat.tau; // tau is just zero
    n = insertion_flat.n;
    src = insertion_flat.src;
    dest = insertion_flat.dest;
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
            kinks_vector[num_kinks] = Kink(tau_new,n_head,src,dest,src,next);
            
            // Save head index
            head_idx = num_kinks;
            
            // Update the number of particles in the initial kink at site i
            kinks_vector[i].n = n_tail;
            
            // Add to Acceptance counter
            insertZero_worm_accepts += 1;
            
            // Worm inserted, add one to tau=0 particle tracker
            N_zero += 1;
        }
        else{ // antiworm
            kinks_vector[num_kinks] = Kink(tau_new,n_tail,src,dest,src,next);
            
            // Save head index
            tail_idx = num_kinks;
            
            // Update number of particles in initial kink of insertion site
            kinks_vector[i].n = n_head;
            
            // Add to Acceptance counter
            insertZero_anti_accepts += 1;
            
            // Antiworm inserted, subtract one to tau=0 particle tracker
            N_zero -= 1;
        }
        
        // "Connect" next of lower bound kink to the new worm end
        kinks_vector[src].next = num_kinks;
        
        // "Connect" prev of next kink to the new worm end
        if (next!=-1){kinks_vector[next].prev = num_kinks;}
        
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

/*----------------------------------------------------------------------------*/

void deleteZero(vector<Kink> &kinks_vector, int &num_kinks, int &head_idx,
                int &tail_idx, int M, int N, float U, float mu, float t,
                float beta, float eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                int &deleteZero_worm_attempts, int &deleteZero_worm_accepts,
                int &deleteZero_anti_attempts, int &deleteZero_anti_accepts){
    
    // Variable declarations
    int k,n,src,dest,prev,next,n_head,n_tail,i,N_b,worm_end_idx;
    double tau,tau_h,tau_t,tau_prev,tau_next,tau_flat,
    l_path,dN,dV,p_iw,p_dw,R,p_type,tau_new,p_wormend,C,W,p_dz,p_iz;
    bool is_worm,delete_head;

    // Cannot delete if there are no worm ends present
    if (head_idx==-1 && tail_idx==-1){return;}

    // Cannot delete if there are no worm ends coming from tau=beta
    if (head_idx!=-1 && tail_idx!=-1){
        if (kinks_vector[kinks_vector[head_idx].prev].tau != 0 and
            kinks_vector[kinks_vector[tail_idx].prev].tau != 0){return;}
    }
    else if (head_idx!=-1){ // only head present
        if (kinks_vector[kinks_vector[head_idx].prev].tau != 0){return;}
    }
    else{ // only tail present
        if (kinks_vector[kinks_vector[tail_idx].prev].tau != 0){return;}
    }

    // Decide which worm end to delete
    boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    if (head_idx!=-1 && tail_idx!=-1){ // both wormends present
        if (kinks_vector[kinks_vector[head_idx].prev].tau == 0 &&
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
    if (delete_head)
        worm_end_idx = head_idx;
    else
        worm_end_idx = tail_idx;

    // Extract worm end attributes
    tau = kinks_vector[worm_end_idx].tau;
    n = kinks_vector[worm_end_idx].n;
    src = kinks_vector[worm_end_idx].src;
    dest = kinks_vector[worm_end_idx].dest;
    prev = kinks_vector[worm_end_idx].prev;
    next = kinks_vector[worm_end_idx].next;

    // Calculate the length of the flat interval
    if (next == -1)
        tau_next = beta;
    else
        tau_next = kinks_vector[next].tau;
    tau_flat = tau_next;

    // No. of particles before,after the worm end to be deleted
    if (delete_head) // delete worm
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
    if (delete_head) // delete head (delete worm)
        deleteZero_worm_attempts += 1;
    else                   // delete tail (delete antiworm)
        deleteZero_anti_attempts += 1;
    
    // Determine the length of path to be modified
    l_path = tau;
    
    // Determine the total particle change based on worm type
    if (delete_head) // delete worm
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
        
        // Update the number of particles in the initial kink of worm end site
        kinks_vector[kinks_vector[worm_end_idx].prev].n = n;
        
        // num_kinks-1 (last available kink) will be swapped. Modify links to it
        if (kinks_vector[num_kinks-1].next!=-1) // avoids access with -1 index
            kinks_vector[kinks_vector[num_kinks-1].next].prev = worm_end_idx;
        kinks_vector[kinks_vector[num_kinks-1].prev].next = worm_end_idx;
        
        swap(kinks_vector[worm_end_idx],kinks_vector[num_kinks-1]);
        
        // Upper bound of flat could've been swapped. Correct if so.
        if (next==num_kinks-1){next=worm_end_idx;}
        
        // The other worm end could've been swapped. Reindex it if so.
        if (delete_head){
            
            if (tail_idx==num_kinks-1){tail_idx=worm_end_idx;}
            
            // If other end was last on its site, also reindex last_kinks array
            if (kinks_vector[tail_idx].next==-1){
                last_kinks[kinks_vector[tail_idx].src] = tail_idx;
            }
        }
        else{
            if (head_idx==num_kinks-1){head_idx=worm_end_idx;}
            
            // If other end was last on its site, also reindex last_kinks array
            if (kinks_vector[head_idx].next==-1){
                last_kinks[kinks_vector[head_idx].src] = head_idx;
            }
        }
        
        // Reconnect the lower,upper bounds of the flat interval.
        kinks_vector[prev].next = next;
        kinks_vector[next].prev = prev;

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
        
        // Deactivate the deleted kink
        kinks_vector[num_kinks-1].tau = -1;
        kinks_vector[num_kinks-1].n = -1;
        kinks_vector[num_kinks-1].src = -1;
        kinks_vector[num_kinks-1].dest = -1;
        kinks_vector[num_kinks-1].prev = -1;
        kinks_vector[num_kinks-1].next = -1;
        
        // Update trackers for: num of active kinks, total particles
        num_kinks -= 1;
        N_tracker += dN;
        
        // If worm end was last kink on site, update last kinks tracker vec.
        if (next==-1){last_kinks[src]=prev;}
        
        return;
    }
    else // reject
        return;
}
/*----------------------------------------------------------------------------*/

void insertBeta(vector<Kink> &kinks_vector, int &num_kinks, int &head_idx,
                int &tail_idx, int M, int N, float U, float mu, float t,
                float beta, float eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                int &insertBeta_worm_attempts, int &insertBeta_worm_accepts,
                int &insertBeta_anti_attempts, int &insertBeta_anti_accepts){
    
    // Variable declarations
    int k,n,src,dest,prev,next,n_head,n_tail,i,N_b;
    double tau,tau_h,tau_t,tau_prev,tau_next,tau_flat,
    l_path,dN,dV,p_iw,p_dw,R,p_type,tau_new,p_wormend,C,W,p_dz,p_iz,
    p_db,p_ib;
    bool is_worm;

    // Cannot insert if there's two worm ends present
    if (head_idx != -1 and tail_idx != -1){return;}
        
    // Randomly select site on which to insert worm/antiworm from tau=0
    boost::random::uniform_int_distribution<> sites(0, M-1);
    i = sites(rng);
    
    // Extract the flat interval where insertion is proposed & its attributes
    Kink insertion_flat = kinks_vector[last_kinks[i]];
    tau_prev = insertion_flat.tau;
    n = insertion_flat.n;
    src = insertion_flat.src;
    dest = insertion_flat.dest;
    prev = insertion_flat.prev;
    next = insertion_flat.next;
    
    // Determine the length of insertion flat interval
    tau_flat = beta - tau_prev;
    
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
    else if (tail_idx!=-1){ // only worm tail present, can insert antiworm only
        if (n==0){
            insertBeta_anti_attempts += 1;
            return; // cannot insert proposed antiworm, no particles present
        }
        else{
            is_worm = false;
            p_type = 1;
        }
    }
    else{ // only head present, can insert worm only
        is_worm = true;
        p_type = 1;
    }
    
    // Add to worm/antiworm insertion attempt counters
    if (is_worm){insertBeta_worm_attempts += 1;}
    else {insertBeta_anti_attempts += 1;}
    
    // Randomly choose where to insert worm end on the flat interval
    tau_new = tau_prev + tau_flat*rnum(rng);
    
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
            if (kinks_vector[kinks_vector[head_idx].next].tau != beta)
                p_wormend = 1; // cannot choose head. was not coming from beta.
            else
                p_wormend = 0.5; // deleteZero could choose either head or tail
        }
        else{ // if insert anti (i.e, a head) the end present was a tail
            if (kinks_vector[kinks_vector[tail_idx].next].tau != beta)
                p_wormend = 1;
            else
                p_wormend = 0.5;
        }
    }
    
    // Determine the length of the path to be modified
    l_path = beta - tau_new;
    
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
    N_b = N_beta;
    
    // Build the weight ratio W'/W
    // C = 1; // C_pre/C_post
    if (is_worm){
        C = sqrt(N_b+1)/sqrt(n+1);
        W = eta * sqrt(n_tail) * C * exp(-dV*(beta-tau_new));
    }
    else{
        C = sqrt(n)/sqrt(N_b);
        W = eta * sqrt(n_tail) * C * exp(-dV*(tau_new-beta));
    }
    
    // Build the Metropolis Ratio (R)
    p_db = 0.5;
    p_ib = 0.5;
    R = W * (p_db/p_ib) * M * p_wormend * tau_flat / p_type;

    // Metropolis sampling
    if (rnum(rng) < R){ // Accept
        
        // Activate the first available kink
        if (is_worm){
            kinks_vector[num_kinks] = Kink (tau_new,n_tail,src,dest,
                                            last_kinks[src],next);
            
            // Save tail index
            tail_idx = num_kinks;
            
            // Add to Acceptance counter
            insertBeta_worm_accepts += 1;
            
            // Worm inserted, add one to tau=beta particle tracker
            N_beta += 1;
        }
        else{ // antiworm
            kinks_vector[num_kinks] = Kink (tau_new,n_head,src,dest,
                                            last_kinks[src],next);
            
            // Save head index
            head_idx = num_kinks;
            
            // Add to Acceptance counter
            insertBeta_anti_accepts += 1;
            
            // Antiworm inserted, subtract one to tau=beta particle tracker
            N_beta -= 1;
        }
        
        // "Connect" next of lower bound kink to the new worm end
        kinks_vector[last_kinks[src]].next = num_kinks;
        
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

/*----------------------------------------------------------------------------*/

void deleteBeta(vector<Kink> &kinks_vector, int &num_kinks, int &head_idx,
                int &tail_idx, int M, int N, float U, float mu, float t,
                float beta, float eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                int &deleteBeta_worm_attempts, int &deleteBeta_worm_accepts,
                int &deleteBeta_anti_attempts, int &deleteBeta_anti_accepts){
    
    // Variable declarations
    int k,n,src,dest,prev,next,n_head,n_tail,i,N_b,worm_end_idx;
    double tau,tau_h,tau_t,tau_prev,tau_next,tau_flat,
    l_path,dN,dV,p_iw,p_dw,R,p_type,tau_new,p_wormend,C,W,p_dz,p_iz,
    p_db,p_ib;
    bool is_worm,delete_head;

    // Cannot delete if there are no worm ends present
    if (head_idx==-1 && tail_idx==-1){return;}

    // Cannot delete if there are no worm ends coming from tau=beta
    if (head_idx!=-1 && tail_idx!=-1){
        if (kinks_vector[head_idx].next != -1 and
            kinks_vector[tail_idx].next != -1){return;}
    }
    else if (head_idx!=-1){ // only head present
        if (kinks_vector[head_idx].next != -1){return;}
    }
    else{ // only tail present
        if (kinks_vector[tail_idx].next != -1){return;}
    }
    
    // Decide which worm end to delete
    boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    if (head_idx!=-1 && tail_idx!=-1){ // both wormends present
        if (kinks_vector[head_idx].next == -1 &&
            kinks_vector[tail_idx].next == -1){ // both last
            if (rnum(rng) < 0.5)
                delete_head = true;  // delete antiworm
            else
                delete_head = false; // delete worm
            p_wormend = 0.5;
        }
        else if (kinks_vector[head_idx].next == -1){
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
    if (delete_head)
        worm_end_idx = head_idx;
    else
        worm_end_idx = tail_idx;

    // Extract worm end attributes
    tau = kinks_vector[worm_end_idx].tau;
    n = kinks_vector[worm_end_idx].n;
    src = kinks_vector[worm_end_idx].src;
    dest = kinks_vector[worm_end_idx].dest;
    prev = kinks_vector[worm_end_idx].prev;
    next = kinks_vector[worm_end_idx].next;

    // Calculate the length of the flat interval
    tau_prev = kinks_vector[prev].tau;
    tau_flat = beta - tau_prev;

    // No. of particles before,after the worm end to be deleted
    if (delete_head) // delete antiworm
        n_tail = n+1;
    else //  delete worm
        n_tail = n;
    n_head = n_tail-1;
    
    // Worm insert (reverse update) probability of choosing worm or antiworm
    if (head_idx!=-1 and tail_idx!=-1) // worm end present before insertion
        p_type = 1;
    else{ // no worm ends present before insertion
        if (kinks_vector[kinks_vector[worm_end_idx].prev].n==0)
            p_type = 1; // only worm can be inserted if no particles on flat
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
        dN = +1 * l_path / beta;
    else // delete worm
        dN = -1 * l_path / beta;
    
    // Canonical simulaton: Restrict updates to interval  N: (N-1,N+1)
    if (canonical)
        if ((N_tracker+dN) <= (N-1) || (N_tracker+dN) >= (N+1)){return;}
    
    // Calculate diagonal energy difference
    dV = (U/2)*(n_tail*(n_tail-1)-n_head*(n_head-1)) - mu*(n_tail-n_head);
    
    // Determine the number of total bosons before the worm/anti was inserted
    if (delete_head)     // delete antiworm
        N_b = N_beta+1;
    else                 // delete worm
        N_b = N_beta-1;
    
    // Build the weigh ratio W'/W
    // C = 1;
    if (delete_head==false){ // delete worm
        C = sqrt(N_b+1)/sqrt(n);
        W = eta * sqrt(n_tail) * C * exp(-dV*(beta-tau));
    }
    else{ // delete antiworm
        C = sqrt(n+1)/sqrt(N_b);
        W = eta * sqrt(n_tail) * C * exp(-dV*(tau-beta));
    }
    
    // Build the Metropolis Ratio  (R)
    p_db = 0.5;
    p_ib = 0.5;
    R = W * (p_db/p_ib) * M * p_wormend * tau_flat / p_type;
    R = 1/R;
    
    // Metropolis sampling
    if (rnum(rng) < R){ // acclast
        
        // num_kinks-1 (last available kink) will be swapped. Modify links to it
        if (kinks_vector[num_kinks-1].next!=-1) // avoids access with -1 index
            kinks_vector[kinks_vector[num_kinks-1].next].prev = worm_end_idx;
        kinks_vector[kinks_vector[num_kinks-1].prev].next = worm_end_idx;
        
        swap(kinks_vector[worm_end_idx],kinks_vector[num_kinks-1]);
        
        // Lower bound of flat could've been swapped. Correct if so.
        if (prev==num_kinks-1){prev=worm_end_idx;}
        
        // The other worm end could've been swapped. Reindex it if so.
        if (delete_head){
            
            if (tail_idx==num_kinks-1){tail_idx=worm_end_idx;}
            
            // If other end was last on its site, also reindex last_kinks array
            if (kinks_vector[tail_idx].next==-1){
                last_kinks[kinks_vector[tail_idx].src] = tail_idx;
            }
        }
        else{
            if (head_idx==num_kinks-1){head_idx=worm_end_idx;}
            
            // If other end was last on its site, also reindex last_kinks array
            if (kinks_vector[head_idx].next==-1){
                last_kinks[kinks_vector[head_idx].src] = head_idx;
            }
        }
        
        // Reconnect the lower,upper bounds of the flat interval.
        kinks_vector[prev].next = next;
        kinks_vector[next].prev = prev;

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
        
        // Deactivate the deleted kink
        kinks_vector[num_kinks-1].tau = -1;
        kinks_vector[num_kinks-1].n = -1;
        kinks_vector[num_kinks-1].src = -1;
        kinks_vector[num_kinks-1].dest = -1;
        kinks_vector[num_kinks-1].prev = -1;
        kinks_vector[num_kinks-1].next = -1;
        
        // Update trackers for: num of active kinks, total particles
        num_kinks -= 1;
        N_tracker += dN;
        
        // Since worm end was last kink on site, update last kinks tracker vec.
        last_kinks[src]=prev;
        
        return;
    }
    else // reject
        return;
}
/*----------------------------------------------------------------------------*/

void timeshift_uniform(vector<Kink> &kinks_vector, int &num_kinks, int &head_idx,
                int &tail_idx, int M, int N, float U, float mu, float t,
                float beta, float eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                int &advance_head_attempts, int &advance_head_accepts,
                int &recede_head_attempts, int &recede_head_accepts,
                int &advance_tail_attempts, int &advance_tail_accepts,
                int &recede_tail_attempts, int &recede_tail_accepts){
    
    // Variable declarations
    int k,n,src,dest,prev,next,n_head,n_tail,i,N_b,worm_end_idx;
    double tau,tau_h,tau_t,tau_prev,tau_next,tau_flat,
    l_path,dN,dV,p_iw,p_dw,R,p_type,tau_new,p_wormend,C,W,p_dz,p_iz,
    p_db,p_ib;
    bool is_worm,delete_head,shift_head;
    
    // Reject update if there are is no worm end present
    if (head_idx==-1 && tail_idx==-1){return;}

    // Choose which worm end to move
    boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    if (head_idx!=-1 && tail_idx!=-1){ // both worm ends present
        tau_h = kinks_vector[head_idx].tau;
        tau_t = kinks_vector[tail_idx].tau;

        // Randomly choose to shift HEAD or TAIL
        if (rnum(rng) < 0.5)
            shift_head = true;
        else
            shift_head = false;
        }
    else if (head_idx!=-1){ // only head present
        tau_h = kinks_vector[head_idx].tau;
        shift_head = true;
    }
    else{ // only tail present
        tau_t = kinks_vector[tail_idx].tau;
        shift_head = false;
    }
    
    // Save the kink index of the end that will be shifted
    if (shift_head){worm_end_idx=head_idx;}
    else {worm_end_idx=tail_idx;}
    
    // Extract worm end attributes
    tau = kinks_vector[worm_end_idx].tau;
    n = kinks_vector[worm_end_idx].n;
    src = kinks_vector[worm_end_idx].src;
    dest = kinks_vector[worm_end_idx].dest;
    prev = kinks_vector[worm_end_idx].prev;
    next = kinks_vector[worm_end_idx].next;
    
    // Measure diagonal energy difference dV
    if (shift_head)
        dV=U*n-mu;
    else
        dV=U*(n-1)-mu;
    
    // Determine the lower and upper bounds of the worm end to be timeshifted
    if (next==-1)
        tau_next = beta;
    else
        tau_next = kinks_vector[next].tau;
    tau_prev = kinks_vector[prev].tau;
    
    // Calculate length of flat interval
    tau_flat = tau_next - tau_prev;

    // Sample the new time of the worm end from uniform distribution
    tau_new = tau_prev + tau_flat*rnum(rng);
    
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
        dN = +1 * l_path/beta;
    else // shift tail
        dN = -1 * l_path/beta;
    
    // Canonical simulaton: Restrict updates to interval  N: (N-1,N+1)
    if (canonical)
        if ((N_tracker+dN) <= (N-1) || (N_tracker+dN) >= (N+1)){return;}
    
    // Build the Metropolis condition (R)
    if (shift_head){
        W = exp(-dV*(tau_new-tau));
    }
    else{
        W = exp(-dV*(tau-tau_new));
    }
    R = W; // recall that timeshift is it's own inverse

    // Metropolis sampling
    if (rnum(rng) < R){
        
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
        kinks_vector[worm_end_idx].tau = tau_new;
        
        // Modify total particle number tracker
        N_tracker += dN;
        
        return;
    }
    else // Reject
        return;
}

/*----------------------------------------------------------------------------*/

void timeshift(vector<Kink> &kinks_vector, int &num_kinks, int &head_idx,
                int &tail_idx, int M, int N, float U, float mu, float t,
                float beta, float eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                int &advance_head_attempts, int &advance_head_accepts,
                int &recede_head_attempts, int &recede_head_accepts,
                int &advance_tail_attempts, int &advance_tail_accepts,
                int &recede_tail_attempts, int &recede_tail_accepts){
    
    // Variable declarations
    int k,n,src,dest,prev,next,n_head,n_tail,i,N_b,worm_end_idx;
    double tau,tau_h,tau_t,tau_prev,tau_next,tau_flat,tau_new,Z,
    l_path,dN,dV,p_iw,p_dw,R,p_type,p_wormend,C,W,p_dz,p_iz,
    p_db,p_ib;
    bool is_worm,delete_head,shift_head;
    
    // Reject update if there is no worm end present
    if (head_idx==-1 && tail_idx==-1){return;}

    // Choose which worm end to move
    boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    if (head_idx!=-1 && tail_idx!=-1){ // both worm ends present
        tau_h = kinks_vector[head_idx].tau;
        tau_t = kinks_vector[tail_idx].tau;

        // Randomly choose to shift HEAD or TAIL
        if (rnum(rng) < 0.5)
            shift_head = true;
        else
            shift_head = false;
        }
    else if (head_idx!=-1){ // only head present
        tau_h = kinks_vector[head_idx].tau;
        shift_head = true;
    }
    else{ // only tail present
        tau_t = kinks_vector[tail_idx].tau;
        shift_head = false;
    }
    
    // Save the kink index of the end that will be shifted
    if (shift_head){worm_end_idx=head_idx;}
    else {worm_end_idx=tail_idx;}
    
    // Extract worm end attributes
    tau = kinks_vector[worm_end_idx].tau;
    n = kinks_vector[worm_end_idx].n;
    src = kinks_vector[worm_end_idx].src;
    dest = kinks_vector[worm_end_idx].dest;
    prev = kinks_vector[worm_end_idx].prev;
    next = kinks_vector[worm_end_idx].next;
    
//    // Measure diagonal energy difference dV=eps_w-eps
//    if (shift_head)
//        dV=U*n-mu;     // dV=eps_w-eps
//    else
//        dV=U*(n-1)-mu; // dV=eps_w-eps
    
    dV=U*(n-!shift_head)-mu;
    
    // To make acceptance ratio unity,shift tail needs to sample w/ dV=eps-eps_w
    if (!shift_head){dV *= -1;} // dV=eps-eps_w
        
    // Determine the lower and upper bounds of the worm end to be timeshifted
    if (next==-1)
        tau_next = beta;
    else
        tau_next = kinks_vector[next].tau;
    tau_prev = kinks_vector[prev].tau;
    
    // Calculate length of flat interval
    tau_flat = tau_next - tau_prev;

    // Sample the new time of the worm end from truncated exponential dist.
    /*:::::::::::::::::::: Truncated Exponential RVS :::::::::::::::::::::::::*/
    Z = 1 - exp(-dV*(tau_next-tau_prev));
    tau_new = tau_prev - log(1-Z*rnum(rng))  / dV;
    /*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
    
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
        dN = +1 * l_path/beta;
    else // shift tail
        dN = -1 * l_path/beta;
    
    // Canonical simulaton: Restrict updates to interval  N: (N-1,N+1)
    if (canonical)
        if ((N_tracker+dN) <= (N-1) || (N_tracker+dN) >= (N+1)){return;}
    
    // Build the Metropolis condition (R)
    R = 1; // Sampling worm end time from truncated exponential makes R unity.

    // Metropolis sampling
    if (rnum(rng) < R){
        
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
        kinks_vector[worm_end_idx].tau = tau_new;
        
        // Modify total particle number tracker
        N_tracker += dN;
        
        return;
    }
    else // Reject
        return;
}

/*----------------------------------------------------------------------------*/

void insert_kink_before_head(vector<Kink> &kinks_vector, int &num_kinks,
                int &head_idx,int &tail_idx,
                int M, int N, float U, float mu, float t,
                vector<vector<int>> &adjacency_matrix, int total_nn,
                float beta, float eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                int &ikbh_attempts, int &ikbh_accepts){
    
    // Variable declarations
    int k,n,src,dest,prev,next,n_head,n_tail,i,N_b,worm_end_idx,j,
    n_i,n_wi,n_j,n_wj,prev_i,prev_j,next_i,next_j,ctr;
    double tau,tau_h,tau_t,tau_prev,tau_next,tau_flat,tau_new,Z,
    l_path,dN,dV,p_iw,p_dw,R,p_type,p_wormend,p_site,C,W,p_dz,p_iz,
    p_db,p_ib,p_dkbh,p_ikbh,tau_prev_i,tau_prev_j,tau_kink,tau_min,dV_i,dV_j;
    bool is_worm,delete_head,shift_head;
        
    // Update not possible if no worm head present
    if (head_idx==-1){return;}
    
    // Need at least two sites to perform a spaceshift
    if (M<2){return;}
        
    // Add to proposal counter
    ikbh_attempts += 1;
    
    // Extract the worm head site
    i = kinks_vector[head_idx].src;
    
    // Randomly choose a nearest neighbor site
    boost::random::uniform_int_distribution<> random_nn(0, total_nn-1);
    j = adjacency_matrix[i][random_nn(rng)];
    p_site = 1.0/total_nn;

    // Retrieve the time of the worm head
    tau_h = kinks_vector[head_idx].tau;
    
    // Determine index of lower/upper kinks of flat where head is (site i)
    prev_i = kinks_vector[head_idx].prev;
    next_i = kinks_vector[head_idx].next;
    
    // Determine index of lower/upper kinks of flat where head jumps to (site j)
    tau = 0;            // tau_prev_j candidate
    prev = j;           // prev_j candidate
    prev_j = j;         // this avoids "variable maybe not initialized" warning
    while (tau<tau_h){
//
//        cout << "HELP ME! I'm STUCK IN A WHILE LOOP AHHHH!!!!!!";
//        cout << tau << " " << tau_h << endl;
//        for (int g=0; g<M; g++){
//            cout << last_kinks[g] << " ";
//        }
//        cout << endl;
        
        // Set the lower bound index
        prev_j = prev;
        
        // Update lower bound index and tau candidates for next iteration
        prev = kinks_vector[prev].next;
        if (prev==-1){break;}
        tau = kinks_vector[prev].tau;
    }
    next_j=prev;
    
    // Determine tau_min by comparing lower bound candidates of src & dest sites
    tau_prev_i = kinks_vector[prev_i].tau;
    tau_prev_j = kinks_vector[prev_j].tau;
    if (tau_prev_i>tau_prev_j){tau_min=tau_prev_i;}
    else {tau_min=tau_prev_j;}
    
    // Extract no. of particles in the flats adjacent to the new kink
    n_wi = kinks_vector[prev_i].n;
    n_i = n_wi-1;
    n_j = kinks_vector[prev_j].n;
    n_wj = n_j+1;                   // "w": segment with the extra particle
    
    // Randomly choose the time of the kink
    boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    tau_kink = tau_min + rnum(rng)*(tau_h-tau_min);
    if (tau_kink == tau_min){return;}
        
    // Calculate the diagonal energy difference on both sites
    dV_i = (U/2.0)*(n_wi*(n_wi-1)-n_i*(n_i-1)) - mu*(n_wi-n_i);
    dV_j = (U/2.0)*(n_wj*(n_wj-1)-n_j*(n_j-1)) - mu*(n_wj-n_j);
    
    // Calculate the weight ratio W'/W
    W = t * n_wj * exp((dV_i-dV_j)*(tau_h-tau_kink));
    
    // Build the Metropolis ratio (R)
    p_dkbh = 0.5;
    p_ikbh = 0.5;
    R = W * (p_dkbh/p_ikbh) * (tau_h-tau_min)/p_site;
    
//    cout << endl << "R: " << R << endl;
//    cout << "p_site: " << p_site << endl;
//    cout << "total_nn: " << total_nn << endl;
//    cout << "tau_h: " << tau_h << endl;
//    cout << "tau_min: " << tau_min << endl;
//    cout << "W: " << W << endl;
//    cout << "t: " << t << endl;
//    cout << "dV_i: " << dV_i << endl;
//    cout << "dV_j: " << dV_j << endl;
//    cout << "prev_j: " << prev_j << endl;
//    cout << "next_j: " << next_j << endl;
//    cout << "prev_j: " << prev_i << endl;
//    cout << "next_j: " << next_i << endl;
//    cout << "i,j: " << i << "," << j << endl;
//    cout << "n_wi,n_j: " << n_wi << "," << n_j << endl;
//    cout << "---------------------------------"<<endl;
    
    // Metropolis Sampling
    if (rnum(rng) < R){ // Accept
        
        // Add to acceptance counter
        ikbh_accepts += 1;
                
        // Change kink that stored head information to a regular kink
        kinks_vector[head_idx].tau = tau_kink;
        kinks_vector[head_idx].n = n_i;
        kinks_vector[head_idx].src = i;
        kinks_vector[head_idx].dest = j;
        kinks_vector[head_idx].prev = prev_i;
        kinks_vector[head_idx].next = next_i;
//        kinks_vector[head_idx].pair = num_kinks;
        
        // Create the kinks on the destination site
        kinks_vector[num_kinks]=Kink(tau_kink,n_wj,j,i,prev_j,num_kinks+1);
        kinks_vector[num_kinks+1]=Kink(tau_h,n_j,j,j,num_kinks,next_j);
        
        // Set new worm head index
        head_idx = num_kinks+1;
                
        // "Connect" next of lower bound kink to new kink
        kinks_vector[prev_j].next = num_kinks;
        
        // "Connect" prev of next kink to worm head
        if(next_j!=-1){kinks_vector[next_j].prev = head_idx;}
                
        // Update number of kinks tracker
        num_kinks += 2;
        
        // If worm head is last kink on site j, update last kinks tracker vector
        if (next_j==-1){last_kinks[j]=head_idx;}
        
        return;
        
    }
    else // Reject
        return;
    }

/*----------------------------------------------------------------------------*/

void delete_kink_before_head(vector<Kink> &kinks_vector, int &num_kinks,
                int &head_idx,int &tail_idx,
                int M, int N, float U, float mu, float t,
                vector<vector<int>> &adjacency_matrix, int total_nn,
                float beta, float eta, bool canonical, double &N_tracker,
                int &N_zero, int &N_beta, vector<int> &last_kinks,
                int &dkbh_attempts, int &dkbh_accepts){

    // Variable declarations
    int k,n,src,dest,prev,next,n_head,n_tail,i,N_b,worm_end_idx,j,
    n_i,n_wi,n_j,n_wj,prev_i,prev_j,next_i,next_j,ctr,kink_idx_i,kink_idx_j;
    double tau,tau_h,tau_t,tau_prev,tau_next,tau_flat,tau_new,Z,
    l_path,dN,dV,p_iw,p_dw,R,p_type,p_wormend,p_site,C,W,p_dz,p_iz,
    p_db,p_ib,p_dkbh,p_ikbh,tau_prev_i,tau_prev_j,tau_kink,tau_min,dV_i,dV_j,
    tau_next_i;
    bool is_worm,delete_head,shift_head;

    // Update not possible if no worm head present
    if (head_idx==-1){return;}
    
    if (kinks_vector[head_idx].prev==tail_idx ||
        kinks_vector[kinks_vector[head_idx].prev].tau==0){return;}

    // Need at least two sites to perform a spaceshift
    if (M<2){return;}

    // Indices of: upper bound kink, kink before head, lower bound kink ; site j
    next_j = kinks_vector[head_idx].next;
    kink_idx_j = kinks_vector[head_idx].prev;
    prev_j = kinks_vector[kink_idx_j].prev;

    // Times of: worm head, kink before head, lower bound kink
    tau_h = kinks_vector[head_idx].tau;
    tau_kink = kinks_vector[kink_idx_j].tau;
    tau_prev_j = kinks_vector[prev_j].tau;

    // Retrieve worm head site (j) and connecting site (i)
    j = kinks_vector[kink_idx_j].src;
    i = kinks_vector[kink_idx_j].dest;

    // Determine index of lower/upper bounds of flat where kink connects to
    tau = 0;            // tau_prev_i candidate
    prev = i;           // prev_i candidate
    prev_i = i;         // this avoids "variable maybe not initialized" warning
    while (tau<tau_kink){
        // Set the lower bound index
        prev_i = prev;

        // Update lower bound index and tau candidates for next iteration
        prev = kinks_vector[prev].next;
        if (prev==-1){break;}
        tau = kinks_vector[prev].tau;
    }
    kink_idx_i = prev;
    next_i=kinks_vector[kink_idx_i].next;

    // Retrieve time of lower,upper bounds on connecting site (i)
    tau_prev_i = kinks_vector[prev_i].tau;
    if (next_i==-1){tau_next_i=beta;}
    else {tau_next_i = kinks_vector[next_i].tau;};
    
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
    n_wi = kinks_vector[prev_i].n;
    n_i = n_wi-1;
    n_j = kinks_vector[prev_j].n;
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
    R = 1/R;
    
//    cout << endl << "R: " << R << endl;
//    cout << "p_site: " << p_site << endl;
//    cout << "total_nn: " << total_nn << endl;
//    cout << "tau_h: " << tau_h << endl;
//    cout << "tau_min: " << tau_min << endl;
//    cout << "W: " << W << endl;
//    cout << "t: " << t << endl;
//    cout << "dV_i: " << dV_i << endl;
//    cout << "dV_j: " << dV_j << endl;
//    cout << "prev_j: " << prev_j << endl;
//    cout << "next_j: " << next_j << endl;
//    cout << "prev_j: " << prev_i << endl;
//    cout << "next_j: " << next_i << endl;
//    cout << "i,j: " << i << "," << j << endl;
//    cout << "n_wi,n_j: " << n_wi << "," << n_j << endl;
//    cout << "---------------------------------"<<endl;

    // Metropolis Sampling
    boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    if (rnum(rng) < R){ // Accept

        // Add to acceptance counter
        dkbh_accepts += 1;
        
        // Change regular kink on i to a worm head
        kinks_vector[kink_idx_i].tau = tau_h;
        kinks_vector[kink_idx_i].n = n_i;
        kinks_vector[kink_idx_i].src = i;
        kinks_vector[kink_idx_i].dest = i;
        kinks_vector[kink_idx_i].prev = prev_i;
        kinks_vector[kink_idx_i].next = next_i;
        
        // num_kinks-1,num_kinks-2 will be swapped. Modify links to these.
        kinks_vector[kinks_vector[num_kinks-1].next].prev = head_idx;
        kinks_vector[kinks_vector[num_kinks-1].prev].next = head_idx;
        
        // Head will be removed from j. Connect upper bound and kink.
        kinks_vector[next_j].prev = kink_idx_j;
        kinks_vector[kink_idx_j].next = next_j;
        
        // Modify last kinks vector (if necessary)
        if (kinks_vector[num_kinks-1].next==-1){
            last_kinks[kinks_vector[num_kinks-1].src] = head_idx;
        }
        
        swap(kinks_vector[head_idx],kinks_vector[num_kinks-1]);
        
        // Some important kinks could've been swapped. Correct if so.
        if (next_j==num_kinks-1){next_j=head_idx;}
        else if (prev_j==num_kinks-1){prev_j=head_idx;}
        else if (next_i==num_kinks-1){next_i=head_idx;}
        else if (prev_i==num_kinks-1){prev_i=head_idx;}
        else if (kink_idx_j==num_kinks-1){kink_idx_j=head_idx;}
        else if (kink_idx_i==num_kinks-1){kink_idx_i=head_idx;}
        else{;}
        
        // Tail could've been swapped. Correct if so.
        if (tail_idx==num_kinks-1){tail_idx=head_idx;}
        
        // num_kinks-1,num_kinks-2 will be swapped. Modify links to these.
        kinks_vector[kinks_vector[num_kinks-2].next].prev = kink_idx_j;
        kinks_vector[kinks_vector[num_kinks-2].prev].next = kink_idx_j;
        
        // Kink will be removed from j. Reconnect upper,lower bound.
        kinks_vector[next_j].prev = prev_j;
        kinks_vector[prev_j].next = next_j;
        
        // Modify last kinks vector (if necessary)
        if (kinks_vector[num_kinks-2].next==-1){
            last_kinks[kinks_vector[num_kinks-2].src] = kink_idx_j;
        }
 
        swap(kinks_vector[kink_idx_j],kinks_vector[num_kinks-2]);
        
        if (next_j==num_kinks-2){next_j=kink_idx_j;}
        else if (prev_j==num_kinks-2){prev_j=kink_idx_j;}
        else if (next_i==num_kinks-2){next_i=kink_idx_j;}
        else if (prev_i==num_kinks-2){prev_i=kink_idx_j;}
        else if (kink_idx_j==num_kinks-2){kink_idx_j=kink_idx_j;}
        else if (kink_idx_i==num_kinks-2){kink_idx_i=kink_idx_j;}
        else{;}
        
        // Tail could've been swapped
        if (tail_idx==num_kinks-2){tail_idx=kink_idx_j;}
        
        // Deactivate the deleted kinks
        kinks_vector[num_kinks-1].tau = -1;
        kinks_vector[num_kinks-1].n = -1;
        kinks_vector[num_kinks-1].src = -1;
        kinks_vector[num_kinks-1].dest = -1;
        kinks_vector[num_kinks-1].prev = -1;
        kinks_vector[num_kinks-1].next = -1;
        
        cout << endl << "first swap" << endl;
        for (int i=0; i<14; i++){
            cout << i << " " << kinks_vector[i] << endl;
        }
        
        // Deactivate the deleted kinks
        kinks_vector[num_kinks-2].tau = -1;
        kinks_vector[num_kinks-2].n = -1;
        kinks_vector[num_kinks-2].src = -1;
        kinks_vector[num_kinks-2].dest = -1;
        kinks_vector[num_kinks-2].prev = -1;
        kinks_vector[num_kinks-2].next = -1;
        
        cout << endl << "second swap" << endl;
        for (int i=0; i<14; i++){
            cout << i << " " << kinks_vector[i] << endl;
        }
            
        // Set new worm head index
        head_idx = kink_idx_i;

        // Update number of kinks tracker
        num_kinks -= 2;

        // Modify the last kinks vector if necessary
        if (next_j==-1){last_kinks[j]=prev_j;}
        
//        if (next_i==-1){last_kinks[i]=head_idx;}
        
        // Print out the indices of each sites last kink
        cout << "Last kink indices: ";
        for (int i=0; i<M ; i++){
            cout << last_kinks[i] << " ";
        }
        
        return;

    }
    else // Reject
        return;
    }

/*---------------------------- Estimators ------------------------------------*/

//double energy(vector<Kink> &kinks_vector, int num_kinks,
//              int M, int N, float U, float mu, float t, float beta){
//
//    float diagonal_energy=0;
//    vector<int> fock_state_half (M,0); // Store the Fock State at tau=beta/2
//    int current_kink_idx, n_at_half, next_kink_idx;
//    float tau_current, tau_next;
//    int n_i, N_half;
//
//    for (int src=0; src<M; src++){
//        current_kink_idx = src;
//        next_kink_idx = kinks_vector[current_kink_idx].next;
//        n_at_half = kinks_vector[current_kink_idx].n;
//        while (next_kink_idx!=-1){
//            tau_current = kinks_vector[current_kink_idx].tau;
//            tau_next = kinks_vector[next_kink_idx].tau;
//            if (tau_next<beta/2){
//                n_at_half = kinks_vector[next_kink_idx].n;
//            }
//            else{
//                if (tau_next-beta/2<beta/2-tau_current){
//                    n_at_half = kinks_vector[next_kink_idx].n;
//                }
//                else{
//                    n_at_half = kinks_vector[current_kink_idx].n;
//                }
//            }
//
//            current_kink_idx++;
//
//        }
//
//        fock_state_half[src] = n_at_half;
//    }
//
//    for (int src=0; src<M; src++){
//        n_i = fock_state_half[src];
//        diagonal_energy += (U/2)*n_i*(n_i-1) - mu*n_i;
//    }
//
//    // Sum the total particles at tau=beta/2
//    N_half=0;
//    for (int i=0; i<M; i++){
//        N_half += fock_state_half[i];
//    }
//
//    return N_half;
//
//    // This bad boy is actually returning the total particles at tau=beta/2
//}

/*----------------------------------------------------------------------------*/




/*----------------------------------------------------------------------------*/

// Main
int main(){
    
    // Time main function execution
    auto start = high_resolution_clock::now();
    
    // Create a uniform distribution with support: [0.0,1.0)
    boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    
    // Bose-Hubbard parameters
    int L = 4, D = 1, N = L;
    float t = 1.0, U = 10, mu = -3.5;
    vector<int> alpha;
    int M = pow(L,D); // total sites
    string boundary_condition = "pbc";
    
    // Simulation parameters
    float eta = 1.0, beta = 1.0;
    bool canonical = true;
    int sweeps=10000000;
    
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
    int insert_worm_attempts=0, insert_worm_accepts=0;
    int delete_worm_attempts=0, delete_worm_accepts=0;

    int insert_anti_attempts=0, insert_anti_accepts=0;
    int delete_anti_attempts=0, delete_anti_accepts=0;
    
    int insertZero_worm_attempts=0, insertZero_worm_accepts=0;
    int deleteZero_worm_attempts=0, deleteZero_worm_accepts=0;

    int insertZero_anti_attempts=0, insertZero_anti_accepts=0;
    int deleteZero_anti_attempts=0, deleteZero_anti_accepts=0;
    
    int insertBeta_worm_attempts=0, insertBeta_worm_accepts=0;
    int deleteBeta_worm_attempts=0, deleteBeta_worm_accepts=0;

    int insertBeta_anti_attempts=0, insertBeta_anti_accepts=0;
    int deleteBeta_anti_attempts=0, deleteBeta_anti_accepts=0;
    
    int advance_head_attempts=0, advance_head_accepts=0;
    int recede_head_attempts=0, recede_head_accepts=0;
    
    int advance_tail_attempts=0, advance_tail_accepts=0;
    int recede_tail_attempts=0, recede_tail_accepts=0;
    
    int ikbh_attempts=0, ikbh_accepts=0;
    int dkbh_attempts=0, dkbh_accepts=0;
    
    // Observables
    double N_sum=0;
    double diagonal_energy=0;
    
    // Non-observables
    int Z_ctr=0;
    
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

    cout << endl;
    // Print out the data structure
    cout << "Initial structure: " << endl;
    for (int i=0;i<8;i++){
        cout << kinks_vector[i] << endl;
    }
    cout << endl;
    
/*---------------------------- Monte Carlo -----------------------------------*/

    boost::random::uniform_int_distribution<> updates(0, 8);
    int label;
    
    sweeps *= (beta*M);
    sweeps = static_cast<int>(sweeps);
    for (int m=0; m < sweeps; m++){
        
//        if (m==57){break;}
        
        label = updates(rng);
        
        if (m > 1){
            cout << endl << "label: " << label << " m: " << m <<
            " num_kinks: " << num_kinks << " head_idx: " << head_idx <<
            " tail_idx: " << tail_idx << endl;
            for (int i=0; i<20; i++){
                cout << i << " " << kinks_vector[i] << endl;
            }
            // Print out the indices of each sites last kink
            cout << "Last kiNk indices before update: ";
            for (int i=0; i<M ; i++){
                cout << last_kinks[i] << " ";
            }
            cout << endl;
        }
                
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
        else{
            // lol
        }

//        // Print out the data structure
//        cout << endl << "label: " << label << endl;
//        for (int i=0;i<10;i++){
//            cout << kinks_vector[i] << endl;
//        }
//
//        // Print out the head and tail indices
//        cout << "head_idx: " << head_idx << endl;
//        cout << "tail_idx: " << tail_idx << endl;
//
//        // Print out the N_tracker
//        cout << "N_tracker: " << N_tracker << endl;
//
//        // Print out number of active kinks
//        cout << "num_kinks: " << num_kinks << endl;
//
//        // Print out the indices of each sites last kink
//        cout << "Last kink indices: ";
//        for (int i=0; i<M ; i++){
//            cout << last_kinks[i] << " ";
//        }
//        cout << endl;
        
        if (last_kinks[0]==last_kinks[1]
           || last_kinks[0]==last_kinks[2]
           || last_kinks[0]==last_kinks[3]
           || last_kinks[1]==last_kinks[2]
           || last_kinks[1]==last_kinks[3]
           || last_kinks[2]==last_kinks[3]){
            cout << "Things bad when : ";
            cout << "label " << label << " m " << m << endl;
            // Print out the indices of each sites last kink
            cout << "Last kiNk indices before update: ";
            for (int i=0; i<M ; i++){
                cout << last_kinks[i] << " ";
            }
            cout << endl;
            break;
        }
        
        // Measure N
        if (m%(static_cast<int>(M*beta))==0 && m>0.2*sweeps){
            if (head_idx==-1 and tail_idx==-1){
                N_sum += N_tracker;
                //diagonal_energy+=energy(kinks_vector,num_kinks,M,N,U,mu,t,beta);
                Z_ctr += 1;
            }
        }
    }

    // Print out the data structure
    cout << endl << "Final structure: " << endl;
    for (int i=0;i<18;i++){
        cout << i << "  " << kinks_vector[i] << endl;
    }
    cout << endl;

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
    
    // Print out the indices of each sites last kink
    cout << "Last kink indices: ";
    for (int i=0; i<M ; i++){
        cout << last_kinks[i] << " ";
    }
    cout << endl;

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
    
    auto end = high_resolution_clock::now();

    auto elapsed_time = duration_cast<nanoseconds>(end - start);
    float duration = elapsed_time.count() * 1e-9;
    
    cout << endl << "<n>: " << (N_sum/M)/Z_ctr << endl;
//    cout << endl << "<V/t>: " << (diagonal_energy/M)/Z_ctr << endl;
//    cout << endl << "<V/t>: " << (diagonal_energy)/Z_ctr + mu*N << endl;
    cout << endl << "Z_ctr: " << Z_ctr << endl;
    
    cout << endl << "Elapsed time: " << duration << " seconds" << endl;
    
//    // Print out adjacency matrix
//    cout << endl << "Adjacency Matrix: " << endl << endl;
//    for (int i=0; i<M; i++){
//        for (int j=0; j<total_nn; j++){
//            cout << adjacency_matrix[i][j] << " ";
//        }
//        cout << endl;
//    }
//    cout << endl;

    cout << endl;
    cout << "Total neighbors: " << total_nn << endl << endl;
    
    return 0;
}
