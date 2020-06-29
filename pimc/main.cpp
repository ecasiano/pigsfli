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
#include<list>
#include<vector>
#include <boost/random.hpp>
using namespace std;
  
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

class Worldlines
{
    public:
    int worm_head_idx,worm_tail_idx;
    int num_kinks;
    
    // Member function declarations
    Worldlines (int,int,int); // constructor
//    bool insert_worm(),delete_worm();
//    bool insert_anti(),delete_anti();
//    bool advance_head(),recede_head();
//    bool advance_tail(),recede_tail();
//    bool insertZero_worm(),deleteZero_worm();
//    bool insertBeta_worm(),deleteBeta_worm();
//    bool insert_kink_before_head(),delete_kink_before_head();
//    bool insert_kink_after_head(),delete_kink_after_head();
//    bool insert_kink_before_tail(),delete_kink_before_tail();
//    bool insert_kink_after_tail(),delete_kink_after_tail();
    
    // Make "<<" a friend of the Worldline class
    friend ostream& operator<<(ostream& os, const Worldlines& dt);
    
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

Worldlines::Worldlines (int L,int D,int N){
    Kink kink (-1,-1,-1,-1,-1,-1);
    vector<Kink> kinks_array (100,kink);

    for (int i=0; i<L; i++){
        
        // Modify kink attributes
        kinks_array[i].tau = 0;
        kinks_array[i].n = 1;
        kinks_array[i].site = i;
        kinks_array[i].dir = 0;
        kinks_array[i].prev = -1;
        kinks_array[i].next = -1;
        }
}

// Overload "<<" operator
ostream& operator<<(ostream& os, const Kink& dt)
{
    os << '<' << dt.tau << ',' << dt.n << ',' << dt.site << ','
    << dt.dir << ',' << dt.prev << ',' << dt.next << '>' << endl;
    
    return os;
}

// Main
int main(){
    
    int L = 4, N = 4;
    int num_kinks = L;
    
    Kink kink (-1,-1,-1,-1,-1,-1);
    
    vector<Kink> kinks_array (100,kink);
                
    for (int i=0; i<100; i++){
        
        if (i<L){
        // Modify kink attributes
        kinks_array[i].tau = 0;
        kinks_array[i].n = 1;
        kinks_array[i].site = i;
        kinks_array[i].dir = 0;
        kinks_array[i].prev = -1;
        kinks_array[i].next = -1;
        }
        
        // print out the initialized kinks
        cout << kinks_array[i] << endl;
        
    }
    
    cout << "num_kinks: " << num_kinks << endl;
    
    // Test random number generation with Boost
    boost::random::mt19937 rng;
    boost::random::uniform_real_distribution<double> rnum(0.0, 1.0);
    for (int i = 0; i < 10; ++i) {
        std::cout << rnum(rng) << "\n";
    }
    
    return 0;
}
