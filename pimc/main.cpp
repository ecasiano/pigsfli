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
using namespace std;
  
// A linked list node
class Kink
{
    public:
    // Attribute declarations
    float tau;
    int n,site,dir,prev,next;
    
    // Member function declarations (prototypes)
    Kink (float,int,int,int,int,int); // Kink constructor
    
    // Make << a friend of the Kink class
    friend ostream& operator<<(ostream& os, const Kink& dt);
};

// Member function definitions: Constructor
Kink::Kink (float a,int b,int c,int d,int e,int f){
    tau = a;
    n = b;
    site = c;
    dir = d;
    prev = e;
    next = f;
}

// Overload << operator
ostream& operator<<(ostream& os, const Kink& dt)
{
    os << '<' << dt.tau << ',' << dt.n << ',' << dt.site << ','
    << dt.dir << ',' << dt.prev << ',' << dt.next << '>' << endl;
    
    return os;
}

// Main
int main(){
    Kink kink (10,1,0,0,-1,-1);
    
    cout << kink.tau << endl;
    
    kink.tau = 26;
    
    cout << kink.tau << endl;
    
    cout << kink << endl;
    
    return 0;
}
