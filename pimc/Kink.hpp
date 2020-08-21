//
//  Kink.hpp
//  pimc
//
//  Created by Emanuel Casiano-Diaz on 8/21/20.
//  Copyright © 2020 Emanuel Casiano-Díaz. All rights reserved.
//

#ifndef Kink_hpp
#define Kink_hpp

#include <stdio.h>
#include <iostream>
using namespace std;

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

#endif /* Kink_hpp */
