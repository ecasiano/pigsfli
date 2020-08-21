//
//  Kink.cpp
//  pimc
//
//  Created by Emanuel Casiano-Diaz on 8/21/20.
//  Copyright © 2020 Emanuel Casiano-Díaz. All rights reserved.
//

#include "Kink.hpp"
#include <iostream>
using namespace std;

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
