#pragma once

struct Kink
{
    double tau = 0.0; // imaginary time
    int site = -1;    // lattice site index
    int type = 0;     // +1 or -1 (arrival or departure)

    // Partner on the OTHER site at the SAME time (i <-> j)
    int partner = -1;

    // Partner on the SAME site at the OTHER time (kink <-> antikink)
    int pair_partner = -1;

    // Per-site time-ordered linked list
    int next = -1;
    int prev = -1;

    // Which bond this event belongs to
    int bond_id = -1;
};
