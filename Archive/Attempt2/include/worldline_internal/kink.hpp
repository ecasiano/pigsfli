#pragma once

struct Kink
{
    double tau = 0.0;
    int site = -1;
    int type = 0;

    // Partner on the OTHER site at the SAME time
    int partner = -1;

    // Partner on the SAME site at the OTHER time
    int pair_partner = -1;

    // Per-site time-ordered linked list
    int next = -1;
    int prev = -1;

    // Which bond this event belongs to
    int bond_id = -1;
};
