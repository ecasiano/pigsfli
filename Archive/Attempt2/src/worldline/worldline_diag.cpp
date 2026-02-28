#include "worldline.hpp"

int Worldline::total_particle_number_from_segments() const
{
    int N = 0;

    for (int site = 0; site < params_.n_sites; ++site)
    {
        std::vector<Segment> segs;
        build_segments_for_site(site, segs);

        if (!segs.empty())
        {
            // Occupation on [0, first kink) is segs[0].occ
            N += segs[0].occ;
        }
        else
        {
            // No kinks on this site: occupation is just left boundary
            N += config_.left_boundary.occupations[site];
        }
    }

    return N;
}
