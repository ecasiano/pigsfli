#ifndef EXPORT_WORLDLINE_HPP
#define EXPORT_WORLDLINE_HPP

#include "pimc.hpp"
#include <fstream>
#include <iostream>

// Export worldline to a text file for Python plotting
inline void export_worldline(const pimc::Worldline& wl,
                             int M,
                             double beta,
                             const std::string& filename)
{
    std::ofstream out(filename);

    if (!out) {
        std::cerr << "Error: could not open " << filename << " for writing\n";
        return;
    }

    for (int site = 0; site < M; site++) {
        int cur = wl.firstKink(site);

        while (cur != -1) {
            const auto& k = wl[cur];
            out << site << " " << k.tau << " " << k.n << "\n";
            cur = k.next;
        }
    }

    out.close();
}

// Print worldline to terminal for debugging
inline void print_worldline(const pimc::Worldline& wl, int M) {
    std::cout << "\nWorldline configuration:\n";
    for (int site = 0; site < M; site++) {
        std::cout << "Site " << site << ": ";
        int cur = wl.firstKink(site);

        while (cur != -1) {
            const auto& k = wl[cur];
            std::cout << "(tau=" << k.tau
                      << ", n=" << k.n
                      << ", partner=" << k.partner << ") ";
            cur = k.next;
        }
        std::cout << "\n";
    }
    std::cout << std::endl;
}

#endif
