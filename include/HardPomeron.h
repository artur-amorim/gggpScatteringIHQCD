#ifndef HARD_POMERON_H
#define HARD_POMERON_H

#include <vector>
#include "Kernel.h"

class HardPomeron : public Kernel{
    public:
        // class constructor
        HardPomeron(const int n_reggeons = 4, const std::vector<double> &pars = {6.491, -4.567, 1.485, 0.653, -0.113});
        // Computes the potential of the gluon kernel
        std::vector<double> computePotential(const double J) const;
};
#endif