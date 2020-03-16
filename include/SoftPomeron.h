#ifndef SOFT_POMERON_H
#define SOFT_POMERON_H

#include <vector>
#include "Kernel.h"

class SoftPomeron : public Kernel{
    public:
        // class constructor
        SoftPomeron(const std::vector<double> &pars = {1/0.178});
        // Computes the potential of the gluon kernel
        std::vector<double> computePotential(const double J) const;
};
#endif