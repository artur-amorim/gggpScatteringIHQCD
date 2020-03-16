#include "IHQCD.h"
#include "SoftPomeron.h"
#include "methods/vectorOperators.hpp"

// SoftKPomeron constructor
SoftPomeron::SoftPomeron(const std::vector<double> &pars):
    Kernel("SoftPomeron", 2, pars) {}

std::vector<double> SoftPomeron::computePotential(const double J) const
{
    // Computes the potential values given J
    // Get the kernel parameters kernelPars = {invls}
    std::vector<double> pars = this->KernelPars();
    const double invls = pars[0];
    // model = IHQCD
    std::vector<double> VSch = ihqcd().getU2() + 2.0 * invls * invls * (J - 2) * ihqcd().getE2As() ;
    return VSch;
}