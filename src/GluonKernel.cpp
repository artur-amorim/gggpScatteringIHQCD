#include "IHQCD.h"
#include "GluonKernel.h"
#include "methods/vectorOperators.hpp"

// GluonKernel constructor
GluonKernel::GluonKernel(const int nreg, const std::vector<double> &pars):
    Kernel("gluon", nreg, pars) {}

std::vector<double> GluonKernel::computePotential(const double J) const
{
    // Computes the potential values given J
    // Get the kernel parameters kernelPars = {invls, a, b ,c , d}
    std::vector<double> pars = this->KernelPars();
    const double invls = pars[0];
    // model = IHQCD
    std::vector<double> VSch = ihqcd().getU2() + (J-2) *  2.0 * invls * invls * ihqcd().getE2As() ;
    return VSch;
}