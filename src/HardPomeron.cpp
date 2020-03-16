#include "IHQCD.h"
#include "HardPomeron.h"
#include "methods/vectorOperators.hpp"

// HardPomeron constructor
HardPomeron::HardPomeron(const int n_reggeons, const std::vector<double> &pars):
    Kernel("HardPomeron", n_reggeons, pars) {}

std::vector<double> HardPomeron::computePotential(const double J) const
{
    // Computes the potential values given J
    // Get the kernel parameters kernelPars = {invls, a, b ,c , d}
    std::vector<double> pars = this->KernelPars();
    const double invls = pars[0];
    const double a = pars[1];
    const double b = pars[2];
    const double c = pars[3];
    const double d = pars[4];
    // model = IHQCD
    std::vector<double> VSch = ihqcd().getU2() + (J-2) * ( 2.0 * invls * invls * ihqcd().getE2As() * (1.0 + d / ihqcd().getl1_2()) + (J+2) * ihqcd().getE2A()  + (a * ihqcd().getaF() + b * ihqcd().getbF() + c * ihqcd().getcF())) ;
    return VSch;
}