#ifndef HQCDP_H
#define HQCDP_H

#include <vector>
#include "Process.h"
#include "Kernel.h"
#include "Spectra.h"

struct NelderMead;

class HQCDP
{
private:
    std::vector<Process*> processes;                                        // Processes relevant for the computation
    std::vector<double> useTVals ;                                          // Values of t
    std::vector<Kernel *> kernels;                                          // Kernels that are going to be used
    std::vector<Spectra> spectrum;                                          // Object that contains the relevant reggeons
    std::vector<double> gns;                                                // Couplings parameters
    double operator()(const std::vector<double> &x);                        // Auxiliary function to fit data
    friend NelderMead;                                                      // Allow NelderMead struct to have access to operator()                                                           
public:
    HQCDP();                                                                                        // HQCDP constructor
    void setGNs(const std::vector<double>& gs);                                                     // Set gns
    std::vector<Kernel *> getKernels() const;                                                       // Get Kernels
    std::vector<Spectra> getSpectrum() const;                                                       // Get Spectrum
    std::vector<double> getUseTVals() const;                                                        // Get the tvals to be used
    void computeNeededTVals();                                                                      // Get all the ts to fit or predict the processes
    void addProcessObservable(Process &photon);                                                     // Add a process to be fitted or predicted
    void addKernel(Kernel &f);                                                                      // Add a kernel to be used in the fit or prediction
    void computeSpectrum(const std::vector< std::vector<double> > &kernelPars = {});                // Vector wish contains vectors with the parameters of each kernel
    int NumberOfDegreesOfFreedom();                                                                 // Computes the number of degrees of freedom of a fit
    double chi2();                                                                                  // Computes the chi2
    void fit(const std::vector<double> &xguess, const double delta);                                // Fit the GNs and GNtildes and the kernel parameters to the data 
};
#endif
