#include "IHQCD.h"
#include "F2.h"
#include "methods/interpolation/Poly_Interp.hpp"
#include "methods/vectorOperators.hpp"

class F2IzNIntegrand
{
    private:
        Poly_Interp<double> func1, func3;
        U1NNMode func2;
    public:
        F2IzNIntegrand(const Poly_Interp<double> &f1, const U1NNMode &f2, const Poly_Interp<double> &f3);
        double operator()(const double x);
};

F2IzNIntegrand::F2IzNIntegrand(const Poly_Interp<double> &f1, const U1NNMode &f2, const Poly_Interp<double> &f3):
                func1(f1), func2(f2), func3(f3) {}

double F2IzNIntegrand::operator()(const double x) {return func1.interp(x) * func2.factor(x) * func3.interp(x);}

double fF2(double * x, void * params)
{
    return ((F2IzNIntegrand *) params)->operator()(*x);
}

extern"C"
{
    void dqags_(double (*f)(double *, void *), void * params,
                double * a, double * b, double * epsabs, double * epsrel, 
                double * result, double * abserr, int * neval,int * ier,
                int * limit, int * lenw, int * last, int * iwork, double * work);
}

F2::F2(std::string file_path) : DeepInelasticScattering(file_path)
{
    std::cout << "Loaded F2." << std::endl;
}

double F2::IzN(const std::vector<double> &kin, const Reggeon &reg)
{
    // Get the kinematical values and J
    const double Q2 = kin[0];
    const double J = reg.getJ();
    std::vector<std::vector<double> > wf = reg.getWf();
    // Search for the correct mode
    std::vector<double> z = ihqcd().getZ();
    std::vector<double> fact1 = exp((1.5-J) * ihqcd().getAstring());
    Poly_Interp<double> f1(z, fact1, 4);
    U1NNMode mode = searchMode(Q2);
    // Interpolate the wavefunction
    Poly_Interp<double> f3(wf[0], wf[1], 4);
    // Define the integrand object
    F2IzNIntegrand integrand(f1, mode, f3);
    void * params = &integrand;
    // Compute the integral
    double izn = 0.0, abserr = 0.0;
    double a = z[0], b = z.back();
    // Absolute and relative tolerances desidered
    double epsabs = 1e-9, epsrel = 1e-9;
    if (Q2 < 2) epsabs = 1e-7;
    // Output variables 
    int neval = 0, ier;
    // Setup the workspace for the method
    int limit = 10000;
    int lenw = 100000;
    int last = 0;
    int * iwork = new int[limit];
    double * work = new double[lenw];
    // Evaluate the integral
    dqags_(fF2, params, &a, &b, &epsabs, &epsrel, &izn, &abserr, &neval, &ier, &limit, &lenw, &last, iwork, work);
    izn = std::pow(Q2, J) * izn;
    // Now we take into account the rescaling of psi_1 and psi_4
    if (reg.getIndex() == 1 || reg.getIndex() == 4) izn = -izn;
    // Free workspace from memory
    delete[] iwork;
    delete[] work;
    return izn;
}