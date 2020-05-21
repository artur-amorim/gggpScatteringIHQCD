#include "IHQCD.h"
#include "FL.h"
#include "methods/interpolation/Poly_Interp.hpp"
#include "methods/vectorOperators.hpp"

class FLIzNIntegrand
{
    private:
        Poly_Interp<double> func1, func3;
        U1NNMode func2;
    public:
        FLIzNIntegrand(const Poly_Interp<double> &f1, const U1NNMode &f2, const Poly_Interp<double> &f3);
        double operator()(const double x);
};

FLIzNIntegrand::FLIzNIntegrand(const Poly_Interp<double> &f1, const U1NNMode &f2, const Poly_Interp<double> &f3):
                func1(f1), func2(f2), func3(f3) {}

double FLIzNIntegrand::operator()(const double x)
{
    double dfq = func2.dfQ(x);
    return func1.interp(x) * (dfq * dfq / func2.Q2()) * func3.interp(x);
}

double fFL(double * x, void * params)
{
    return ((FLIzNIntegrand *) params)->operator()(*x);
}


extern"C"
{
    void dqags_(double (*f)(double *, void *),  void * params,
                double * a, double * b, double * epsabs, double * epsrel, 
                double * result, double * abserr, int * neval,int * ier,
                int * limit, int * lenw, int * last, int * iwork, double * work);
}


FL::FL(const std::string file_path) : DeepInelasticScattering(file_path)
{
    std::cout << "Loaded FL." << std::endl;
}

double FL::IzN(const std::vector<double> &kin, const Reggeon &reg)
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
    FLIzNIntegrand integrand(f1, mode, f3);
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
    dqags_(fFL, params, &a, &b, &epsabs, &epsrel, &izn, &abserr, &neval, &ier, &limit, &lenw, &last, iwork, work);
    izn = std::pow(Q2, J) * izn;
    // Free workspace from memory
    delete[] iwork;
    delete[] work;
    return izn;
}