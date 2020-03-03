#ifndef F2_H
#define F2_H

#include <vector>
#include "DeepInelasticScattering.h"

class F2_alphaQED : public DeepInelasticScattering {
    public :
        F2_alphaQED(const bool rrsslog = false, std::string file_path = "expdata/F2_data.txt");                                         // Constructor
        double IzN(const std::vector<double> &kin, const Reggeon &reg);
};

class F2IzNIntegrand: public IzNIntegrand
{
    public:
        F2IzNIntegrand(const Poly_Interp<double> &f1, const U1NNMode &f2, const Poly_Interp<double> &f3);
        double operator()(const double x);
};

#endif