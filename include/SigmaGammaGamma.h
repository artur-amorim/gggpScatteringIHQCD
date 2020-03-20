#ifndef SIGGAMMAGAMMA_H
#define SIGGAMMAGAMMA_H

#include <vector>
#include "PhotonScattering.h"

class SigmaGammaGamma: public PhotonScattering
{
    private:
        void  loadData(std::string file_path = "");
    public :
        SigmaGammaGamma(std::string data_path = "");
        SigmaGammaGamma(const SigmaGammaGamma &sigma);
        std::vector<double> expVal();
        std::vector<double>  expErr();
        std::vector<std::vector<double> >  expKinematics();
        double IzN(const std::vector<double> &kin, const Reggeon &reg);
        double IzNBar(const std::vector<double> &kin, const Reggeon &reg, const std::vector<double> &gs);
        std::vector<kinStruct>  getIzs(const std::vector< std::vector<double> > &points, const std::vector<Spectra> &spec);
        std::vector<kinStruct>  getIzsBar(const std::vector< std::vector<double> >  &points, const std::vector<Spectra> &spec,
                                          const std::vector<double>  &gs);
        std::vector<double>  diffObsWeighted(const std::vector<kinStruct>  &Izs, const std::vector<kinStruct>  &IzsBar, const std::vector< std::vector<double> >  &points);
        std::vector<double>  predict(const std::vector<kinStruct>  &Izs, const std::vector<kinStruct>  &IzsBar,
                                     const std::vector< std::vector<double> >  &points, const bool savePredictions);
        SigmaGammaGamma& operator= (const SigmaGammaGamma &procs);
        ~SigmaGammaGamma();
};

#endif