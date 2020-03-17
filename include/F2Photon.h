#ifndef F2PHOTON_H
#define F2PHOTON_H

#include <vector>
#include "PhotonScattering.h"
#include "U1NNMode.h"

class F2Photon: public PhotonScattering
{
    private:
        std::vector<U1NNMode> modes;
        // Computes the necessary U1 nonnormalizable modes
        void computeU1NNModes();
        // Function that searches for the mode with virtuality Q2
        U1NNMode searchMode(const double Q2);
        void loadData(std::string file_path = "");
        void copy(const F2Photon & rhs);
    public :
        F2Photon(std::string data_path = "");
        F2Photon(const F2Photon &f2);
        std::vector<double> expVal();
        std::vector<double>  expErr();
        std::vector<std::vector<double> >  expKinematics();
        double IzN(const std::vector<double> &kin, const Reggeon &reg);
        double IzNBar(const std::vector<double> &kin, const Reggeon &reg, const std::vector<double> &gs);
        std::vector<kinStruct>  getIzs(const std::vector< std::vector<double> > &points, const std::vector<Spectra> &spec);
        std::vector<kinStruct>  getIzsBar(const std::vector< std::vector<double> >  &points, const std::vector<Spectra> &spec,
                                          const std::vector<double>  &gs);
        std::vector<double>  predict(const std::vector<kinStruct>  &Izs, const std::vector<kinStruct>  &IzsBar,
                                     const std::vector< std::vector<double> >  &points, const bool savePredictions);
        F2Photon& operator= (const F2Photon &procs);
        ~F2Photon();
};

#endif