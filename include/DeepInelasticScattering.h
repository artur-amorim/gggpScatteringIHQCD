#ifndef DIS_H
#define DIS_H

#include <vector>
#include "ProcessObservable.h"
#include "U1NNMode.h"

class DeepInelasticScattering: public ProcessObservable
{
    private:
        std::vector<U1NNMode> modes;                                                                        // Vector containing the necessary U1 nonnormalizable modes
        void loadData(std::string file_path = "");
        void copy(const DeepInelasticScattering & rhs);                                                     // copy function of DIS
    protected:
        U1NNMode searchMode(const double Q2);                                                               // Function that searches for the mode with virtuality Q2
    public:
        DeepInelasticScattering(const bool rrsslog = false, std::string file_path = "");                    // Class constructor
        DeepInelasticScattering(const DeepInelasticScattering &dis);                                        // Class copy constructor
        std::vector<double> expVal();                                                                       // Returns a vector with the experimental values of DIS Observable (F2 or FL)
        std::vector<double> expErr();                                                                       // Returns a vector with the experimental errors of DIS Observable (F2 or FL)
        std::vector< std::vector<double> > expKinematics();                                                 // Returns a vector of vectors with Q2, x as elements
        std::vector<double> getNeededTVals();                                                               // Returns a vector with the single element 0;
        void computeU1NNModes();                                                                            // Computes the necessary U1 nonnormalizable modes
        double IzNBar(std::vector<double> kin, const Reggeon &reg, const std::vector<double> &gs);
        std::vector<U1NNMode> getModes();                                                                   // Returns the relevant nonnormalizable modes
        std::vector<kinStruct> getIzs(const std::vector< std::vector<double> > &points, 
                                                  const std::vector<Spectra> &spec);                        // Gets all the Izs
        std::vector<kinStruct> getIzsBar(const std::vector< std::vector<double> > &points,
                                                     const std::vector<Spectra> &spec,
                                                     const std::vector < double > &gs,
                                                     const std::vector<double> &gtildes = {});              // Gets all the Izsbar
        std::vector<double> predictObs(const std::vector<kinStruct> &Izs,
                                       const std::vector<kinStruct> &IzsBar,
                                       const std::vector< std::vector<double> > &points, const bool savePredictions = false);    // Predicts the DIS observable
        DeepInelasticScattering& operator= (const DeepInelasticScattering &rhs);                            // Class assignment operator
        ~DeepInelasticScattering();                                                                         // Class destructor
};

class IzNIntegrand
{
    protected:
        Poly_Interp<double> func1, func3;
        U1NNMode func2;
    public:
        IzNIntegrand(const Poly_Interp<double> &f1, const U1NNMode &f2, const Poly_Interp<double> &f3);
        virtual double operator()(const double x);
};

double f(double * x, void * params);

#endif