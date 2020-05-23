#include <iostream>
#include <vector>
#include "IHQCD.h"
#include "F2Photon.h"
#include "SigmaGammaGamma.h"
#include "HardPomeron.h"
#include "HQCDP.h"
#include "schrodinger/schrodinger.h"

using namespace std;

int main(int argc, char ** argv)
{
    double g1, g2, g3, g4;
    string data_path_f2, data_path_sigma;
    if (argc < 7)
    {
        g1 = 0.0; g2 = 0.0; g3 = 0.0; g4 = 0.0;
        data_path_f2 = "expdata/F2_photon/F2_photon_xmax_0.01.txt";
        data_path_sigma = "expdata/GammaGamma/gammagammaScatteringL3Processed.txt"; 
    }
    else
    {
        data_path_f2 = argv[1]; data_path_sigma = argv[2];
        g1 = stod(argv[3]); g2 = stod(argv[4]); g3 = stod(argv[5]); g4 = stod(argv[6]);
    }
    cout << "Using values:" << endl;
    cout << "g1: " << g1 << " g2: " << g2 <<  " g3: " << g3 << " g4: " << g4 << endl;  
    // Setup sigma(gamma gamma -> hadrons) object
    F2Photon photon(data_path_f2);
    SigmaGammaGamma sigma(data_path_sigma);

    // Setup HardPomeron Kernel and GNs vector
    HardPomeron hard;
    vector<double> GNs = {g1, g2, g3, g4};

    // Setup HQCDP object
    HQCDP hqcdp;
    hqcdp.addProcessObservable(photon);
    hqcdp.addProcessObservable(sigma);
    hqcdp.addKernel(hard);
    hqcdp.setGNs(GNs);

    cout << "Number of degrees of freedom: " << hqcdp.NumberOfDegreesOfFreedom() << endl;
    
    // Compute the spectrum to check Reggeon properties
    chebSetN(1000);
    hqcdp.computeSpectrum();

    vector<double> gns_guess = {g1, g2, g3, g4};
    double gns_delta = 0.1;
    hqcdp.fit(gns_guess, gns_delta);

    return 0;
}