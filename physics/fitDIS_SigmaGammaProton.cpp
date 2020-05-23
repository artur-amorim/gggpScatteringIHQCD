#include <iostream>
#include <vector>
#include "IHQCD.h"
#include "F2.h"
#include "FL.h"
#include "SigmaGammaProton.h"
#include "HardPomeron.h"
#include "HQCDP.h"
#include "schrodinger/schrodinger.h"

using namespace std;

int main(int argc, char ** argv)
{
    double g1, g2, g3, g4;
    string data_path_f2, data_path_fl, data_path_sigma;
    if (argc < 8)
    {
        data_path_f2 = "expdata/DIS/F2_data.txt";
        data_path_fl = "expdata/DIS/FL_data.txt";
        data_path_sigma = "expdata/GammaP/SigmaGammaProton.txt";
        g1 = 0.0; g2 = 0.0; g3 = 0.0; g4 = 0.0;
    }
    else
    {
        data_path_f2 = argv[1]; data_path_fl = argv[2]; data_path_sigma = argv[3];
        g1 = stod(argv[4]); g2 = stod(argv[5]); g3 = stod(argv[6]); g4 = stod(argv[7]);
    }
    cout << "Starting fit with values:" << endl;
    cout << "g1: " << g1 << " g2: " << g2 << " g3: " << g3 << " g4: " << g4 << endl;  
    // Setup sigma(gamma gamma -> hadrons) object
    F2 f2(data_path_f2);
    FL fl(data_path_fl);
    SigmaGammaProton sigma(data_path_sigma);

    // Setup SoftPomeron Kernel and GNs vector
    HardPomeron hard;
    vector<double> GNs = {g1, g2, g3, g4};

    // Setup HQCDP object
    HQCDP hqcdp;
    hqcdp.addProcessObservable(f2);
    hqcdp.addProcessObservable(fl);
    hqcdp.addProcessObservable(sigma);
    hqcdp.addKernel(hard);
    hqcdp.setGNs(GNs);

    cout << "Number of degrees of freedom: " << hqcdp.NumberOfDegreesOfFreedom() << endl;

    
    // Compute the spectrum to check Reggeon properties
    chebSetN(1000);
    hqcdp.computeSpectrum();

    vector<double> gns_guess = {g1, g2, g3, g4};
    double gns_delta = 1;
    hqcdp.fit(gns_guess, gns_delta);

    return 0;
}