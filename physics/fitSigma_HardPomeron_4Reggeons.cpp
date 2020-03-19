#include <iostream>
#include <vector>
#include "IHQCD.h"
#include "SigmaGammaGamma.h"
#include "HardPomeron.h"
#include "HQCDP.h"
#include "schrodinger/schrodinger.h"

using namespace std;

int main(int argc, char ** argv)
{
    double g1, g2, g3, g4;
    string data_path;
    if (argc < 6)
    {
        g1 = 0.0; g2 = 0.0; g3 = 0.0; g4 = 0.0;
        data_path = "expdata/gammagammaScatteringL3Processed.txt";
    }
    else
    {
        data_path = argv[1];
        g1 = stod(argv[2]); g2 = stod(argv[3]); g3 = stod(argv[4]); g4 = stod(argv[5]);
    }
    cout << "Starting fit with values:" << endl;
    cout << "g1: " << g1 << " g2: " << g2 << " g3: " << g3 << " g4: " << g4 << endl;  
    // Setup sigma(gamma gamma -> hadrons) object
    SigmaGammaGamma sigma(data_path);

    // Setup SoftPomeron Kernel and GNs vector
    HardPomeron hard;
    vector<double> GNs = {g1, g2, g3, g4};

    // Setup HQCDP object
    HQCDP hqcdp;
    hqcdp.addProcessObservable(sigma);
    hqcdp.addKernel(hard);
    hqcdp.setGNs(GNs);

    cout << "Number of degrees of freedom: " << hqcdp.NumberOfDegreesOfFreedom() << endl;

    
    // Compute the spectrum to check Reggeon properties
    chebSetN(800);
    hqcdp.computeSpectrum();

    vector<double> gns_guess = {g1, g2, g3, g4};
    double gns_delta = 0.1;
    hqcdp.fit(gns_guess, gns_delta);

    return 0;
}