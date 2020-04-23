#include <iostream>
#include <vector>
#include "IHQCD.h"
#include "SigmaGammaGamma.h"
#include "SoftPomeron.h"
#include "HQCDP.h"
#include "schrodinger/schrodinger.h"

using namespace std;

int main(int argc, char ** argv)
{
    double g1, g2;
    string data_path;
    if (argc < 3)
    {
        g1 = 0.0; g2 = 0.0;
        data_path = "expdata/gammagammaScatteringL3Processed.txt";
    }
    else
    {
        data_path = argv[1];
        g1 = stod(argv[2]); g2 = stod(argv[3]);
    }
    cout << "Starting fit with values:" << endl;
    cout << "g1: " << g1 << " g2: " << g2 << endl;  
    // Setup sigma(gamma gamma -> hadrons) object
    SigmaGammaGamma sigma(data_path);

    // Setup SoftPomeron Kernel and GNs vector
    SoftPomeron soft;
    vector<double> GNs = {g1, g2};

    // Setup HQCDP object
    HQCDP hqcdp;
    hqcdp.addProcessObservable(sigma);
    hqcdp.addKernel(soft);
    hqcdp.setGNs(GNs);

    cout << "Number of degrees of freedom: " << hqcdp.NumberOfDegreesOfFreedom() << endl;

    
    // Compute the spectrum to check Reggeon properties
    chebSetN(800);
    hqcdp.computeSpectrum();

    vector<double> gns_guess = {g1, g2};
    double gns_delta = 10;
    hqcdp.fit(gns_guess, gns_delta);

    return 0;
}