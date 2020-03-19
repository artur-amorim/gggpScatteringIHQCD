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
    cout << "Using values:" << endl;
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

    // Test computeNeededTVals
    hqcdp.computeNeededTVals();
    vector<double> tvals = hqcdp.getUseTVals();
    cout << "Displaying the needed t values" << endl;
    for(int i = 0; i < tvals.size(); i++) cout << tvals[i] << " ";
    cout << endl;
    
    // Compute the spectrum to check Reggeon properties
    chebSetN(800);
    hqcdp.computeSpectrum();

    vector<Spectra> spectrum = hqcdp.getSpectrum();
    vector<Reggeon> reggeons = spectrum[0].getReggeons();

    cout << "Reggeon index\tReggeon J for t = " << spectrum[0].getT() << endl;
    for(int i = 0; i < reggeons.size(); i++)
    {
        cout << reggeons[i].getIndex() << '\t' << reggeons[i].getJ() << endl;
    }

    // Get chi2
    cout << "Computing chi2" << endl;
    double CHI2 = hqcdp.chi2();
    cout << "The value of chi2 for the gn values provided is " << CHI2 << endl;

    return 0;
}