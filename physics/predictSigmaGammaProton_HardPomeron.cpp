#include <iostream>
#include <fstream>
#include <string>
#include "HQCDP.h"
#include "IHQCD.h"
#include "HardPomeron.h"
#include "Reggeon.h"
#include "SigmaGammaProton.h"
#include "schrodinger/schrodinger.h"

using namespace std;

int main(int argc, char ** argv)
{
    double g1, g2, g3, g4;
    string data_path;
    if (argc < 6)
    {
        data_path = "expdata/GammaP/SigmaGammaProton.txt";
        g1 = 0.508067; g2 = 1.23958; g3 = 0.981595; g4 = -2.16804;
    }
    else
    {
        data_path = argv[1];
        g1 = stod(argv[2]); g2 = stod(argv[3]); g3 = stod(argv[4]); g4 = stod(argv[5]);
    }

    cout << "Using values:" << endl;
    cout << "g1: " << g1 << " g2: " << g2 << " g3: " << g3 << " g4: " << g4 << endl;
    
    // Create SigmaGammaProton object
    SigmaGammaProton sigma(data_path);
    vector<vector<double> > sigma_points = sigma.expKinematics();

    /// Setup HardPomeron Kernel and GNs vector
    HardPomeron hard;
    vector<double> GNs = {g1, g2, g3, g4};

    // Setup HQCDP object
    HQCDP hqcdp;
    hqcdp.addProcessObservable(sigma);
    hqcdp.addKernel(hard);
    hqcdp.setGNs(GNs);

    // Compute the spectrum to geth the reggeons
    chebSetN(1000);
    hqcdp.computeSpectrum();
    vector<Spectra> spectrum = hqcdp.getSpectrum();

    // Compute IzNBars
    cout << "Computing sigma(gamma p -> hadrons) IzNBars" << endl;
    vector<kinStruct> SigmaGammaPIzNBars = sigma.getIzsBar(sigma_points, spectrum, GNs);

    // Compute IzNs
    cout << "Computing sigma(gamma p -> hadrons) IzNs" << endl;
    vector<kinStruct> SigmaGammaPIzNs = sigma.getIzs(sigma_points, spectrum);
    
    // Compute sigma(gamma p -> hadrons) chi2
    const int nPoints = sigma_points[0].size();
    double sigma_chi2 = sigma.chi2(SigmaGammaPIzNs, SigmaGammaPIzNBars, sigma_points);
    cout << "Number of points: " << nPoints << endl;
    cout << "The sigma(gamma p -> hadrons) chi2 / points is " << sigma_chi2 / nPoints << endl;
    cout << "The sigma(gamma p -> hadrons) chi2 / Ndof is " << sigma_chi2 / (nPoints -4 ) << endl;

    vector<double> Ws;
    for(double W = 1.5; W < 300; W += 0.1) Ws.push_back(W);
    vector<double> WPlus(Ws.size(), 0.0), WMinus(Ws.size(), 0.0);
    vector<vector<double>> kinPts = {Ws, WPlus, WMinus};
    // Compute new IzNBars
    SigmaGammaPIzNBars = sigma.getIzsBar(kinPts, spectrum, GNs);
    // Compute new IzNs
    SigmaGammaPIzNs = sigma.getIzs(kinPts, spectrum);
    // Compute  predictions of sigma(gamma p -> hadrons)
    std::cout << "Predicting sigma(gamma p -> hadrons) for the given values of sqrt(s)" << std::endl;
    vector<double> sigma_pred = sigma.predict(SigmaGammaPIzNs, SigmaGammaPIzNBars, kinPts, true);

    return 0;
}