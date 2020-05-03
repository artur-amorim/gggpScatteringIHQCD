#include <iostream>
#include <fstream>
#include <string>
#include "HQCDP.h"
#include "IHQCD.h"
#include "SoftPomeron.h"
#include "Reggeon.h"
#include "SigmaGammaProton.h"
#include "schrodinger/schrodinger.h"

using namespace std;

int main(int argc, char ** argv)
{
    double g1, g2;
    string data_path;
    if (argc < 3)
    {
        data_path = "expdata/GammaP/SigmaGammaProton.txt";
        g1 = 0; g2 = 0;
    }
    else
    {
        data_path = argv[1];
        g1 = stod(argv[2]); g2 = stod(argv[3]);
    }

    cout << "Using values:" << endl;
    cout << "g1: " << g1 << " g2: " << g2 << endl;
    
    // Create SigmaGammaProton object
    SigmaGammaProton sigma(data_path);
    vector<vector<double> > sigma_points = sigma.expKinematics();

    /// Setup SoftPomeron Kernel and GNs vector
    SoftPomeron soft;
    vector<double> GNs = {g1, g2};

    // Setup HQCDP object
    HQCDP hqcdp;
    hqcdp.addProcessObservable(sigma);
    hqcdp.addKernel(soft);
    hqcdp.setGNs(GNs);

    // Compute the spectrum to geth the reggeons
    chebSetN(1000);
    hqcdp.computeSpectrum();
    vector<Spectra> spectrum = hqcdp.getSpectrum();

    // Compute IzNBars
    cout << "Computing SigmaGammaP IzNBars" << endl;
    vector<kinStruct> SigmaGammaPIzNBars = sigma.getIzsBar(sigma_points, spectrum, GNs);

    // Compute IzNs
    cout << "Computing SigmaGammaP IzNs" << endl;
    vector<kinStruct> SigmaGammaPIzNs = sigma.getIzs(sigma_points, spectrum);

    // Compute sigma(gamma p -> hadrons)
    vector<double> sigma_pred = sigma.predict(SigmaGammaPIzNs, SigmaGammaPIzNBars, sigma_points, false);
    
    // Compute sigma(gamma p -> hadrons) chi2
    double sigma_chi2 = sigma.rss(SigmaGammaPIzNs, SigmaGammaPIzNBars, sigma_points);
    cout << "The sigma(gamma p -> hadrons) chi2 is " << sigma_chi2 / (sigma_points[0].size() -4 ) << endl;

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
    sigma_pred = sigma.predict(SigmaGammaPIzNs, SigmaGammaPIzNBars, kinPts, true);

    return 0;
}